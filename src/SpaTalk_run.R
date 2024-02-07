loadSpaTalkPkg <- function(){
    #require(tidyverse)
    require(dplyr)
    require(data.table)
    require(glue)
    require(SpaTalk)
    require(Seurat)
    source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/readH5ad.R"))
}
msg_previous <- function(strF){
    message("\t\tPrevious results is found and used: ",strF,
            "\n\t\t---> Please remove or rename the above file if a new run is wanted!)")
}
batchKey <- "library_id"

SpaTalk_sp <- function(X,L,species="Mouse",spot_max_cell){
    colnames(X) <- gsub("\\-","_",colnames(X))
    st_meta_qc <- L %>% tibble::rownames_to_column("sampleid_spot") %>% 
        dplyr::mutate(spot=gsub("\\-","_",sampleid_spot)) %>%
        dplyr::select(spot,x=array_col,y=array_row)
    st_count_qc_rev <- rev_gene(data = as.matrix(X),
                                data_type = "count",
                                species = species,
                                geneinfo = geneinfo)# geneinfo provided by SpaTalk package
    # create ST SpaTalk data
    obj <- createSpaTalk(st_data = st_count_qc_rev,
                         st_meta = st_meta_qc, 
                         species = species,
                         if_st_is_sc = F,
                         spot_max_cell = spot_max_cell)
    return(obj)
}
SpaTalk_sc <- function(scRef_count,scRef_ct,rm_ct,ct_rm_min=200,sel_ct=NULL,downsample_rate=1,downsample_min=3000){
    ct_count <- table(scRef_ct)
    rm_ct <- unique(c(rm_ct,names(ct_count)[ct_count<ct_rm_min]))
    if(is.null(sel_ct)) sel_ct <- unique(scRef_ct)
    sel <- !scRef_ct%in%rm_ct & scRef_ct%in%sel_ct
    if(downsample_rate<1 && downsample_rate>0){
        message("\t\tDownsampleing is enabled with rate: ",downsample_rate)
        ct_count <- table(scRef_ct[sel])
        ct_keep <- names(ct_count)[ct_count<downsample_min]
        sel[sample(seq_along(sel)[sel],round(sum(sel)*(1-downsample_rate)))] <- F
        sel <- sel | scRef_ct%in%ct_keep
    }
    ct_count <- table(scRef_ct[sel])
    message("\t\tCell type number: ",length(ct_count),"\n\t\t",paste(paste0(names(ct_count),":",ct_count),collapse="; "))
    
    return(list(umi=methods::as(scRef_count[,sel], "dgCMatrix"),ct=scRef_ct[sel]))
}
SpaTalk_c2l <- function(obj,C2L,sel_ct){
    if(is.null(C2L)) return(C2L)
    rownames(C2L) <- gsub("\\-","_",rownames(C2L))
    if(sum(colnames(C2L)%in%sel_ct)<3){
        stop("Error: Less than 3 common cell type names between sc reference annotation and C2L deconverlution")
    }
    return(as.matrix(C2L[rownames(C2L)%in%colnames(obj@data$rawdata),colnames(C2L)%in%sel_ct,drop=F]))
}
SpaTalk_plot <- function(obj,lri_celltype_df,out_dir,n_core){
    message("\tPloting lri_celltype")
    a <- mclapply(1:nrow(lri_celltype_df), function(x, obj_input=obj, lri_celltype_df_input=lri_celltype_df){
        tryCatch({
            #Plot the made-up spatial relationships that spatalk created
            ccdist_plot <- SpaTalk::plot_ccdist(object = obj_input, 
                                       celltype_sender = lri_celltype_df_input$sender[[x]], 
                                       celltype_receiver = lri_celltype_df_input$receiver[[x]],            
                                       size = 0.2,
                                       arrow_length = 0.02)
            ggsave(plot = ccdist_plot, 
                   filename = file.path(out_dir, glue("{lri_celltype_df_input$sender[[x]]}_{lri_celltype_df_input$receiver[[x]]}_ccdist.png") ), 
                   width=8,height=6,units="in",dpi=300)},
            error=function(e){
                message("\t\t",x,": Could not write out distance plot")
                print(e)
            })
        #If a top 20 plot exists, print it
        tryCatch({
            top20_lri <- plot_cci_lrpairs(object = obj_input,
                                          celltype_sender = lri_celltype_df_input$sender[[x]],
                                          celltype_receiver = lri_celltype_df_input$receiver[[x]],
                                          type = "sig")
            ggsave(plot = top20_lri, 
                   filename = file.path(out_dir, glue("{lri_celltype_df_input$sender[[x]]}_{lri_celltype_df_input$receiver[[x]]}_top20_cci_lrpairs.png") ), 
                   width=5,height=5,units="in",dpi=300)
        }, error=function(e){
            message("\t\t",x,": Could not write out top 20 plot.")
            print(e)
        })
    }, mc.cores = n_core)
}
SpaTalk_one <- function(X,L,scRef_umi,scRef_ct,strOut,
                        species="Mouse",spot_max_cell=30,
                        rm_ct=NULL,rm_ct_min_cell=100,sel_ct=NULL,down_ct_rate,down_ct_min,
                        C2L=NULL,ct_senders=NULL,ct_receivers=NULL,
                        n_perm=1000,min_per=0.2,n_mc_core=2,n_core=2,
                        n_neighbor=10,min_pairs=5,min_pairs_ratio=0,co_exp_ratio=0.1){#
    strF1 <- file.path(strOut,"dec_cci_all.rds")
    if(file.exists(strF1)){
        msg_previous(strF1)
        obj <- readRDS(strF1)
    }else{
        strF2 <- file.path(strOut,"dec_celltype.rds")
        if(file.exists(strF2)){
            msg_previous(strF2)
            obj <- readRDS(strF2)
        }else{
            # process one slide
            message("\tPreprocessing ...")
            obj <- SpaTalk_sp(X,L,species="Mouse",spot_max_cell)
            rm(X,L)
            scRef <- SpaTalk_sc(scRef_umi,scRef_ct,rm_ct,rm_ct_min_cell,sel_ct,down_ct_rate,down_ct_min)
            rm(scRef_umi,scRef_ct,rm_ct,sel_ct)
            gc()
            C2L <- SpaTalk_c2l(obj,C2L,unique(scRef$ct))
            message("\tStarting spatalk dec_celltype")
            st <- Sys.time()
            obj <- SpaTalk::dec_celltype(obj,
                                         sc_data = scRef$umi,
                                         sc_celltype = scRef$ct,
                                         iter_num = n_perm,
                                         use_n_cores = n_core,
                                         min_percent = min_per,
                                         method=ifelse(is.null(C2L),7,1),
                                         env=ifelse(is.null(C2L),sapply(strsplit(Sys.getenv('condaEnv_C2L'),' '),tail,1),'base'),
                                         anaconda_path=dirname(dirname(Sys.getenv('CONDA_EXE'))),
                                         dec_result = C2L)
            message("\tFinished spatalk dec_celltype with ",format(Sys.time()-st))
            saveRDS(obj,strF2)
        }
        message("\tStarting spatalk::dec_cci_all or dec_cci")
        st <- Sys.time() 
        data("lrpairs")
        data("pathways")
        obj <- find_lr_path(object = obj, lrpairs = lrpairs, pathways = pathways)
        
        if(F){
            # reduce n_core to avoid
            # slurmstepd: error: Detected 1 oom_kill event... Some of the step tasks have been OOM Killed.
            obj <- SpaTalk::dec_cci_all(obj,
                                        n_neighbor = n_neighbor,
                                        min_pairs = min_pairs,
                                        min_pairs_ratio = min_pairs_ratio,
                                        per_num = n_perm, 
                                        co_exp_ratio = co_exp_ratio,
                                        use_n_cores = n_mc_core)
        }else{
            message("run on individual pair")
            ct <- colnames(obj@coef)
            if(is.null(ct_senders)) ct_senders <- ct
            if(is.null(ct_receivers)) ct_receivers <- ct
            sender_receiver_df <- setNames(expand.grid(ct_senders,ct_receivers),c('sender','receiver'))
            sender_receiver_df <- sender_receiver_df[apply(sender_receiver_df,1,function(x)return(length(unique(x))))>1,]
            for(i in 1:nrow(sender_receiver_df)){
                #if(sender_receiver_df$sender[i]==sender_receiver_df$receiver[i]) next
                message("\t\tAnalyzing sender-receiver pair (",i,"/",nrow(sender_receiver_df),"): ",
                        sender_receiver_df$sender[i],"->",sender_receiver_df$receiver[i])
                tryCatch({
                    obj <- SpaTalk::dec_cci(obj,
                                            celltype_sender = sender_receiver_df$sender[i],
                                            celltype_receiver = sender_receiver_df$receiver[i],
                                            n_neighbor = n_neighbor, 
                                            min_pairs = min_pairs,
                                            min_pairs_ratio = min_pairs_ratio,
                                            per_num = n_perm,
                                            co_exp_ratio = co_exp_ratio,
                                            use_n_cores = n_core)}, 
                    error=function(e, obj_err=obj){
                        message("\t\t\tSomething went wrong for this pairing.")
                        print(e)
                        return()
                        #return(obj_err)
                    }
                )
            }
        }
        message("\tFinished spatalk::dec_cci_all with ",format(Sys.time()-st))
        saveRDS(obj,strF1)
    }
    lri <- obj@lrpair
    lri_celltype_df <- obj@lrpair %>%
        dplyr::select(sender=celltype_sender,receiver=celltype_receiver) %>% distinct()
    SpaTalk_plot(obj,lri_celltype_df,strOut,n_core)
    
    strF3 <- file.path(strOut,"get_lr_path.rds")
    if(file.exists(strF3)){
        msg_previous(strF3)
        lr_path_list <- readRDS(strF3)
    }else{
        obj@dist <- matrix()
        message("\tStarting performing pathway analysis")
        st <- Sys.time()
        lr_path_list <- mclapply(1:nrow(lri), function(x,obj_input=obj){
            #Isolate significant LRs
            lri_input <- obj_input@lrpair
            #Perform pathway analysis
            res <- get_lr_path(object = obj_input,
                               celltype_sender = lri_input$celltype_sender[[x]],
                               celltype_receiver = lri_input$celltype_receiver[[x]],
                               ligand = lri_input$ligand[[x]],
                               receptor = lri_input$receptor[[x]],
                               min_gene_num = 3 
            )
            return(res)
        },mc.cores=n_mc_core)
        saveRDS(lr_path_list,strF3)
        message("\tFinished pathway analysis with ",format(Sys.time()-st))
    }
    message("\tStarting extract pathway results")
    tf_df <- rbindlist( lapply(lr_path_list, "[[", 1) ) %>% distinct()
    lr_path_df <- rbindlist( lapply(lr_path_list, "[[", 2) )
    
    return(list(lri=lri,tf_df=tf_df,lr_path_df=lr_path_df))
}

SpaTalk_run <- function(config,strH5ad,strOut,strFinal,strC2L=""){
    config <- sapply(config,function(x){
        if(is.list(x) && length(x)==0)
            return(NULL)
        return(x)})
    if(file.exists(strFinal)){
        msg_previous(strFinal)
    }else{
        message("Loading ...")
        obs <- getobs(strH5ad)
        X <- getX(strH5ad)
        scRef_X <- getX(config$st_scH5ad)
        scRef_obs <- getobs(config$st_scH5ad)
        
        #remove cell tyoe with less cells
        ct_count <- table(scRef_obs[[config$tg_annotation_obs]])
        config$st_rm_ct <- unique(c(config$st_rm_ct,names(ct_count)[ct_count<config$st_rm_ct_min_cell]))
        c2l <- NULL
        if(nchar(strC2L)>3 && file.exists(strC2L)) c2l <- t(getX(strC2L))
        lri <- tf_df <- lr_path_df <- list()
        for(one in unique(obs[[batchKey]])){
            message("Processing ",one)
            selC <- rep(T,nrow(scRef_obs))
            if(!is.null(config$st_matchColumn) && config$st_matchColumn%in%colnames(scRef_obs)){
                matchKey <- unique(obs[obs[[batchKey]]==one,config$st_matchColumn])
                selC <- scRef_obs[[config$st_matchColumn]]%in% matchKey
                message("\tMatching scRef: ",paste(matchKey,collapse=", "))
            }
            spotIDs <- rownames(obs)[obs[[batchKey]]==one]
            strTmp <- file.path(strOut,one)
            dir.create(strTmp,showWarnings=F)
            res <- SpaTalk_one(X[,spotIDs],obs[spotIDs,],scRef_X[,selC],scRef_obs[[config$st_annotation_obs]][selC],strTmp,
                               config$st_species,config$st_spot_max_cell,
                               config$st_rm_ct,config$st_rm_ct_min_cell,config$st_sel_ct,
                               config$st_downsample_ct_rate,config$st_downsample_ct_min,
                               c2l[spotIDs,],config$st_senders,config$st_receivers,
                               config$st_iter_num,config$st_min_percent,config$st_mc_core,config$core,
                               config$st_n_neighbor,config$st_min_pairs,config$st_min_pairs_ratio,config$st_co_exp_ratio)
            lri[[one]] <- res$lri
            tf_df[[one]] <- res$tf_df
            lr_path_df[[one]] <- res$lr_path_df
        }
        
        spatalk_obj <- list(lri=dplyr::bind_rows(lri,.id=batchKey),
                            tf_df=dplyr::bind_rows(tf_df,.id=batchKey),
                            lr_path_df=dplyr::bind_rows(lr_path_df,.id=batchKey))
        saveRDS(spatalk_obj,strFinal)
    }
}

main <- function(){
    suppressMessages(suppressWarnings(loadSpaTalkPkg()))
    args = commandArgs(trailingOnly=TRUE)
    
    SpaTalk_run(sapply(yaml::read_yaml(args[1]),function(x){
        if(is.list(x) && length(x)==0)
            return(NULL)
        return(x)}),
                args[2],args[3],args[4],args[5])
    
}

main()


