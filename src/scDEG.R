require(R6)
require(cli)
require(BiocParallel)

allcell <- "spDEG_all_cells"
cside_ct_anno <- "CSIDE_CT_anno"

sparse_batch_apply <- function(X,MARGIN,FUN,shortLab="",...){
  step <- 500
  while(step>10){
    Xres <- c()
    ii <- 0
    #saveRDS(X,'test.rds')
    #print(class(X))
    Xlength <- ifelse(MARGIN==1,nrow(X),ncol(X))
    tryCatch({
      steps <- unique(c(seq(0,Xlength,step),Xlength))
      stepPair <- sapply(2:length(steps),function(i)return(list(c(steps[i-1]+1,steps[i]))))
      message(shortLab,": step ",step)
      Xres <- bplapply(stepPair,function(step,marg,oneFun){
        if(marg==1){
          A <- apply(X[step[1]:step[2],],MARGIN,FUN)#,...
        }else if(marg==2){
          A = apply(X[,step[1]:step[2]],MARGIN,FUN)#,...
        }else{
          stop("unknown MARGIN: ",MARGIN)
        }
        return(A)
      },marg=MARGIN,oneFun=FUN,#...,
      BPPARAM = MulticoreParam(workers=min(30,length(stepPair),max(1,parallelly::availableCores()-1)),
                               tasks=length(stepPair)))#,progressbar=T
      if(sum(!is.null(unlist(sapply(Xres,ncol))))==0){
        Xres <- unlist(Xres)
        #if(MARGIN==1) Xres <- Xres[rownames(X)]
        #else Xres <- Xres[colnames(X)]
      }else{
        Xres <- dplyr::bind_cols(Xres)
        #if(MARGIN==1) Xres <- Xres[,rownames(X)]
        #else Xres <- Xres[,colnames(X)]
      }
      step <- -1
    },
    error=function(cond){
      step <- floor(step/2)
      message("reduce step: ",step)
    })
  }
  #print(head(Xres))
  if(step==-1) return(Xres)
  stop("sparse_batch_apply")
}
# R6 class, save memory by keep one copy of the data matrix
# print out message for user to remove data outside R6 once the initialization process is finished
scDEG <- R6Class("scDEG",
                 public=list(
                   initialize=function(X=NULL,meta=NULL,
                                       id_col=NULL,cluster_col=NULL,grp_col=NULL,ctrl_value=NULL,alt_value=NULL,
                                       NA_str=NULL,
                                       strX=NULL,strMeta=NULL,
                                       pipelinePath=""){
                     message("Loading ...")
                     stopifnot(!is.null(id_col))
                     stopifnot(!is.null(cluster_col))
                     stopifnot(!is.null(grp_col))
                     stopifnot(!is.null(ctrl_value))
                     stopifnot(!is.null(alt_value))
                     stopifnot(!is.null(X) || (!is.null(strX) && file.exists(strX)))
                     stopifnot(!is.null(meta) || (!is.null(strMeta) && file.exists(strMeta)))
                     suppressMessages(suppressWarnings(private$loadPKG(pipelinePath)))
                     
                     message("Reading ...")
                     private$id_col <- id_col
                     private$grp_col <- grp_col
                     private$ctrl_value <- ctrl_value
                     private$alt_value <- alt_value
                     private$NA_str <- NA_str
                     if(is.null(meta)) private$getMeta(strMeta)
                     else private$meta <- meta
                     if(is.na(cluster_col)){## if all cells are included
                       cluster_col <- allcell
                       private$meta[,allcell]<- allcell
                     }
                     private$cluster_col <- cluster_col
                     private$check_meta()
                     
                     if(is.null(X)){
                       private$getUMI(strX)
                     }else{
                       private$X <- X
                       message("***** Please considering remove the original data matrix before 'run'! *****")
                     }
                     
                     if(sum(!rownames(private$meta)%in%colnames(private$X))>0){
                       stop("There are ",sum(!rownames(private$meta)%in%colnames(private$X))," cells defined in the meta but not in UMI!")
                     }
                     private$X <- private$X[,rownames(private$meta)]
                     
                     message("The scDEG initialization is completed successfully!")
                   },
                   run=function(alg,fliter_list,clusters=NULL,covar=NULL,prefix=NULL,method="HL",...){
                     if(T){
                       list2env(fliter_list,envir=environment())
                       #The following should be in the fitler_list
                       stopifnot(exists('rmGene'))
                       stopifnot(exists('min.cells.per.gene'))
                       stopifnot(exists('min.perc.cells.per.gene'))
                       stopifnot(exists('min.cells.per.gene.type'))
                       stopifnot(exists('perc_threshold'))
                       stopifnot(exists('perc.filter.type'))
                       stopifnot(exists('lib_size_low'))
                       stopifnot(exists('lib_size_high'))
                       stopifnot(exists('min.genes.per.cell'))
                       stopifnot(exists('min.cells.per.subj'))
                       stopifnot(exists('min.ave.pseudo.bulk.cpm'))
                     }
                     if(!is.null(covar) && (is.na(covar) || nchar(covar))) covar<-NULL
                     if(is.null(clusters)) clusters <- unique(private$meta[,private$cluster_col])
                     
                     allDEG <- NULL
                     for(one in clusters){
                       private$g_index <- rep(T,nrow(private$X))
                       private$c_index <- private$meta[,private$cluster_col]==one
                       private$scDEG_apply_filter_contrast(rmGene,
                                                   min.cells.per.gene,min.perc.cells.per.gene,min.cells.per.gene.type,
                                                   #set both number to be 0 to avoid any gene filtering
                                                   perc_threshold, perc.filter.type,
                                                   lib_size_low, lib_size_high,
                                                   min.genes.per.cell, min.cells.per.subj)
                       sInfo <- private$scDEG_check_model(covar)
                       if(!is.null(prefix)) data.table::fwrite(sInfo,paste0(prefix,"_",one,"_sampleInfo.csv"))
                       
                       if(alg %in%c('DESeq2')){
                         private$createPseudo()
                         private$scDEG_apply_filter_pseudoBulk(min.ave.pseudo.bulk.cpm)
                         message("pseudo bulk model: ",paste(colnames(private$sInfo),collapse="+"))
                       }
                       if(sum(private$c_index)<100) stop("Too few cells after filtering!")
                       oneDEG <- function(){
                         de <- switch (alg,
                           'NEBULA' = private$scDEG_Nebula(covar,method),
                           'DESeq2' = private$scDEG_DESeq2(...),
                           'CSIDE' = private$scDEG_CSDIE(prefix,method=method,...)
                         )
                         return(de)
                       }
                       message("Cluster: ",one)
                       de <- cbind(oneDEG(),cluster=one)
                       if(is.null(allDEG)) allDEG <- de
                       else allDEG <- rbind(allDEG,de)
                     }
                     return(allDEG)

                   }
                 ),
                 private=list(
                   X=NULL,meta=NULL,id_col=NULL,cluster_col=NULL,grp_col=NULL,ctrl_value=NULL,alt_value=NULL,NA_str=NULL,
                   strH5ad=NULL,cell_row='array_row',cell_col='array_col',
                   c_index=NULL,g_index=NULL,pseudoX=NULL,sInfo=NULL,
                   # preprocess -----
                   loadPKG=function(strPath){
                     stopifnot(require(data.table))
                     stopifnot(require(dplyr))
                     #stopifnot(require(glue))
                     #stopifnot(require(SingleCellExperiment))
                     #stopifnot(require(edgeR))
                     stopifnot(require(DESeq2))
                     stopifnot(require(apeglm))
                     #stopifnot(require(glmmTMB))
                     #stopifnot(require(BiocParallel))
                     #stopifnot(require(Seurat))
                     #stopifnot(require(MAST))
                     stopifnot(require(Matrix))
                     #stopifnot(require(slam))
                     #stopifnot(require(foreach))
                     #stopifnot(require(doMC))
                     #stopifnot(require(biomaRt))
                     #stopifnot(require(emmeans))
                     stopifnot(require(nebula))
                     #stopifnot(require(scater))
                     #stopifnot(require(peakRAM))
                     stopifnot(require(spacexr))
                     stopifnot(require(doParallel))
                     source(file.path(strPath,"readH5ad.R"))
                   },
                   getUMI=function(strUMI){
                     message('Getting UMI from ',strUMI)
                     if(grepl("rds$",strUMI)){
                       private$X <- readRDS(strUMI)
                     }else if(grepl("h5ad$",strUMI)){
                       private$X <- getX(strUMI)
                       private$strH5ad <- strUMI
                     }else{
                       stop(paste("unknown UMI format:",strUMI))
                     }
                     gc()
                   },
                   getMeta=function(strMeta){
                     if(grepl("rds$",strMeta)){
                       private$meta <- readRDS(strMeta)
                     }else if(grepl("h5ad$",strMeta)){
                       private$meta <- getobs(strMeta)
                     }else{
                       stop(paste("unknown meta format:",strMeta))
                     }
                   },
                   check_meta=function(){
                     if(!private$id_col %in% colnames(private$meta)) stop(paste(private$id_col,"NOT in meta"))
                     if(!private$cluster_col %in% colnames(private$meta)) stop(paste(private$cluster_col,"NOT in meta"))
                     if(!private$grp_col %in% colnames(private$meta)) stop(paste(private$grp_col,"NOT in meta"))
                     if(!private$ctrl_value %in% private$meta[,private$grp_col]) stop(paste(private$ctrl_value,"NOT in",private$grp_col))
                     if(!private$alt_value %in% private$meta[,private$grp_col]) stop(paste(private$alt_value,"NOT in",private$grp_col))
                     private$meta <- private$meta[order(private$meta[,private$id_col]),]
                   },
                   scDEG_apply_filter=function(rmGene=c("Mt-","MT-","mt-"), lib_size_low = 200, lib_size_high = 20*10^6,
                                               min.cells.per.gene = 50,min.perc.cells.per.gene = 0.01,
                                               min.genes.per.cell = 500){
                     # remove low expressed genes
                     if(length(rmGene)){
                       message("Filtering genes which start with ",paste(rmGene,collapse=", "))
                       private$g_index <- private$g_index & !grepl(paste(paste0("^",rmGene),collapse="|"),rownames(private$X))
                     }
                     if(min.perc.cells.per.gene>0 && min.perc.cells.per.gene<1){
                       message("\tCalculating minimal cell number of a gene based on  min.perc.cells.per.gene: ",min.perc.cells.per.gene)
                       min.cells.per.gene <- max(min.cells.per.gene,ceiling(min.perc.cells.per.gene*ncol(private$X)))
                     }
                     message("Filtering genes by minimal cell number: ",min.cells.per.gene)
                     private$g_index <- private$g_index & sparse_batch_apply(private$X,1,function(x)return(sum(x>0)),shortLab="CellN per gene")>min.cells.per.gene
                     message("--- Total of ",sum(private$g_index)," genes")
                     # remove low sequence depth cells
                     libSize <- colSums(private$X[as.vector(private$g_index),])
                     nGene <- sparse_batch_apply(private$X[as.vector(private$g_index),],2,function(x)return(sum(x>0,na.rm=T)),shortLab="GeneN per cell")
                     private$cellFiltering(libSize>lib_size_low & libSize<lib_size_high,
                                           paste0("Filtering cells by sequence depth: ",lib_size_low," ~ ",lib_size_high))
                     private$cellFiltering(nGene>min.genes.per.cell,
                                           paste0("Filtering cells by min genes: ",min.genes.per.cell))
                     message("--- Total of ",sum(private$c_index)," cells")
                   },
                   scDEG_apply_filter_contrast=function(rmGene=c("Mt-","MT-","mt-"),
                                                        min.cells.per.gene = 50,min.perc.cells.per.gene = 0.10,min.cells.per.gene.type = "and",
                                                        #set both number to be 0 to avoid any gene filtering
                                                        perc_threshold = 0.9, perc.filter.type = "and",
                                                        lib_size_low = 200, lib_size_high = 20*10^6,
                                                        min.genes.per.cell = 500, min.cells.per.subj = 5){
                     if(length(rmGene)>0){
                       message("Filtering genes which start with ",paste(rmGene,collapse=", "))
                       private$g_index <- private$g_index & !grepl(paste(paste0("^",rmGene),collapse="|"),rownames(private$X))
                     }
                     ctrl <- private$c_index & private$meta[,private$grp_col]==private$ctrl_value
                     alter <- private$c_index & private$meta[,private$grp_col]==private$alt_value
                     
                     ctrlG <- sparse_batch_apply(private$X[,as.vector(ctrl)],1,function(x)return(sum(x>0,na.rm=T)),shortLab="Ctrl cellN per gene")
                     alterG <- sparse_batch_apply(private$X[,as.vector(alter)],1,function(x)return(sum(x>0)),shortLab="Alt cellN per gene")
                     if(min.cells.per.gene>0){
                       message("Filtering genes by minimal cell number: ",min.cells.per.gene)
                       if(min.perc.cells.per.gene>0 && min.perc.cells.per.gene<1){
                         message("\tCalculating minimal cell number of a gene based on min.perc.cells.per.gene: ",min.perc.cells.per.gene)
                         min.cells.per.gene <- max(min.cells.per.gene,ceiling(min.perc.cells.per.gene*min(sum(ctrl),sum(alter))))
                       }
                       if(grepl('and',min.cells.per.gene.type,ignore.case=T)){
                         private$g_index <- private$g_index & (ctrlG>min.cells.per.gene & alterG>min.cells.per.gene)
                       }else{
                         private$g_index <- private$g_index & (ctrlG>min.cells.per.gene | alterG>min.cells.per.gene)
                       }
                       message("\t",sum(private$g_index)," genes")
                     }
                     
                     if(perc_threshold<1){
                       message("Filtering genes by maximal zero-expression percentile within a group: ",perc_threshold)
                       if(grepl('and',perc.filter.type,ignore.case=T)){
                         private$g_index <- private$g_index & (ctrlG/sum(ctrl) > (1-perc_threshold) & 
                                                                 alterG/sum(alter) > (1-perc_threshold))
                       }else{
                         private$g_index <- private$g_index & (ctrlG/sum(ctrl) > (1-perc_threshold) | 
                                                                 alterG/sum(alter) > (1-perc_threshold))
                       }
                     }
                     message("--- Total of ",sum(private$g_index)," genes")
                     stopifnot(sum(private$g_index)>100)
                     
                     # remove low sequence depth cells
                     message("remove low sequence depth cells")
                     libSize <- nGene <- rep(0,length(private$c_index))
                     libSize[ctrl|alter] <- colSums(private$X[as.vector(private$g_index),as.vector(ctrl|alter)])
                     nGene[ctrl|alter] <- sparse_batch_apply(private$X[as.vector(private$g_index),as.vector(ctrl|alter)],2,function(x)return(sum(x>0,na.rm=T)),shortLab="GeneN per cells")
                     private$cellFiltering(libSize>lib_size_low & libSize<lib_size_high,
                                   paste0("Filtering cells by sequence depth: ",lib_size_low," ~ ",lib_size_high))
                     private$cellFiltering(nGene>min.genes.per.cell,
                                   paste0("Filtering cells by min genes: ",min.genes.per.cell))
                     
                     # remove the samples with low cell number
                     cNum <- table(private$meta[private$c_index,private$id_col])
                     #if(sum(cNum<min.cells.per.subj)>0) message("removing: ",paste(names(cNum)[cNum<min.cells.per.subj],collapse=","))
                     private$cellFiltering(private$meta[,private$id_col]%in%names(cNum[cNum>=min.cells.per.subj]),
                                   paste0("Filtering samples with low cell number: ",min.cells.per.subj))
                     
                     message("--- Total of ",sum(private$c_index)," cells")
                   },
                   scDEG_check_model=function(covar){
                     if(is.list(covar)) covar <- names(covar)
                     if(sum(!covar%in%colnames(private$meta))>0) stop("Undefinied covar in meta!")
                     meta <- private$meta[private$c_index,c(private$id_col,private$cluster_col,private$grp_col,covar)] %>% 
                       group_by_at(private$id_col) %>% 
                       group_modify(function(x,k){
                         oneInfo <- data.frame(row.names=k,cellN=nrow(x))
                         endCol <- c()
                         for(one in colnames(x)){
                           if(is.numeric(x[[one]])){
                             oneInfo <- cbind(oneInfo,setNames(data.frame(mean(x[[one]])),paste0(one,"_mean")))
                             oneInfo <- cbind(oneInfo,setNames(data.frame(sd(x[[one]])),paste0(one,"_sd")))
                           }else{
                             nC <- table(x[[one]])
                             if(length(nC)==1){
                               oneInfo <- cbind(oneInfo,setNames(data.frame(names(nC)),one))
                             }else if(length(nC)>(nrow(x)/4)){
                               warning("\tSkip! There are more than a quater of unique values of ",one," in ",k)
                               oneInfo <- cbind(oneInfo,setNames(data.frame("Skip"),one))
                             }else{
                               oneInfo <- cbind(oneInfo,setNames(data.frame(paste(paste(names(nC),nC,sep=": "),collapse="; ")),one))
                               endCol <- c(endCol,one)
                             }
                           }
                         }
                         if(length(endCol)>0) oneInfo <- cbind(oneInfo[,!colnames(oneInfo)%in%endCol],oneInfo[,endCol])
                         return(oneInfo)
                       })
                     message("Please check the following sample meta information:")
                     print(data.frame(meta))
                     
                     private$sInfo <- as.data.frame(meta[,!(grepl('_sd$',colnames(meta)))])
                     colnames(private$sInfo) <- gsub("_mean$","",colnames(private$sInfo))
                     rownames(private$sInfo) <- private$sInfo[,private$id_col]
                     if(sum(!c(private$id_col,private$cluster_col,private$grp_col,covar)%in%colnames(private$sInfo))>0){
                       message("WARNING: not all defined column in sample information table!")
                     }else{
                       private$sInfo <- private$sInfo[,c(private$grp_col,covar),drop=F]
                       private$sInfo <- private$sInfo[,sapply(private$sInfo,function(x)return(sum(grepl("^Skip$",x)|(grepl(":",x)&grepl(";",x)))==0)),drop=F]
                     }
                     return(meta)
                   },
                   scDEG_apply_filter_pseudoBulk=function(min.ave.pseudo.bulk.cpm = 1){
                     sTotal <- colSums(private$pseudoX)
                     g_avg_cpm <- apply(private$pseudoX,1,function(x)return(mean(x/sTotal*1e6)))
                     message("Filtering genes by minimal average pseudo CPM: ",min.ave.pseudo.bulk.cpm)
                     private$pseudoX <- private$pseudoX[g_avg_cpm>min.ave.pseudo.bulk.cpm,]
                     message("--- Total of ",nrow(private$pseudoX)," genes")
                   },
                   # DEG by NEBULA pipeline -------
                   scDEG_Nebula=function(covar,method){
                     meta <-  private$meta[private$c_index,c(private$grp_col,covar),drop=F]
                     message(private$grp_col)
                     print(head(meta))
                     meta[,private$grp_col] <- factor(meta[,private$grp_col],levels=c(private$ctrl_value,private$alt_value))
                     df = model.matrix(as.formula(paste0("~",paste(colnames(meta),collapse="+"))),data=meta)
                     re = nebula(private$X[as.vector(private$g_index),as.vector(private$c_index)],
                                 private$meta[private$c_index,private$id_col],
                                 pred=df,offset=colSums(private$X[as.vector(private$g_index),as.vector(private$c_index)]),method=method)
                     final_table <- data.frame(ID=re$summary[,"gene"],re$summary[,grep(private$grp_col,colnames(re$summary))])
                     colnames(final_table) <- c("ID","logFC","se","Pvalue")
                     final_table$log2FC <- log2(exp(final_table$logFC))
                     final_table$FDR <- p.adjust(final_table[,"Pvalue"],method="fdr")
                     final_table$method <- "NEBULA"
                     final_table$algorithm <- re$algorithm
                     first4 <- c("ID","log2FC","Pvalue","FDR")
                     final_table <- final_table[,c(first4,colnames(final_table)[!colnames(final_table)%in%first4])]
                     return(final_table)
                   },
                   # DEG by DESeq2 ---------
                   scDEG_DESeq2=function(...){
                     private$X <- NULL
                     gc()
                     register(MulticoreParam(max(1,parallelly::availableCores()-1)))
                     meta <- private$sInfo
                     meta[,private$grp_col] <- factor(meta[,private$grp_col],levels=c(private$ctrl_value,private$alt_value))
                     #X <- private$pseudoX
                     #save(X,meta,file="test.RData")
                     dds <- DESeqDataSetFromMatrix(countData=private$pseudoX,
                                                   colData=meta,
                                                   design=as.formula(paste0("~",paste(colnames(meta),collapse="+"))))
                     dds <- DESeq(dds,parallel=T)
                     res = lfcShrink(dds, coef = paste(private$grp_col, private$alt_value, "vs", private$ctrl_value, sep = "_"),
                                     format="DataFrame",type = "apeglm", parallel = T)
                     final_table <- data.frame(row.names=NULL,
                                               ID=rownames(res),
                                               log2FC=res$log2FoldChange,
                                               Pvalue=res$pvalue,
                                               FDR=res$padj,
                                               se=res$lfcSE,
                                               method="DESeq2",
                                               algorithm='lfcShrink')
                     return(final_table)
                   },
                   # DEG by C-SIDE -----
                   scDEG_CSDIE=function(prefix=NULL,csideParam=list(),cside_clusters=NULL,method=NA,...){
                     list2env(csideParam,envir=environment())
                     stopifnot(exists('cside_strRef'))
                     stopifnot(exists('cside_ct_col'))
                     stopifnot(exists('cside_ct_N'))
                     stopifnot(exists('cside_grp_col'))
                     stopifnot(exists('cside_wt'))
                     stopifnot(exists('cside_ct_count'))
                     stopifnot(exists('cside_weight_threshold'))
                     stopifnot(exists('cside_gene_threshold'))
                     stopifnot(exists('cside_log_fc'))
                     stopifnot(exists('cside_fdr'))
                     stopifnot(exists('cside_core'))
                     message("CSIDE model: ",method)
                     cside_meta_regression<- grepl('regression',method,ignore.case=T)
                     cside_merge <- grepl('merge',method,ignore.case=T)
                     strRCTD <- paste0(prefix,"_cside.rds")
                     if(file.exists(strRCTD)){
                       message("Using previous results: ",strRCTD)
                       message("\tIf a new run is wanted, please rename/remove the above file!")
                       RCTD <- readRDS(strRCTD)
                     }else{
                       #core <- max(1,parallelly::availableCores()-5)
                       message("Getting reference")
                       sc_ref <- private$get_cside_ref(cside_strRef,cside_ct_col,cside_ct_N)
                       w <- private$get_cside_weight(cside_wt)
                       u_cell_types <- unique(sc_ref@cell_types)
                       if(!is.null(w) && sum(!u_cell_types%in%colnames(w))>0){
                           message("The following cell types in reference did NOT include in weight matrix:")
                           message("\t",paste(u_cell_types[!u_cell_types%in%colnames(w)],collapse=", "))
                       }
                       if(is.null(cside_clusters)) cside_clusters <- u_cell_types
                       if(cside_merge){
                           message("Creating RCTD ...")
                           suppressMessages({
                               RCTD <- create.RCTD(spatialRNA=private$get_cside_merged_puck(),
                                                   reference=sc_ref,
                                                   max_cores = cside_core,
                                                   keep_reference = T,
                                                   CELL_MIN_INSTANCE = 0)#Don't drop celltypes on the basis of minimum number
                           })
                           message("Updating weight matrix ...")
                           if(is.null(w)){
                               RCTD <- run.RCTD(RCTD,doublet_mode='doublet')
                           }else{
                               RCTD <- private$update_cside_weight(w,RCTD,cside_merge)
                           }
                           puck_exvar <- private$get_cside_exvar(RCTD,cside_merge)
                           suppressMessages({
                               RCTD <- run.CSIDE.single(RCTD,
                                                        explanatory.variable=puck_exvar,
                                                        cell_types=cside_clusters,
                                                        weight_threshold = cside_weight_threshold,
                                                        gene_threshold = cside_gene_threshold,
                                                        doublet_mode=F,
                                                        cell_type_threshold=cside_ct_count,
                                                        log_fc_thresh = cside_log_fc,
                                                        fdr = cside_fdr)
                           })
                       }else{
                           if(cside_meta_regression) private$check_cside_model()
                           grp_ids <- private$get_cside_repgroup(cside_grp_col)#cside_grp_col
                           message("Getting puck replicates")
                           puck_list <- private$get_cside_puck()
                           message("Creating RCTD ...")
                           suppressMessages({
                               RCTD_reps <- create.RCTD.replicates(spatialRNA.replicates = puck_list,
                                                                   reference = sc_ref,
                                                                   replicate_names = names(puck_list),
                                                                   group_ids = grp_ids[names(puck_list)],
                                                                   max_cores = cside_core,
                                                                   keep_reference = T,
                                                                   CELL_MIN_INSTANCE = 0) #Don't drop celltypes on the basis of minimum number
                           })
                           message("Updating weight matrix ...")
                           if(is.null(w)){
                               RCTD_reps <- run.RCTD.replicates(RCTD_reps, doublet_mode='doublet')
                           }else{
                               RCTD_reps <- private$update_cside_weight(w,RCTD_reps,cside_merge)
                           }
                           puck_exvar <- private$get_cside_exvar(RCTD_reps,cside_merge,cside_meta_regression)
                           ## run cside
                           message("Running cside ...")
                           suppressMessages({
                               RCTD_reps <- run.CSIDE.replicates(RCTD_reps,
                                                                 cell_types=cside_clusters,
                                                                 explanatory.variable.replicates=puck_exvar,
                                                                 weight_threshold = cside_weight_threshold,
                                                                 log_fc_thresh = cside_log_fc,
                                                                 doublet_mode = F, #Ran with full weights, not doublet
                                                                 gene_threshold = cside_gene_threshold,
                                                                 population_de = F, # Get DE across all cell types, run this later
                                                                 cell_type_threshold = cside_ct_count,
                                                                 fdr = cside_fdr)
                           })
                           ## Run meta regression (covariates) to compare groups
                           # this is not tested, not including
                           message("Inferencing DE results ...")
                           RCTD <- private$get_cside_inference(RCTD_reps,cside_log_fc,cside_fdr,cside_meta_regression)
                       }
                       saveRDS(RCTD,strRCTD)
                     }
                     return(private$extract_cside_deg(RCTD,cside_merge))
                   },
                   check_cside_model=function(){
                       sel_grp <- table(private$meta[private$c_index,c(private$grp_col,private$id_col)])
                       if(sum(apply(sel_grp,2,function(x)return(sum(x>0)))!=1)>0) stop(paste("CSIDE model is incorrect: 'group' column is not associated with 'sample' column"))
                   },
                   get_cside_ref=function(strH5ad,anno_col,N){
                     message("Creating reference from ",strH5ad)
                     message("\tPlease make sure doublets are annotated in column ",anno_col," as defined in 'NA_string'")
                     X <- getX(strH5ad)
                     meta <- getobs(strH5ad)
                     selCell <- c()
                     UMI <- colSums(X)
                     message("\tDownsampling for each annotation from total of ",length(UMI)," cells: top cells by UMI")
                     for(one in unique(meta[,anno_col])){
                       if(one%in%private$NA_str) next
                       if(sum(meta[,anno_col]==one)<=N){
                         selCell <- c(selCell,rownames(meta)[meta[,anno_col]==one])
                       }else{
                         selCell <- c(selCell,names(sort(UMI[meta[,anno_col]==one],decreasing=T))[1:N])
                       }
                     }
                     message("\t",length(selCell)," cells selected")
                     return(spacexr::Reference(X[,selCell], setNames(as.factor(meta[selCell,anno_col]),selCell), UMI[selCell]))
                   },
                   get_cside_puck=function(){
                     if(!private$cell_row %in%colnames(private$meta) || !private$cell_col%in%colnames(private$meta))
                       stop(private$cell_row," and ",private$cell_col," are required in cell meta information")
                     nUMI <- colSums(private$X)
                     sIDs <- unique(private$meta[,private$id_col][private$c_index])
                     message("Creating RCTD puck list")
                     puck_list <- lapply(setNames(sIDs,sIDs),function(sid){
                       message("\t",sid)
                       sel <- private$c_index & private$meta[,private$id_col]==sid
                       return(spacexr::SpatialRNA(private$meta[sel,c(private$cell_row,private$cell_col)],
                                                  private$X[,as.vector(sel)],
                                                  nUMI[sel]))
                     })
                     return(puck_list)
                   },
                   get_cside_merged_puck=function(){
                       if(!private$cell_row %in%colnames(private$meta) || !private$cell_col%in%colnames(private$meta))
                           stop(private$cell_row," and ",private$cell_col," are required in cell meta information")
                       nUMI <- colSums(private$X)
                       sIDs <- unique(private$meta[,private$id_col][private$c_index])
                       message("Creating merged RCTD puck")
                       step_row <- max(private$meta$array_row)+1
                       step_col <- max(private$meta$array_col)+1
                       colN <- ceiling(sqrt(length(sIDs)))
                       coords <- private$meta[,c(private$cell_row,private$cell_col)]
                       for(i in 1:length(sIDs)){
                           sel <- private$meta[,private$id_col]==sIDs[i]
                           coords[sel,] <- coords[sel,] + rep(c(step_row*(ceiling(i/colN)-1),step_col*((i-1)%%colN)),each=sum(sel))
                       }
                       return(spacexr::SpatialRNA(coords[private$c_index,],
                                                  private$X[,as.vector(private$c_index)],
                                                  nUMI[private$c_index]))
                   },
                   get_cside_repgroup=function(grp_col){
                     if(is.null(grp_col)) return(NULL)
                     if(!grp_col%in%colnames(private$meta)){
                       message(grp_col," is not in the cell meta column!")
                       return(NULL)
                     }
                     grp <- unique(private$meta[private$c_index,c(private$id_col,grp_col)])
                     return(setNames(as.numeric(factor(grp[,grp_col])),
                                     grp[,private$id_col]))
                     
                   },
                   get_cside_weight=function(wt_col){
                     if(is.null(wt_col)) return(NULL)
                     message("get annotation weight matrix: ",wt_col)
                     if(wt_col=='c2l'){
                       w <- private$meta[private$c_index,grepl("^c2l_",colnames(private$meta))]
                       colnames(w) <- gsub("^c2l_","",colnames(w))
                     }else if(wt_col=='tg'){
                       if(is.null(private$strH5ad) || !file.exists(private$strH5ad))
                         stop("Missing input h5ad (UMI) for tangram weights")
                       if(!"tangram_ct_count" %in% getobsmKey(strWeight))
                         stop("Missing 'tangram_ct_count' in obsm in ",private$strH5ad)
                       obsm <- getobsm(private$strH5ad,"tangram_ct_count")
                       w <- obsm[,!colnames(obsm)%in%c('x','y','centroids','cell_n')]
                       w <- cbind(w,tg_unknow=obsm$cell_n-apply(w,1,sum))
                     }else if(wt_col=='RCTD'){
                         w <- private$meta[private$c_index,grepl("^RCTD_",colnames(private$meta))]
                         colnames(w) <- gsub("^RCTD_","",colnames(w))
                     }else{
                       message("Unknown weight matrix setting (wt_col): ",wt_col)
                       message("\tUsing RCTD from spacexr")
                       return(NULL)
                     }
                     w <- w[,!colnames(w)%in%private$NA_str]
                     w <- t(apply(w,1,function(x)return(x/sum(x))))
                     return(w)
                   },
                   update_cside_weight=function(w,RCTD,cside_merge=F){
                     if(!is.null(w)){
                         if(cside_merge){
                             barcodes <- RCTD@spatialRNA@counts@Dimnames[[2]]
                             RCTD <- import_weights(RCTD, w[barcodes,])
                             RCTD@config$doublet_mode <- "full"
                             RCTD@config$RCTDmode <- "full"
                         }else{
                             for (i in 1:length(RCTD@RCTD.reps)) {
                                 # Subset each replicate barcodes
                                 barcodes <- RCTD@RCTD.reps[[i]]@spatialRNA@counts@Dimnames[[2]]
                                 RCTD@RCTD.reps[[i]] <- import_weights(RCTD@RCTD.reps[[i]], w[barcodes,])
                                 # Need to set the doublet_mode as "full" so that the code will use this weights in slot: myRCTD@results$weights
                                 RCTD@RCTD.reps[[i]]@config$doublet_mode <- "full"
                                 RCTD@RCTD.reps[[i]]@config$RCTDmode <- "full"
                                 #message(i,": ",RCTD@RCTD.reps[[i]]@config$doublet_mode)
                             }
                         }
                     }
                     #sapply(RCTD@RCTD.reps,function(one)print(one@config$doublet_mode))
                     return(RCTD)
                   },
                   get_cside_exvar=function(RCTD,cside_merge,cside_meta_regression){
                       if(cside_merge){
                           barcodes <- RCTD@spatialRNA@counts@Dimnames[[2]]
                           puck_exvar <- setNames(as.integer(private$meta[barcodes,private$grp_col] == private$alt_value),
                                                  barcodes)
                           if(length(unique(puck_exvar))<2) stop('Only one comparison group exists!')
                       }else{
                           puck_exvar <- list()
                           nUMI <- colSums(private$X)
                           for(i in 1:length(RCTD@RCTD.reps)){
                               barcodes <- RCTD@RCTD.reps[[i]]@spatialRNA@counts@Dimnames[[2]]
                               if(cside_meta_regression){
                                   exvar <- as.integer(nUMI[barcodes]>mean(nUMI[barcodes]))
                               }else{
                                   exvar <- as.integer(private$meta[barcodes,private$grp_col] == private$alt_value)
                               }
                               if(length(unique(exvar))<2) stop('Only one comparison group exists in replicate: ',
                                                              names(RCTD@group_ids)[i],
                                                              "\nCould consider either 'regression' or 'merge' in CSIDE model")
                               names(exvar) <- barcodes
                               puck_exvar[[i]] <-  exvar
                           }
                       }
                     return(puck_exvar)
                   },
                   get_cside_covar=function(covar){
                     if(is.null(covar)) return(NULL)
                     meta <- private$meta[private$c_index,] %>% 
                       dplyr::select(any_of(c(private$id_col,names(covar)))) %>% 
                       dplyr::distinct()
                     if(nrow(meta) != length(unique(meta[,private$id_col]))){
                       message("Some covar columns are NOT slice-specific:")
                       print(meta)
                       return(NULL)
                     }
                     rownames(meta) <- meta[,private$id_col]
                     meta <- meta[,!colnames(meta)==private$id_col,drop=F]
                     for(one in names(covar)){
                       if(is.character(meta[[one]]) || is.factor(meta[[one]])){
                         if(is.null(covar[[one]])) covar[[one]] <- unique(meta[,one])
                         meta[[one]] <- factor(meta[one],levels=covar[[one]])
                       }
                     }
                     return(model.matrix(as.formula(paste("~",paste(colnames(meta),collapse="+"))),meta))
                   },
                   get_cside_inference=function(RCTD_reps,cside_log_fc,cside_fdr,cside_meta_regression){
                       if(cside_meta_regression){
                           sel_grp <- table(private$meta[private$c_index,c(private$grp_col,private$id_col)])[,names(RCTD_reps@group_ids)]
                           meta.design <- data.frame(alt=rep(0,ncol(sel_grp)))
                           meta.design[sel_grp[private$alt_value,]>0,'alt'] <- 1
                           return(CSIDE.population.inference(RCTD_reps,log_fc_thresh=cside_log_fc,fdr=cside_fdr,
                                                             meta = TRUE,
                                                             meta.design.matrix = meta.design,
                                                             meta.test_var = 'alt'))
                       }else{
                           return(CSIDE.population.inference(RCTD_reps,log_fc_thresh=cside_log_fc,fdr=cside_fdr))
                       }
                   },
                   extract_cside_deg=function(RCTD,cside_merge){
                     first4 <- c('gene',"log2FC","p","q_val")
                     if(cside_merge){
                         res_df <- purrr::map(RCTD@de_results$all_gene_list, function(x){
                                x %>% 
                                 tibble::rownames_to_column("gene") %>%
                                 dplyr::mutate(q_val=p.adjust(p_val,method='fdr'))
                             }) %>%
                             dplyr::bind_rows(.id=cside_ct_anno) %>% 
                             dplyr::rename(p=p_val) %>%
                             dplyr::mutate(log10p=-log10(p),
                                           log2FC= log2(exp(log_fc)),
                                           !!cside_ct_anno := gsub(" ","_",.data[[cside_ct_anno]]))
                     }else{
                         res_df <- purrr::map(RCTD@population_de_results, function(x) x %>% tibble::rownames_to_column("gene")) %>%
                             dplyr::bind_rows(.id=cside_ct_anno) %>% 
                             dplyr::mutate(log10p=-log10(p),
                                    log2FC= log2(exp(log_fc_est)),
                                    !!cside_ct_anno := gsub(" ","_",.data[[cside_ct_anno]]))
                     }
                     res_df <- res_df[,c(first4,colnames(res_df)[!colnames(res_df)%in%first4])]
                     return(res_df)
                   },
                   #create pseudo -------
                   createPseudo=function(){
                     #private$pseudoX <- NULL
                     message("Creating pseudo bulk")
                     private$pseudoX <- bplapply(rownames(private$sInfo),function(one){
                       return(setNames(data.frame(rowSums(private$X[as.vector(private$g_index),as.vector(private$c_index&private$meta[,private$id_col]==one)])),
                                       one))
                     },BPPARAM = MulticoreParam(workers=min(30,nrow(private$sInfo),max(1,parallelly::availableCores()-1)),
                                              tasks=nrow(private$sInfo),progressbar=T))
                     private$pseudoX <- as.matrix(dplyr::bind_cols(private$pseudoX))
                     #for(one in rownames(private$sInfo)) private$pseudoX <- cbind(private$pseudoX,rowSums(private$X[private$g_index,private$c_index&private$meta[,private$id_col]==one]))
                     #colnames(private$pseudoX) <- rownames(private$sInfo)
                     mode(private$pseudoX) <- 'integer'
                   },
                   #cell filtering -------
                   filter_subset=function(subsetList){
                     if(is.null(subsetList)) return()
                     message("C-SIDE covars filtering")
                     for(one in names(subsetList)){
                       if(is.null(subsetList[[one]])) next
                       private$c_index <- private$c_index & private$meta[,one]%in%subsetList[[one]]
                       message("\t",one,": ",sum(private$c_index)," cells")
                     }
                   },
                   cellFiltering=function(sel,msg){
                     message(msg)
                     private$c_index <- private$c_index & sel
                     message("\t",sum(private$c_index)," cells")
                   }
                 ))


testFun <- function(){
  strH5ad <- "/mnt/depts/dept04/compbio/users/zouyang/tools/scRNAsequest/example/testFull/TST11837_oyoung_batch_raw_obsAdd.h5ad"
  strH5ad <- "/mnt/depts/dept04/compbio/projects/htvc/data/single_cell/human/AD_sc_atlas_2023rerun/scDEG/07Fujita_ExNeuron/scAD_2023rerun_raw_obsAdd_labelFilter.h5ad"
  #sce <- scDEG$new(id_col='library_id',cluster_col='predicted.A.celltype1',grp_col='Genotype',ctrl_value='WT',alt_value='DMSXL',
  #                 strX=strH5ad,strMeta=strH5ad,pipelinePath="/mnt/depts/dept04/compbio/users/zouyang/tools/scRNAsequest/src")
  sce <- scDEG$new(id_col='subjectID',cluster_col='celltype_coarse',grp_col='AD_status',ctrl_value='Control',alt_value='AD',
                   strX=strH5ad,strMeta=strH5ad,pipelinePath="/mnt/depts/dept04/compbio/users/zouyang/tools/scRNAsequest/src")
  filters <- list(rmGene=c("Mt-","MT-","mt-"),
                  min.cells.per.gene = 50,min.perc.cells.per.gene = 0.10,min.cells.per.gene.type = "and",
                  perc_threshold = 0.9, perc.filter.type = "and",
                  lib_size_low = 200, lib_size_high = 20*10^6,
                  min.genes.per.cell = 500, min.cells.per.subj = 5,
                  min.ave.pseudo.bulk.cpm=1)
  return(sce$run('DESeq2',filters,'ExNeuron',c('Sex','age_death'),method='HL'))#NEBULA
}
#fTable <- testFun()
#source("src/scDEG.R");testFun()