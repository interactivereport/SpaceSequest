loadRCTDPkg <- function(){
    require(spacexr)
    require(dplyr)
    source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/utility.R"))
    source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/readH5ad.R"))
}

RCTD_process <- function(config,strH5ad,strOut){
    sc_ref <- RCTD_getRef(config$rctd_scH5ad,config$rctd_annotation_obs,config$rctd_downsample,config$rctd_ignore_annotation)
    message("\tRCTD: processing")
    X <- getX(strH5ad)
    meta <- getobs(strH5ad)
    nUMI <- colSums(X)
    dir.create(file.path(dirname(strOut),"tmp"),showWarnings=F)
    RCTD <- c()
    for(sid in unique(meta[[batchKey]])){
        message("\t\t",sid)
        strOne <- file.path(dirname(strOut),"tmp",paste0(sid,".rds"))
        if(file.exists(strOne)){
            message("\t\tUsing previous results: ",strOne,"\n\t\t\tremove or rename the above file to rerun!")
            oneRCTD <- readRDS(strOne)
        }else{
            sel <- as.vector(meta[[batchKey]]==sid)
            oneD <- spacexr::SpatialRNA(meta[sel,c("array_row","array_col")],X[,sel],nUMI[sel])
            oneRCTD <- create.RCTD(oneD,sc_ref,max_cores=config$core)
            oneRCTD <- run.RCTD(oneRCTD,doublet_mode=config$rctd_doublet_mode)#'doublet'
            saveRDS(oneRCTD,strOne)
        }
        if(config$rctd_doublet_mode=="full"){
            RCTD <- rbind(RCTD,normalize_weights(oneRCTD@results$weights))
        }else if(config$rctd_doublet_mode=="doublet"){
            RCTD <- rbind(RCTD,get_doublet_weights(oneRCTD))
        }else if(config$rctd_doublet_mode=="multi"){
            RCTD <- rbind(RCTD,reshape2::dcast(
                dplyr::bind_rows(lapply(setNames(oneRCTD@results,colnames(oneRCTD@spatialRNA@counts)),function(one){
                    if(length(one$cell_type_list)==1){
                        data.frame(cell_type=one$cell_type_list,weights=one$sub_weights,row.names=NULL)
                    }else{
                        data.frame(cell_type=names(one$sub_weights),weights=one$sub_weights,row.names=NULL)
                    }
                }),.id='barcode'),
                barcode ~ cell_type, value.var = "weights", fill = 0) %>%
                    tibble::column_to_rownames(var = "barcode"))
        }else{
            stop("\tRCTD: Unknown double_mode: ",config$rctd_doublet_mode)
        }
    }
    colnames(RCTD) <- paste0("RCTD_",colnames(RCTD))
    writeFeather(RCTD,strOut)
}
RCTD_getRef <- function(strRef,anno_col,N=2000,NA_str=c("NA","NULL","unknown")){
    message("\tRCTD: build reference")
    X <- getX(strRef)
    meta <- getobs(strRef)
    selCell <- c()
    UMI <- colSums(X)
    message("\tRCTD: Downsampling for each annotation from total of ",length(UMI)," cells: top cells by UMI")
    for(one in unique(meta[,anno_col])){
        if(one%in%NA_str) next
        if(sum(meta[,anno_col]==one)<=N){
            selCell <- c(selCell,rownames(meta)[meta[,anno_col]==one])
        }else{
            selCell <- c(selCell,names(sort(UMI[meta[,anno_col]==one],decreasing=T))[1:N])
        }
    }
    message("\t",length(selCell)," cells selected")
    return(spacexr::Reference(X[,selCell], setNames(as.factor(meta[selCell,anno_col]),selCell), UMI[selCell]))
}

main <- function(){
    suppressWarnings(suppressPackageStartupMessages(loadRCTDPkg()))
    args = commandArgs(trailingOnly=TRUE)
    strConfig <- args[1]
    strH5ad <- args[2]
    strOut <- args[3]
    
    #list[config,meta] <- getConfig(strConfig)
    list2env(getConfig(strConfig), envir = environment())
    sr_checkMeta(config,meta)
    dir.create(dirname(strOut),showWarnings=F)
    
    sce <- RCTD_process(config,strH5ad,strOut)
    
    #strPDF <- gsub("rds$","pdf",strOut)
    #BayesSpace_plot(sce,strPDF)
    message("\n\n")
}

main()