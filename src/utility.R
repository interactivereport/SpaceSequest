loadUtilityPkg <- function(){
    require(BayesSpace)
    require(gsubfn)
}

msgError <- function(...){
    message("ERROR: ",...)
    q()
}

getConfig <- function(strF){
    if(!file.exists(strF))
        msgError("Config file dose NOT exist: ",strF)
    config <- yaml::read_yaml(strF)
    if(is.null(config$sample_meta) || !file.exists(config$sample_meta))
        msgError("Sample meta file is required in the config file")
    if(is.null(config$output))
        msgError("The output path is required in the config file")
    if(is.null(config$prj_name))
        msgError("The project name is required in the config file")
    suppressWarnings(clusterN <- as.numeric(config$clusterN))
    clusterN <- clusterN[!is.na(clusterN)]
    if(length(clusterN)==0)
        msgError("Numeric 'clusterN' is required!")
    config$clusterN <- clusterN
    meta <- data.frame(data.table::fread(config$sample_meta,header=T))
    return(list(config=config,meta=meta))
}

writeHDF <- function(D,strF){
    suppressPackageStartupMessages(require(rhdf5))
    if(attr(class(D),"package")=='Matrix') D <- as.matrix(D)
    D <- as.data.frame(D)
    colnames(D) <- gsub(" ","_",colnames(D))
    a <- suppressWarnings(file.remove(strF))
    h5createFile(strF)
    h5write(rownames(D),strF,"axis_row")
    h5write(colnames(D),strF,"axis_col")
    for(one in colnames(D)){
        if(class(D[[one]][1])%in%c("numeric","character")){
            h5write(D[[one]],strF,one)
        }else if(class(D[[one]][1])%in%c("factor")){
            h5createGroup(strF,one)
            h5write(levels(D[[one]]),strF,paste0(one,"/value"))
            h5write(as.integer(D[[one]]),strF,paste0(one,"/key"))
        }else if(class(D[[one]][1])%in%c("logical")){
            h5write(as.integer(D[[one]]),strF,one)
        }else{
            stop("writeHDF (R) unknown type: ",class(D[[one]][1]))
        }
    }
    H5close()
}
readHDF <- function(strF){
    suppressPackageStartupMessages(require(rhdf5))
    hInfo <- h5ls(strF)
    rName <- h5read(strF,"axis_row")
    cName <- h5read(strF,"axis_col")
    gName <- hInfo[hInfo$otype=='H5I_GROUP','name']
    df <- list()
    for(one in gName){
        val <- h5read(strF,paste0(one,"/value"))
        df[[one]] <- val[h5read(strF,paste0(one,"/key"))]
    }
    for(one in hInfo$name){
        if(one %in% c("axis_row",'axis_col','key','value',gName)) next
        df[[one]] <- h5read(strF,one)
        if(sum(!df[[one]] %in%c(0,1))==0)
            df[[one]]  <- df[[one]] ==1
    }
    H5close()
    return(data.frame(df,row.names=rName)[,cName])
}
writeFeather <- function(D,strF){
    suppressPackageStartupMessages(require(arrow))
    if(attr(class(D),"package")=='Matrix') D <- as.matrix(D)
    D <- data.frame(rowNames=rownames(D),D)
    colnames(D) <- gsub(" ","_",colnames(D))
    a <- arrow::write_feather(D,strF)
}
readFeather <- function(strF){
    suppressPackageStartupMessages(require(arrow))
    D <- arrow::read_feather(strF)
    if("rowNames" %in% colnames(D)){
        D <- data.frame(row.names=D[['rowNames']],D[-which(colnames(D)=='rowNames')])
    }
    return(D)
}

# reading visium as singleCellExperiment
sr_path_column="SpaceRanger_path"
batchKey="library_id"
sr_checkMeta <- function(config,meta){
    sName = config$sample_name
    if(is.null(sName))
        msgError("sample_name is required in the config file!")
    if(!sName %in% colnames(meta))
        msgError(sName," specified in config is not a column in the sample sheet!")
    if(!sr_path_column %in% colnames(meta))
        msgError(sr_path_column," is required column in the sample sheet!")
    
}
sr_read <- function(meta,sName,strRDS=NULL){
    suppressMessages(suppressWarnings(loadUtilityPkg()))
    message("Reading samples ...")
    if(!is.null(strRDS) && file.exists(strRDS)){
        message("\tThe raw file exists: ",strRDS,"\n\tIf a new read is required, please rename/remove the above file!")
        return()
    }
    all_slices <- list()
    for(i in 1:nrow(meta)){
        sID <- as.character(meta[i,sName])
        message("\t",sID)
        all_slices[[sID]] <- readVisium(meta[i,sr_path_column])
        colData(all_slices[[sID]]) <- cbind(colData(all_slices[[sID]]),
                                                      setNames(data.frame(sID),batchKey))
    }
    sce <- do.call(cbind,all_slices)
    colnames(sce) <- paste(colnames(sce),colData(sce)[[batchKey]],sep="_")
    colData(sce)[[batchKey]] <- as.factor(colData(sce)[[batchKey]])
    if(!is.null(strRDS))
        saveRDS(sce,file=strRDS)
    else
        return(sce)
    return()
}

