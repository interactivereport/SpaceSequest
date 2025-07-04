## wrapper function for pipeline_class.R
## to have a better user experience, easier usage
PKGloading <- function(){
  require(rhdf5)
  require(Matrix)
  require(peakRAM)
  options(future.globals.maxSize=16*1024^3,stringsAsFactors=F)
}
#
checkFileExist <- function(strF,msg="file"){
  if(!file.exists(strF)){
    message("ERROR: ",msg," does not exist:\n",strF,"\n")
    q()
  }
  return(strF)
}
main <- function(strConfig){
  suppressMessages(suppressWarnings(PKGloading()))
  message("spDEG starting ...")
  config <- yaml::read_yaml(strConfig)
  strDEG <- checkFileExist(config$DEG_desp,"DEG description file")
  if(!is.null(config$prj_name)){
    prefix <- paste0(config$output,"/",config$prj_name)
    strH5adraw <- checkFileExist(paste0(file.path(config$output,"raw",config$prj_name),".h5ad"),"raw count h5ad file")
    strH5ad <- checkFileExist(paste0(file.path(config$output,config$prj_name),".h5ad"),"pipeline output h5ad contains cell annotation file")
  }else{
    prefix <- paste0(config$output,"/",config$DBname)
    strH5adraw <- checkFileExist(config$UMI,"raw count h5ad file")
    strH5ad <- checkFileExist(config$meta,"cell annotation file")
  }
  
  
  #sample,cluster,group,alt,ref,covars[comma separated],method[default NEBULA]
  #compInfo <- as.data.frame(data.table::fread(strDEG))
  suppressMessages(compInfo <- as.data.frame(data.table::fread(paste("grep -v '^#'",strDEG))))
  if(nrow(compInfo)==0){
    message("Empty comparison description file!\nDONE")
    q()
  }
  colnames(compInfo) <- sapply(strsplit(colnames(compInfo),"[[:punct:]]"),head,1)
  #meta <- getMeta(strH5ad)
  ## check compInfo against meta in checkInput function
  message("\tComparison description file checking ...")
  if("comparisonName"%in%colnames(compInfo) && length(unique(compInfo$comparisonName))!=nrow(compInfo))
    stop("Unique comparison names are required!")

  ## create spDEG folder
  message("\tcreating spDEG folder ...")
  spDEGpath <- paste0(prefix,"_spDEG/")
  dir.create(spDEGpath,showWarnings=F)

  ## enumerate all comparisons
  message("\tcreating spDEG tasks ...")
  cmds <- apply(compInfo,1,function(x){
    if(is.na(x["method"]) || nchar(x["method"])<1) x["method"] <- "NEBULA"
    if(grepl("nebula",x['method'],ignore.case=T) && (is.na(x["model"]) || nchar(x["model"])<1))
      x["model"] <- "HL"
    cluster_col <- x['cluster']
    if(grepl("cside",x['method'],ignore.case=T)) cluster_col=NULL
    coV <- NULL
    if(!is.na(x["covars"]) && nchar(x["covars"])>2){
      if(grepl("^\\{",x["covars"])) coV <- yaml::yaml.load(x["covars"])
      else coV <- unlist(strsplit(x["covars"],"\\+"),use.name=F)
    }
    if("comparisonName"%in%names(x)){
      message(x["comparisonName"])
      strOut <- file.path(spDEGpath,x["comparisonName"])
    }else{
      stop("comparisonName is required column in the DEG definition file")
    }
    if(dir.exists(strOut)) message("\tUsing existing: ",strOut)
    pipelineV <- ifelse(is.null(config$DEG_pipeline_version),'v1',config$DEG_pipeline_version)
    #system(paste("rm -rf",strOut))
    return(scRNAseq_DE(strH5adraw,strH5ad,strOut,x["method"],
                x["sample"],cluster_col,
                x["group"],x["ref"],x["alt"],
                x["model"],coV,
                NA_str=config$NAstring,

                R6_min.genes.per.cell = config$R6_min.genes.per.cell,
                R6_min.cells.per.gene = config$R6_min.cells.per.gene,
                R6_min.perc.cells.per.gene = config$R6_min.perc.cells.per.gene,
                R6_min.cells.per.gene.type = config$R6_min.cells.per.gene.type,
                R6_perc_threshold = config$R6_perc_threshold,
                R6_perc.filter.type = config$R6_perc.filter.type,
                R6_min.ave.pseudo.bulk.cpm = config$R6_min.ave.pseudo.bulk.cpm,
                R6_min.cells.per.subj = config$R6_min.cells.per.subj,
                
                cside_strRef=config$cside_strRef,
                cside_ct_col=config$cside_ct_col,
                cside_grp_col=config$cside_grp_col,
                cside_ct_N=config$cside_ct_N,
                cside_wt=config$cside_wt,
                cside_ct_count=config$cside_ct_count,
                cside_weight_threshold=config$cside_weight_threshold,
                cside_gene_threshold=config$cside_gene_threshold,
                cside_log_fc=config$cside_log_fc,
                cside_fdr=config$cside_fdr,
                cside_core = config$cside_core,
                
                #core=floor(as.numeric(gsub("G","",config$memory))/16),
                ver=pipelineV))
  })
  cmds <- unlist(cmds)#,use.names=F
  #cmds <- cmds[grepl("Astro",names(cmds))]
  #message(paste(paste(names(cmds),cmds,sep=": "),collapse="\n"))
  write(rjson::toJSON(cmds),paste0(prefix,"_spDEG.cmd.json"))
  #print(head(cmds))
  #writeLines(paste(cmds,collapse="\n"),paste0(prefix,"_spDEG.cmd"))
  cat("spDEG task creation completed")
  
}

scRNAseq_DE <- function(
    strCount,
    strMeta,
    output,
    method,
    column_sample,
    column_cluster,
    column_group=NULL,
    grp_ref=NULL,
    grp_alt=NULL,
    method_model=NULL,
    column_covars=NULL,
    NA_str=NULL,
    
    R6_min.genes.per.cell = 250,
    R6_min.cells.per.gene = 3,
    R6_min.perc.cells.per.gene = 0.1,
    R6_min.cells.per.gene.type = "or",
    R6_perc_threshold = 0.90,
    R6_perc.filter.type = "and",
    R6_min.ave.pseudo.bulk.cpm = 1,
    R6_min.cells.per.subj = 3,
    
    cside_strRef=NULL,
    cside_ct_col=NULL,
    cside_ct_N=500,
    cside_grp_col=NULL,
    cside_wt="c2l",
    cside_ct_count=0,
    cside_weight_threshold=0,
    cside_gene_threshold=5e-5,
    cside_log_fc=0.4,
    cside_fdr=0.05,
    cside_core=4,
    
    ver='v1',
    core=1,
    parallel=FALSE,
    addSRC=NULL
){
    env <- as.list(environment())
    checkInput(env)
    saveRDS(env,file=file.path(output,"env.rds"))
    if(grepl("^cside$",method,ignore.case=T)){
        cmd <- setNames(paste0("cd ",output,"; Rscript ",scRNAseq_DE_path,"/scRNAseq_DE.R ",
                               scRNAseq_DE_path,' CSIDE'),
                        paste(basename(output),"CSIDE",grp_alt,grp_ref,sep="_"))
    }else{
        meta <- getMeta(strMeta)
        cmd <- c()
        for(one in unique(meta[,column_cluster])){
            #"Rscript cmdPath/scRNAseq_DE.R cmdPath grpInterest" paste(jID,one,addSRC)
            cmd <- c(cmd,paste0("cd ",output,"; Rscript ",scRNAseq_DE_path,"/scRNAseq_DE.R ",
                                scRNAseq_DE_path,' "',one,'"'))
        }
        names(cmd) <- gsub(" ",".",paste(basename(output),unique(meta[,column_cluster]),grp_alt,grp_ref,sep="_"))
    }
    return(list(cmd))
}

getMeta <- function(strMeta){
  if(grepl("rds$",strMeta)){
    meta <- readRDS(strMeta)
  }else if(grepl("h5ad$",strMeta)){
    meta <- getobs(strMeta)
  }else{
    stop(paste("unknown meta format:",strMeta))
  }
  return(meta)
}
getUMI <- function(strUMI){
  if(grepl("rds$",strUMI)){
    UMI <- readRDS(strUMI)
  }else if(grepl("h5ad$",strUMI)){
    UMI <- getX(strUMI)
  }else{
    stop(paste("unknown UMI format:",strUMI))
  }
  return(UMI)
}

checkInput <- function(env){
    if(!file.exists(env$strCount) || !file.exists(env$strMeta)){
        stop("Either count RDS file or meta RDS file is missing!")
    }
    meta <- getMeta(env$strMeta)
    for(one in c(env$column_sample,env$column_cluster,env$column_group,env$column_covars)){
        if(!one%in%colnames(meta))
            stop(paste(one,"is not in the sample meta table!"))
    }
    if(!is.null(env$column_group)){
        if(is.null(env$grp_ref) || is.null(env$grp_alt))stop(paste("grp_ref and grp_alt are required for",env$column_group))
        for(one in c(env$grp_ref,env$grp_alt)){
            if(!one%in%unique(meta[,env$column_group]))
                stop(paste(one,"is not in the column",env$column_group,"from sample meta table!"))
        }
    }
    allMethods <- c("DESeq2","NEBULA","CSIDE")
    if(!env$method%in%allMethods)
        stop(paste0("method (",env$method,")is not supported!\nSelect from ",
                    paste(allMethods,collapse=", ")))
    if(env$method=="NEBULA"){
        if(!env$method_model%in%c("LN", "HL")) stop("method_model has to be LN or HL for NEBULA method!")
        if(env$method_model!="HL") warning("method_model is recommended to be HL for NEBULA method!")
    }
    else if(env$method=="glmmTMB"){
        if(!env$method_model%in% c("nbinom2", "nbinom1", "poisson", "nbinom2zi", "nbinom1zi"))
            stop("method_model has to be nbinom2, nbinom1, poisson, nbinom2zi or nbinom1zi for glmmTMB method!")
        if(env$method_model!="nbinom2") warning("method_model is recommended to be nbinom2 for NEBULA method!")
    }else if(env$method=="CSIDE"){
        if(is.null(env$cside_strRef) || !file.exists(env$cside_strRef))
            stop("Reference h5ad (cside_strRef) with raw UMI is required!")
        if(is.null(env$cside_ct_col))
            stop("Annotation column header (cside_ct_col) in reference cell meta is required!")
    }
    system(paste("mkdir -p",env$output))
}
checkCellMeta <- function(meta,env,cluster_interest){
  metaCol <- c(env$column_sample,env$column_group,env$column_covars)
  sel <- meta[,env$column_cluster]==cluster_interest
  if(!is.null(env$NA_str) && length(env$NA_str)>0){
    message(paste0("\tRemoving cell with '",paste(env$NA_str,collapse="' or '"),
                   " in columns of '",paste(metaCol,collapse="','"),"'"))
    sel <- sel & apply(meta[,metaCol],1,
                       function(x)return(length(intersect(x,env$NA_str))==0))
  }
  # check total cells
  if(sum(sel)<10){
    message("\tSkip: Too few cells (",sum(sel),")")
    return()
  }
  # check if any subject in both groups
  subjGroup <- ftable(meta[sel,c(env$column_sample,env$column_group)])
  print(subjGroup)
  if(sum(c(env$grp_ref,env$grp_alt)%in%data.frame(subjGroup)[[env$column_group]])!=2){
    message("Skip: missing at least one group!")
    return()
  }
  # check if at least 2 subjects in each group
  message("Checking minimal ",env$R6_min.cells.per.subj," cells per subject")
  subjGroupCell <- apply(subjGroup,2,function(x)return(sum(x>=env$R6_min.cells.per.subj)))
  if(sum(subjGroupCell[c(env$grp_ref,env$grp_alt)]>1)<2){
    message("Skip: One group contains less than 2 subjects!")
    return()
  }
  # check if pseudo bulk test, the cell level cov is not supported
  if(env$method %in% c("t_test","u_test","edgeR","limma","DESeq2")){
    if(length(unique(meta[,env$column_sample])) != nrow(unique(meta[,metaCol]))){
      message("Skip: cell level meta (e.g. pct_MT, library_size) is NOT supported in pseudo bulk comparison!")
      return()
    }
  }
  subjCell <- table(meta[sel,env$column_sample])
  sel <- sel & (meta[,env$column_sample]%in%names(subjCell)[subjCell>=env$R6_min.cells.per.subj])
  subjMeta <- meta[sel,c(env$column_sample,env$column_group,env$column_covars)]
  for(one in colnames(subjMeta)){
    if(sum(is.na(subjMeta[,one]))>0){
      message("Skip: NA found in ",one)
      return()
    }
  }
  subjMeta <- subjMeta[!duplicated(subjMeta[,env$column_sample]),]
  rownames(subjMeta) <- NULL
  print(subjMeta)
  return(as.vector(sel))
}

scRNAseq_DE_one_v1 <- function(env,cluster_interest,strSrc=NULL){
  if(!is.null(strSrc)) source(paste0(strSrc,"/scDEG.R"))
  strOut <- env$output
  filters <- list(rmGene=c("Mt-","MT-","mt-"),
                  min.cells.per.gene = env$R6_min.cells.per.gene,min.perc.cells.per.gene = env$R6_min.perc.cells.per.gene,
                  min.cells.per.gene.type = env$R6_min.cells.per.gene.type,
                  perc_threshold = env$R6_perc_threshold, perc.filter.type = env$R6_perc.filter.type,
                  lib_size_low = 200, lib_size_high = 20*10^6,
                  min.genes.per.cell = env$R6_min.genes.per.cell, min.cells.per.subj = env$R6_min.cells.per.subj,
                  min.ave.pseudo.bulk.cpm=env$R6_min.ave.pseudo.bulk.cpm)
  
  if(grepl("^cside$",env$method,ignore.case=T) && cluster_interest=="CSIDE"){
    env$column_cluster <- NA
    cluster_interest <- NULL
    strF <- file.path(strOut,paste0(env$grp_alt,".vs.",env$grp_ref,"_allDEG.rds"))
  }else{
    if(!is.null(env$column_group)){
      strF <- file.path(strOut,paste0(gsub("__","_",paste0(env$grp_alt,".vs.",env$grp_ref)),"__",
                                      gsub("__","_",paste(env$column_cluster,cluster_interest,sep=":")),"__",
                                      gsub("__","_",env$method),
                                      #"__cov:",paste(env)
                                      ".csv"))
      if(!is.null(env$column_covars)){
        strF <- gsub("\\.csv$",paste0("__cov:",paste(gsub("__","_",env$column_covars),collapse="+"),".csv"),strF)
      }
    }else{
      strF <- file.path(strOut,paste0(cluster_interest,".vs.Rest","__",gsub("_",".",env$column_cluster),".csv"))
      intrestGrp <- as.character(allMeta[,env$column_cluster])
      intrestGrp[intrestGrp!=cluster_interest] <- "Rest"
      allMeta <- cbind(allMeta,all="all",intrestGrp=intrestGrp)
      env$column_cluster <- "all"
      env$column_group <- "intrestGrp"
      env$grp_ref <- "Rest"
      env$grp_alt <- cluster_interest
      cluster_interest <- "all"
    }
  }
  if(file.exists(strF)){
    message("\tSkip: ",strF," exists!")
    return()
  }
  message("===== ",env$method,": ",cluster_interest," =====")
  sce <- scDEG$new(id_col=env$column_sample,
                   cluster_col=env$column_cluster,
                   grp_col=env$column_group,ctrl_value=env$grp_ref,alt_value=env$grp_alt,
                   strX=env$strCount,strMeta=env$strMeta,
                   NA_str=env$NA_str,
                   pipelinePath=strSrc)
  de <- sce$run(alg=env$method,fliter_list=filters,clusters=cluster_interest,
                covar=env$column_covars,method=env$method_model,
                csideParam=env[grepl("^cside",names(env))],
                prefix=fs::path_ext_remove(strF))
  if(grepl("^cside$",env$method,ignore.case=T)){
    saveRDS(de,strF)
    for(one in unique(de[[cside_ct_anno]])){
      strF <- file.path(strOut,paste0(gsub("__","_",paste0(env$grp_alt,".vs.",env$grp_ref)),"__",
                                      gsub("__","_",paste(env$cside_ct_col,one,sep=":")),"__",
                                      gsub("__","_",env$method),
                                      ".csv"))
      write.csv(de[de[[cside_ct_anno]]==one,],file=strF,row.names = FALSE)
    }
  }else{
      write.csv(de,file=strF,row.names = FALSE)
  }
  return()
}

scRNAseq_DE_dist <- function(env,cluster_interest,strSrc=NULL){
  if(env$ver=='v1'){
    message("=== spDEG pipeline v1 ===")
    return(scRNAseq_DE_one_v1(env,cluster_interest,strSrc))
  }else{
    stop("Unknown spDEG pipeline version: ",env$ver)
  }
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args)==1){
  a <- commandArgs()
  strPath <- gsub("^--file=","",grep("^--file",a,value=T)[1])
  scRNAseq_DE_path <<- dirname(normalizePath(strPath))
  source(paste0(scRNAseq_DE_path,"/readH5ad.R"))
  main(checkFileExist(args[1],"config file"))
}
if(length(args)>1){
  selGrp <- paste(args[-1],collapse=" ")
  #message("\n\n\n",args[-1],": ",selGrp,"\n\n\n")
  suppressMessages(suppressWarnings(PKGloading()))
  strAppPath <- dirname(gsub("--file=","",grep("file=",commandArgs(),value=T)))
  source(file.path(strAppPath,"readH5ad.R"))
  memUse <- peakRAM(Finfo <- scRNAseq_DE_dist(readRDS("env.rds"),
                                            selGrp,
                                            args[1]))
  memUse[1] <- selGrp
  print(memUse)
  data.table::fwrite(memUse,file="timeMem.csv",append=T)
  message("Successful!")
}

