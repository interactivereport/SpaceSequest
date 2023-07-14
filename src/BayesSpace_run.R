loadBayesPkg <- function(){
    require(scater)
    require(BayesSpace)
    require(harmony)
    require(viridis)
    require(mclust)
    require(ggplot2)
    require(BiocParallel)
    source(paste0(dirname(gsub("--file=","",grep("file=",commandArgs(),value=T))),"/utility.R"))
}
BayesSpace_clusterKey <- "BayesSpace_clusterN"
BayesSpace_process <- function(config,strRaw,strOut){
    strProcess <- gsub("\\.rds","_preprocess.rds",strRaw)
    if(file.exists(strProcess)){
        message("Prcessed file exists: ",strProcess,"\n\t If a new process is wanted, please remove/rename the above file")
        sce <- readRDS(strProcess)
    }else{
        message("Preprocessing ...")
        sce <- readRDS(strRaw)
        # QC
        message("\tQC")
        if(!is.null(config$rmGeneStart)){
            isRM <- grepl(paste(paste0("^",config$rmGeneStart),collapse="|"),rowData(sce)$gene_name)
            if(sum(isRM)>0){
                message("\t\tThe following ",sum(isRM)," genes will be removed from BayesSpace analysis:")
                message("\t\t",paste(rowData(sce)$gene_name[isRM],collapse=", "))
                sce <- sce[!isRM,]
            }else{
                message("No rmGene is found!")
            }
        }
        sizeFactors = scuttle::librarySizeFactors(sce)
        sce = sce[,sizeFactors!=0] #remove spots with no expression
        
        # preprocess update the layout to avoid overlapping
        message("\tUpdating layouts")
        sce$row.orig <- sce$row
        sce$col.orig <- sce$col
        sID <- levels(colData(sce)[[batchKey]])
        rStep <- 10*(ceiling(max(sce$row.orig,na.rm=T)/10)+1)
        cStep <- 10*(ceiling(max(sce$col.orig,na.rm=T)/10)+1)
        sampleLayout <- matrix(c(sID), nrow=round(sqrt(length(sID))),ncol=ceiling(sqrt(length(sID))),byrow=T)
        sampleLayout[duplicated(as.vector(sampleLayout))] <- NA
        for(one in sID){
            message("\t\t",one)
            rc <- which(sampleLayout==one,arr.ind=T)
            sel <- colData(sce)[[batchKey]]==one
            sce$row[sel] <- rStep*rc[1]-sce$row.orig[sel] # flip the vertical to keep consistency with SpaGCN
            sce$col[sel] <- cStep*(rc[2]-1)+sce$col.orig[sel]
        }
        
        # normalization
        message("\tlognormalize & PCA")
        sce <- spatialPreprocess(sce,platform="Visium", n.PCs = 50) #lognormalize, PCA
        message("\tUMAP")
        sce <- runUMAP(sce, dimred = "PCA")
        colnames(reducedDim(sce, "UMAP")) = c("UMAP1", "UMAP2")
        
        # Harmonization
        message("\tHarmonization")
        sce = RunHarmony(sce, batchKey, verbose = T)
        sce = runUMAP(sce, dimred = "HARMONY", name = "HARMONY_umap")
        colnames(reducedDim(sce, "HARMONY_umap")) = c("UMAP1", "UMAP2")
        
        #BayesSpace tune clusters
        message("\tBayesSpace tune clusters")
        mclust1 = Mclust(reducedDim(sce, "HARMONY")[,1:20], G = 1:20, modelNames = "EEE")
        sce = qTune(sce, qs = seq(2,15), use.dimred = "HARMONY", nrep = 2000)
        pdf(gsub("rds","BIC.pdf",strRaw))
        plot(mclust1$BIC)
        a <- dev.off()
        saveRDS(sce,strProcess)
    }
    
    clusters <- BayesSpace_cluster(sce,config,gsub("\\.rds$","_cluster",strOut))
    saveRDS(clusters,strOut)
    colData(sce) <- cbind(colData(sce),clusters)
    return(sce)
}
BayesSpace_cluster <- function(sce,config,strPrefix){
    message("BayesSpace cluster ...")
    clusterN <- config$clusterN
    clusters <- bplapply(clusterN,function(q){
        strF <- paste0(strPrefix,q,".rds")
        if(file.exists(strF)){
            message("\tUsing previous cluster results: ",strF,"\n\t If a new process is wanted, please remove/rename the above file")
            oneCluster <- readRDS(strF)
        }else{
            message("\tClusterN: ",q)
            oneSCE <- spatialCluster(sce,q=q,use.dimred="HARMONY",nrep=25000)
            oneCluster <- setNames(data.frame(row.names=colnames(oneSCE),A=paste0("C",oneSCE$spatial.cluster)),paste0(BayesSpace_clusterKey,q))
            saveRDS(oneCluster,strF)
        }
        return(oneCluster)
    },BPPARAM = MulticoreParam(workers=min(5,length(clusterN),max(1,parallelly::availableCores()-2)),
                               tasks=length(clusterN)))
    return(do.call(cbind,clusters))
}

BayesSpace_plot <- function(sce,strPDF){
    sID <- levels(colData(sce)[[batchKey]])
    pdf(strPDF,width=3*ceiling(sqrt(length(sID)))+2,height=3*round(sqrt(length(sID))))
    print(clusterPlot(sce,color=NA,label=batchKey)+labs(fill = "Capture\narea"))

    print(ggplot(data.frame(reducedDim(sce, "UMAP")),
           aes(x = UMAP1, y = UMAP2, color = colData(sce)[[batchKey]])) +
        geom_point(alpha = 0.1) +
        ggtitle("Uncorrected")+
        labs(color = batchKey) +
        guides(color = guide_legend(override.aes = list(alpha = 1))) +
        theme_bw())
    
    print(ggplot(data.frame(reducedDim(sce, "HARMONY_umap")),
           aes(x = UMAP1, y = UMAP2, color = colData(sce)[[batchKey]])) +
        geom_point(alpha = 0.1) +
        ggtitle("Harmony")+
        labs(color = batchKey) +
        guides(color = guide_legend(override.aes = list(alpha = 1))) +
        theme_bw())
    qPlot(sce)
    for(one in grep(paste0("^",BayesSpace_clusterKey),colnames(colData(sce)),value=T)){
        print(clusterPlot(sce,label=one,color=NA)+ggtitle(one))
    }
    a<-dev.off()
}

main <- function(){
    suppressMessages(suppressWarnings(loadBayesPkg()))
    args = commandArgs(trailingOnly=TRUE)
    strConfig <- args[1]
    strOut <- args[2]
    
    list[config,meta] <- getConfig(strConfig)
    sr_checkMeta(config,meta)
    dir.create(dirname(strOut),showWarnings=F)
    
    strRaw <- gsub("\\.rds$","_raw.rds",strOut)
    sr_read(meta,config$sample_name,strRaw)
    
    sce <- BayesSpace_process(config,strRaw,strOut)
    
    strPDF <- gsub("rds$","pdf",strOut)
    BayesSpace_plot(sce,strPDF)
    message("\n\n")
}

main()