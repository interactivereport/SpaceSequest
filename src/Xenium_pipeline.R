#Xenium_pipe.R
#Xenium processing pipeline
#Input: config.yml file
#Output: various results

# CurrPath <- dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(CurrPath)

library(yaml)
library(Seurat)
library(progressr)
library(glmGamPoi)
library(Azimuth)    #remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
library(ggplot2)

args <- commandArgs(TRUE)
Config <- args[1]
i <- as.integer(args[2])

yaml <- yaml.load_file(Config)
RunAnno <- FALSE

#Preparation
sampleMeta <- read.csv(yaml$sampleMeta)
if (! is.null(yaml$reference)){
  RunAnno = TRUE
  if (is.null(yaml$reference_name)){
    print("No reference_name specified. Turn off reference mapping.")
    RunAnno = FALSE
  }
}
Sample_number <- dim(sampleMeta)[1]
dir.create(yaml$output_dir, showWarnings = FALSE)

print("##################################################################")
print("Starting Xenium data processing pipeline.")
print("")

# for (i in 1:Sample_number){}  #for loop disabled because 'Data <- subset(Data, subset = nCount_Xenium > 0)' won't run in the second for loop.

Current_meta <- sampleMeta[i,]
Current_name <- Current_meta$Sample
Current_dir <- Current_meta$Directory
if(dim(Current_meta)[2] > 2){
  Current_meta_rest <- Current_meta[,3:ncol(Current_meta)]
}

print(paste0("Processing sample: ",i, ", ", Current_name))

print("Data loading...")
Data <- LoadXenium(Current_dir, fov = "fov")
print("Data loaded")
print("Filterig data using nCount_Xenium > 0...")
Data <- subset(Data, subset = nCount_Xenium > 0)
CellNumber <- dim(Data@meta.data)[1]
Current_meta_rest_dup <- Current_meta_rest[rep(1, CellNumber), ]
Data@meta.data <- cbind(Data@meta.data, Current_meta_rest_dup)

plot1 <- VlnPlot(Data, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.nFeature_nCount_Vln.pdf"), plot1)
png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.nFeature_nCount_Vln.png"), width = 800, height = 800)
print(plot1)
dev.off()
print("QC figures (nFeature_Xenium, nCount_Xenium) generated.")

print("Running SCTransform...")
Data <- SCTransform(Data, assay = "Xenium")
print("Running PCA...")
Data <- RunPCA(Data, npcs = 30, features = rownames(Data))
print("Running UMAP...")
Data <- RunUMAP(Data, dims = 1:30)
print("Identifying clusters...")
Data <- FindNeighbors(Data, reduction = "pca", dims = 1:30)
Data <- FindClusters(Data, resolution = yaml$cluster_resolution)
print("Plotting clusters...")

plot2 <- DimPlot(Data)
ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.original_cluster.UMAP.pdf"), plot2)
png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.original_cluster.UMAP.png"), width = 800, height = 800)
print(plot2)
dev.off()

plot3 <- ImageDimPlot(Data, cols = "polychrome", size = 0.75)
ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.original_cluster.tissue.pdf"), plot3)
png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.original_cluster.tissue.png"), width = 800, height = 800)
print(plot3)
dev.off()

if (RunAnno == TRUE){
  print("Reference was assigned in the config file.")
  print("Performing label transfer to annotate cell types...")
  print(paste0("Using reference: ", yaml$reference))
  print(paste0("Cell type annotation column name: ", yaml$reference_name))
  print("Running cell type label transfer...")
  tryCatch(
    expr = {Data <- RunAzimuth(Data, reference = yaml$reference)}, #ref: https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html
    error=function(e) {print("Error occurred in label transfer. Please consider updating seurat-data and azimuth. They have versions that are compatible with Seurat5")}
  )
  RefName <- paste0("predicted.",yaml$reference_name)
  if (RefName %in% colnames(Data@meta.data)){
    print("Cell type label transfer was successful.")
    Idents(Data) <- Data@meta.data[[RefName]]
    
    print("Plotting figures using annotated cell types...")
    plot4 <- DimPlot(Data)
    ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.annotated.UMAP.pdf"), plot4)
    png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.annotated.UMAP.png"), width = 800, height = 800)
    print(plot4)
    dev.off()
    
    plot5 <- ImageDimPlot(Data, cols = "polychrome", size = 0.75)
    ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.annotated.tissue.pdf"), plot5)
    png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.annotated.tissue.png"), width = 800, height = 800)
    print(plot5)
    dev.off()
  }
}
print("Saving Rdata file...")
save(Data, file = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".RData"))
remove(Data, plot1, plot2, plot3, plot4, plot5)
print("Processing finished.")
print("******************************************************************")
print("")

print("##################################################################")