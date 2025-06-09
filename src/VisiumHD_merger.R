#VisiumHD_merger.R
#VisiumHD data merging script
#Input: config.yml file
#Output: merged Seurat object, and associated plots

# CurrPath <- dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(CurrPath)

library(yaml)
library(Seurat)
library(ggplot2)
library(harmony)
# library(progressr)
# library(glmGamPoi)
# library(Azimuth)    #remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)

args <- commandArgs(TRUE)
Config <- args[1]
use_harmony=TRUE

yaml <- yaml.load_file(Config)
IntegrateData <- yaml$integrate_data
RunAnno = FALSE

if (IntegrateData == FALSE){
  print("Skip data integrating step. Workflow ends.")
  print("##################################################################")
}

if (yaml$integrate_with_harmony == FALSE){
  use_harmony=FALSE
}

#Preparation
Cluster_resolution=as.numeric(yaml$cluster_resolution)
sampleMeta <- read.csv(yaml$sampleMeta)
if (! is.null(yaml$reference)){
  RunAnno = TRUE
  if (is.null(yaml$reference_name)){
    print("No reference_name specified. Turn off reference mapping.")
    RunAnno = FALSE
  }
}
if (yaml$bin_resolution == "8um"){
  binsize=8
  spatial_assay="Spatial.008um"
  nCount_assay="nCount_Spatial.008um"
  nFeature_assay="nFeature_Spatial.008um"
}else if (yaml$bin_resolution == "16um"){
  binsize=16
  spatial_assay="Spatial.016um"
  nCount_assay="nCount_Spatial.016um"
  nFeature_assay="nFeature_Spatial.016um"
}else if (yaml$bin_resolution == "2um"){
  binsize=2
  spatial_assay="Spatial.002um"
  nCount_assay="nCount_Spatial.002um"
  nFeature_assay="nFeature_Spatial.002um"
}else{
  print("Error in selecting bin resolutions. Please use 2um, 8um, or 16um. Using 8um as default.")
  binsize=8
  spatial_assay="Spatial.008um"
  nCount_assay="nCount_Spatial.008um"
  nFeature_assay="nFeature_Spatial.008um"
}
Sample_number <- dim(sampleMeta)[1]
dir.create(yaml$output_dir, showWarnings = FALSE)

print("##################################################################")
print("Starting data integration workflow.")

DataList=list()
print(paste0("Total number of data to integrate: ", Sample_number))
for (i in 1:Sample_number){
  print(paste0("Processing individual data: ",i))
  Current_meta <- sampleMeta[i,]
  Current_name <- Current_meta$Sample
  Current_dir <- Current_meta$Directory
  
  load(paste0(yaml$output_dir, "/", yaml$project_ID, ".", Current_name, ".RData"))
  Data@meta.data$dataset <- Current_name
  DefaultAssay(Data) <- spatial_assay
  DataList[i] = Data
  names(DataList)[i] <- Current_name
  rm(Data)
}

#The first parameter of merge should be a Seurat object, the second (y) can be one Seurat object or a list of several.
options(future.globals.maxSize = 800000 * 1024^2)
MergedData <- merge(x = DataList[[1]], y = DataList[-1])

print("Running Normalization...")
MergedData <- NormalizeData(MergedData)
print("Running FindVariableFeatures...")
MergedData <- FindVariableFeatures(MergedData)
print("Running ScaleData...")
MergedData <- ScaleData(MergedData)

print("Selecting 50,000 cells for the sketch assay...")
MergedData <- SketchData(object = MergedData, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(MergedData) <- "sketch"

print("Running sketch FindVariableFeatures...")
MergedData <- FindVariableFeatures(MergedData)
print("Running sketch ScaleData...")
MergedData <- ScaleData(MergedData)
print("Running sketch PCA...")
MergedData <- RunPCA(MergedData, assay = "sketch", reduction.name = "pca.sketch")
if (use_harmony == TRUE && Sample_number>1){
  print("Running Harmony...")
  MergedData <- IntegrateLayers(
    object = MergedData, method = HarmonyIntegration,
    orig.reduction = "pca.sketch", new.reduction = "harmony",
    verbose = TRUE
  )
  print("Running sketch FindNeighbors...")
  MergedData <- FindNeighbors(MergedData, assay = "sketch", reduction = "harmony", dims = 1:50)
  print("Running sketch FindClusters...")
  MergedData <- FindClusters(MergedData, cluster.name = "seurat_cluster.sketched", reduction = "harmony", resolution = Cluster_resolution)
  print("Running sketch UMAP...")
  MergedData <- RunUMAP(MergedData, reduction = "harmony", reduction.name = "umap.sketch", return.model = T, dims = 1:50)
}else{
  print("Running sketch FindNeighbors...")
  MergedData <- FindNeighbors(MergedData, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
  print("Running sketch FindClusters...")
  MergedData <- FindClusters(MergedData, cluster.name = "seurat_cluster.sketched", resolution = Cluster_resolution)
  print("Running sketch UMAP...")
  MergedData <- RunUMAP(MergedData, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)
  
}

print("Projecting sketch results to the full dataset...")
MergedData <- ProjectData(
  object = MergedData,
  assay = spatial_assay,
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

DefaultAssay(MergedData) <- "sketch"
Idents(MergedData) <- "seurat_cluster.sketched"

print("Plotting UMAP using the integrated sketch dataset...")
plot1 <- DimPlot(MergedData, reduction = "umap.sketch", label = F) + ggtitle(paste0("MergedData, ", ", Sketched clustering (50,000 cells)")) + theme(legend.position = "bottom")
ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.plot.UMAP.original_cluster.sketch.pdf"), plot1)
png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.plot.UMAP.original_cluster.sketch.png"), width = 800, height = 800)
print(plot1)
dev.off()

print("Plotting UMAP using the full dataset...")
DefaultAssay(MergedData) <- spatial_assay
Idents(MergedData) <- "seurat_cluster.projected"
plot2 <- DimPlot(MergedData, reduction = "full.umap.sketch", label = F) + ggtitle(paste0("MergedData, ", ", Projected clustering (full dataset)")) + theme(legend.position = "bottom")
ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.plot.UMAP.original_cluster.full.pdf"), plot2)
png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.plot.UMAP.original_cluster.full.png"), width = 800, height = 800)
print(plot2)
dev.off()

print("Plotting UMAP using the full dataset, showing dataset...")
DefaultAssay(MergedData) <- spatial_assay
plot3 <- DimPlot(MergedData, reduction = "full.umap.sketch", group.by="dataset") + ggtitle(paste0("MergedData, ", ", colored by dataset")) + theme(legend.position = "bottom")
ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.plot.UMAP.dataset.full.pdf"), plot3)
png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.plot.UMAP.dataset.full.png"), width = 800, height = 800)
print(plot3)
dev.off()

if ("cell_type" %in% colnames(MergedData@meta.data)){
  print("Found cell type labels. Plotting UMAP using the cell type labels...")
  DefaultAssay(MergedData) <- spatial_assay
  Idents(MergedData) <- "cell_type"
  plot4 <- DimPlot(MergedData, reduction = "full.umap.sketch", label = F) + ggtitle(paste0(Current_name, ", Projected clustering (full dataset)")) + theme(legend.position = "bottom")
  ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.plot.UMAP.annotated.full.pdf"), plot4)
  png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.plot.UMAP.annotated.full.png"), width = 800, height = 800)
  print(plot4)
  dev.off()
}

print("Saving the integrated Rdata file...")
save(MergedData, file = paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.RData"))
remove(MergedData, plot1, plot2, plot3)
print("Processing finished.")
print("******************************************************************")
print("")

print("##################################################################")