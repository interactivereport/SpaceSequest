#VisiumHD_pipeline.R
#VisiumHD processing pipeline
#Input: config.yml file
#Output: various results, figures, and an Rdata file for each input sample.

# CurrPath <- dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(CurrPath)

library(yaml)
library(Seurat)
library(ggplot2)
library(Matrix)
library(Azimuth)
library(spacexr)

args <- commandArgs(TRUE)
Config <- args[1]
srcdir <- args[2]
i <- as.integer(args[3])

yaml <- yaml.load_file(Config)
# yaml <- yaml.load_file("/edgehpc/dept/compbio/projects/SpaceSequest/VisiumHD_public/config.yml")
RunAnno <- FALSE
binsize=8
hires <- FALSE
Cluster_resolution=as.numeric(yaml$cluster_resolution)

#Preparation
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
print("Starting Visium HD data processing pipeline.")
print(paste0("There are total ",Sample_number, " samples to process."))
print("")

Current_meta <- sampleMeta[i,]
Current_name <- Current_meta$Sample
Current_dir <- Current_meta$Directory
Hires_img <- Current_meta$Hires_image
if (is.null(Hires_img)){
  hires=FALSE
  print("No high-resolution image provided")
}else{
  hires=TRUE
}
print(paste0("Processing sample ", i, ": ", Current_name))

print("Data loading...")
if (hires){
  Data <- Load10X_Spatial(data.dir = Current_dir, bin.size = c(binsize), image = hires)
}else{
  Data <- Load10X_Spatial(data.dir = Current_dir, bin.size = c(binsize))
}
options(future.globals.maxSize = 800000 * 1024^2)
print("Loaded assays:")
Assays(Data)
DefaultAssay(Data) <- spatial_assay
print("Plotting QC figures for nCount and nFeature")
plot1 <- SpatialFeaturePlot(Data, features = nCount_assay)
ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.counts.pdf"), plot1)
png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.counts.png"), width = 2500, height = 2500)
print(plot1)
dev.off()

plot2 <- VlnPlot(Data, features = c(nFeature_assay, nCount_assay), pt.size = 0)
ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.count_feature.pdf"), plot2)
png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.count_feature.png"), width = 800, height = 800)
print(plot2)
dev.off()

print("Running Normalization...")
Data <- NormalizeData(Data)
print("Running FindVariableFeatures...")
Data <- FindVariableFeatures(Data)
print("Running ScaleData...")
Data <- ScaleData(Data)

print("Selecting 50,000 cells for the sketch assay...")
Data <- SketchData(object = Data, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(Data) <- "sketch"

print("Running sketch FindVariableFeatures...")
Data <- FindVariableFeatures(Data)
print("Running sketch ScaleData...")
Data <- ScaleData(Data)
print("Running sketch PCA...")
Data <- RunPCA(Data, assay = "sketch", reduction.name = "pca.sketch")
print("Running sketch FindNeighbors...")
Data <- FindNeighbors(Data, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
print("Running sketch FindClusters...")
Data <- FindClusters(Data, cluster.name = "seurat_cluster.sketched", resolution = Cluster_resolution)
print("Running sketch UMAP...")
Data <- RunUMAP(Data, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

print("Projecting sketch results to the full dataset...")
Data <- ProjectData(
  object = Data,
  assay = spatial_assay,
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

print("Plotting UMAP using the sketch dataset...")
DefaultAssay(Data) <- "sketch"
Idents(Data) <- "seurat_cluster.sketched"
plot3 <- DimPlot(Data, reduction = "umap.sketch", label = F) + ggtitle(paste0(Current_name, ", Sketched clustering (50,000 cells)")) + theme(legend.position = "bottom")
ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.UMAP.sketch.pdf"), plot3)
png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.UMAP.sketch.png"), width = 800, height = 800)
print(plot3)
dev.off()

print("Plotting UMAP using the full dataset...")
DefaultAssay(Data) <- spatial_assay
Idents(Data) <- "seurat_cluster.projected"
plot4 <- DimPlot(Data, reduction = "full.umap.sketch", label = F) + ggtitle(paste0(Current_name, ", Projected clustering (full dataset)")) + theme(legend.position = "bottom")
ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.UMAP.full.pdf"), plot4)
png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.UMAP.full.png"), width = 800, height = 800)
print(plot4)
dev.off()

plot5 <- SpatialDimPlot(Data, label = T, repel = T, label.size = 4)
ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.Spatial.cluster.pdf"), plot5)
png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.Spatial.cluster.png"), width = 1500, height = 1500)
print(plot5)
dev.off()

#################Reference annotation
if (RunAnno == TRUE){
  print("Starting reference-based cell type annotation")
  print("Loading refernce...")
  refdata <- readRDS(yaml$reference)
  
  #Filtering refdata. Each cell type needs to have at least 25 cells.
  print("Filtering data to make sure each cell type has at least 25 cells.")
  Cell_counts <- table(refdata@meta.data[[yaml$reference_name]])
  Keep_cell_type <- names(Cell_counts[Cell_counts>25])
  refdata@meta.data$keep_anno <- refdata@meta.data[[yaml$reference_name]] %in% Keep_cell_type
  tf<- refdata@meta.data[[yaml$reference_name]] %in% Keep_cell_type
  refdata <- subset(refdata, subset=keep_anno==TRUE)
  
  print("Formatting reference...")
  Idents(refdata) <- yaml$reference_name
  counts <- refdata[["RNA"]]$counts
  cluster <- as.factor(refdata@meta.data[[yaml$reference_name]])
  nUMI <- refdata$nCount_RNA
  levels(cluster) <- gsub("/", "-", levels(cluster))
  cluster <- droplevels(cluster)
  names(cluster) <- names(nUMI)
  
  # create the RCTD reference object
  print("Preparing RCTD reference object...")
  reference <- Reference(counts, cluster, nUMI)
  
  DefaultAssay(Data) <- "sketch"
  counts_hd <- Data[["sketch"]]$counts
  data_cells_hd <- colnames(Data[["sketch"]])
  coords <- GetTissueCoordinates(Data)[data_cells_hd, 1:2]
  
  # create the RCTD query object
  print("Preparing RCTD query object...")
  query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
  
  # run RCTD
  print("Running RCTD, using doublet mode...")
  print("This step may take a long time...")
  RCTD <- create.RCTD(query, reference, max_cores = 28)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  # add results back to Seurat object
  
  print("Formatting RCTD results...")
  Data <- AddMetaData(Data, metadata = RCTD@results$results_df)

  # project RCTD labels from sketched cells to the full dataset
  print("Projecting RCTD labels on the full dataset")
  Data$cell_type <- as.character(Data$first_type)
  Data$cell_type[is.na(Data$cell_type)] <- "Unknown"
  Data <- ProjectData(
    object = Data,
    assay = spatial_assay,
    full.reduction = "full.pca.sketch",
    sketched.assay = "sketch",
    sketched.reduction = "pca.sketch",
    umap.model = "umap.sketch",
    dims = 1:50,
    refdata = list(cell_type = "cell_type")
  )
  
  if ("cell_type" %in% colnames(Data@meta.data)){
    print("Cell type label transfer was successful.")
    Idents(Data) <- Data$cell_type
    plot6 <- SpatialDimPlot(Data, label = T, repel = T, label.size = 4)
    ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.Spatial.annotated_cluster.pdf"), plot6)
    png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.Spatial.annotated_cluster.png"), width = 1500, height = 1500)
    print(plot6)
    dev.off()
  }
}
print("Saving Rdata file...")
save(Data, file = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".RData"))
remove(Data, plot1, plot2, plot3, plot4, plot5)
print("Processing finished")
print("******************************************************************")
print("")

print("##################################################################")