#CosMx_merger.R
#CosMx data merging script
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
print("Starting data integration workflow.")

DataList=list()
for (i in 1:Sample_number){
  print(paste0("Processing individual data: ",i))
  Current_meta <- sampleMeta[i,]
  Current_name <- Current_meta$Sample
  Current_dir <- Current_meta$Directory
  
  load(paste0(yaml$output_dir, "/", yaml$project_ID, ".", Current_name, ".RData"))
  Data@meta.data$dataset <- Current_name
  DataList[i] = Data
  names(DataList)[i] <- Current_name
  rm(Data)
}

#The first parameter of merge should be a Seurat object, the second (y) can be one Seurat object or a list of several.
options(future.globals.maxSize = 800000 * 1024^2)
MergedData <- merge(x = DataList[[1]], y = DataList[-1])

print("Running SCTransform...")
MergedData <- SCTransform(MergedData, assay = "Nanostring", clip.range = c(-10, 10), verbose = TRUE)
print("Running PCA...")
MergedData <- RunPCA(MergedData, npcs = 30, features = rownames(MergedData))
if (use_harmony == TRUE && Sample_number>1){
  print("Running Harmony...")
  MergedData <- IntegrateLayers(
    object = MergedData, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    verbose = TRUE
  )
  print("Running UMAP...")
  MergedData <- RunUMAP(MergedData, dims = 1:30, reduction = "harmony")
  print("Identifying clusters...")
  MergedData <- FindNeighbors(MergedData, reduction = "harmony", dims = 1:30)
}else{
  print("Running UMAP...")
  MergedData <- RunUMAP(MergedData, dims = 1:30)
  print("Identifying clusters...")
  MergedData <- FindNeighbors(MergedData, reduction = "pca", dims = 1:30)
}

MergedData <- FindClusters(MergedData, resolution = yaml$cluster_resolution)
print("Plotting clusters...")

plot1 <- DimPlot(MergedData, group.by = "dataset")
ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.plot.dataset.pdf"), plot1)
png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.plot.dataset.png"), width = 800, height = 800)
print(plot1)
dev.off()

if (RunAnno == TRUE){
  RefName <- paste0("predicted.",yaml$reference_name)
  if (RefName %in% colnames(MergedData@meta.data)){
    print("Use cell type labels to plot UMAP.")
    Idents(MergedData) <- MergedData@meta.data[[RefName]]
    
    plot2 <- DimPlot(MergedData)
    ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.plot.annotated.pdf"), plot2)
    png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.plot.annotated.png"), width = 800, height = 800)
    print(plot2)
    dev.off()
  }
}else{
  plot2 <- DimPlot(MergedData)
  ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.plot.original_cluster.pdf"), plot2)
  png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.plot.original_cluster.png"), width = 800, height = 800)
  print(plot2)
  dev.off()
}

print("Saving the integrated Rdata file...")
save(MergedData, file = paste0(yaml$output_dir,"/",yaml$project_ID,".integrated.RData"))
remove(MergedData, plot1, plot2)
print("Processing finished.")
print("******************************************************************")
print("")

print("##################################################################")