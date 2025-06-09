#CosMx_pipe.R
#CosMx processing pipeline
#Input: config.yml file
#Output: various results including integrated data, etc.

# CurrPath <- dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(CurrPath)

library(yaml)
library(Seurat)
library(progressr)
library(glmGamPoi)
library(Azimuth)    #remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
library(ggplot2)

# detachAllPackages <- function() {
#   basic.packages <- c("package:yaml","package:Seurat","package:progressr","package:glmGamPoi","package:Azimuth","package:ggplot2","package:base")
#   package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
#   package.list <- setdiff(package.list,basic.packages)
#   if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
# }

args <- commandArgs(TRUE)
Config <- args[1]
srcdir <- args[2]
i <- as.integer(args[3])
source(paste0(srcdir, "/LoadNanostring_fix.R"))

yaml <- yaml.load_file(Config)
RunAnno <- FALSE
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
Sample_number <- dim(sampleMeta)[1]
dir.create(yaml$output_dir, showWarnings = FALSE)

print("##################################################################")
print("Starting CosMx data processing pipeline.")
print(paste0("There are total ",Sample_number, " samples to process."))
print("")

Current_meta <- sampleMeta[i,]
Current_name <- Current_meta$Sample
Current_dir <- Current_meta$Directory
Current_exprMat <- Current_meta$exprMat_file
Current_fov <- Current_meta$fov_file
Current_metafile <- Current_meta$metadata_file
Current_tx <- Current_meta$tx_file
Current_polygon <- Current_meta$polygon_file
print(paste0("Processing sample ", i, ": ", Current_name))

print("Step 1: Data loading...")
Data <- tryCatch(
  expr = {Data <- LoadNanostring(data.dir = Current_dir, fov = Current_fov, assay = "Nanostring")},
  error=function(e) LoadNanostring_fix(data.dir = Current_dir, fov = Current_fov, assay = "Nanostring")
)

print("Data loaded.")
print("Running SCTransform...")
options(future.globals.maxSize = 800000 * 1024^2) #This is to avoid global msxSize error during label transfer. See: https://github.com/satijalab/seurat/issues/1845, https://satijalab.org/seurat/archive/v3.0/future_vignette.html

Data <- SCTransform(Data, assay="Nanostring", clip.range = c(-10, 10), verbose = TRUE)
print("Running PCA...")
Data <- RunPCA(Data, npcs=50)
print("Running UMAP...")
Data <- RunUMAP(Data, dim=1:50)
print("Identifying clusters...")
Data <- FindNeighbors(Data, reduction="pca", dims=1:50)
Data <- FindClusters(Data, resolution=Cluster_resolution)
print("Plotting clusters...")

plot1 <- ImageDimPlot(Data, fov = Current_fov, axes = TRUE, cols = "glasbey")
ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.original_cluster.pdf"), plot1)
png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.original_cluster.png"), width = 800, height = 800)
print(plot1)
dev.off()

plot2 <- DimPlot(Data, label=TRUE)
ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.original_cluster.UMAP.pdf"), plot2)
png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.original_cluster.UMAP.png"), width = 800, height = 800)
print(plot2)
dev.off()

#Cell type label transfer
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
    
    plot3 <- ImageDimPlot(Data, fov = Current_fov, axes = TRUE, cols = "glasbey")
    ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.annotated.pdf"), plot3)
    png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.annotated.png"), width = 800, height = 800)
    print(plot3) #Use print to allow figure generation inside if/else
    dev.off()
    
    plot4 <- DimPlot(Data, label=TRUE)
    ggsave(paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.annotated.UMAP.pdf"), plot4)
    png(filename = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".plot.annotated.UMAP.png"), width = 800, height = 800)
    print(plot4)
    dev.off()
  }
}
print("Saving Rdata file...")
save(Data, file = paste0(yaml$output_dir,"/",yaml$project_ID,".",Current_name,".RData"))
remove(Data, plot1, plot2, plot3, plot4)
print("Processing finished")
print("******************************************************************")
print("")

print("##################################################################")