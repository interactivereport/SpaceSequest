#GeoMx processing pipeline
#Input: config.yml file
#Output: various results including Quickomics files, DEG results, etc.

# CurrPath <- dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(CurrPath)

library(yaml)
library(readxl)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(knitr)
library(scales)
library(reshape2)
library(cowplot)
library(umap)
library(Rtsne)
library(pheatmap)
library(preprocessCore)

args <- commandArgs(TRUE)
yaml <- yaml.load_file(args[1])
RunDE <- FALSE

#Preparation
comparison <- read.csv(yaml$comparison)
if (dim(comparison)[1] > 0){
  RunDE = TRUE
}

#Data loading
print("Starting the GeoMx pipeline")
dcc_files <- dir(yaml$data_path, pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
pkc_files <- yaml$pkcs_file
data_annotation <- yaml$data_annotation
sheet <- yaml$annotation_sheet
sheet_excel <- read_xlsx(yaml$data_annotation, sheet = sheet)
dir.create(yaml$output_dir, showWarnings = FALSE)

#Data loading
print("Loading dcc files")
Data <-
  readNanoStringGeoMxSet(dccFiles = dcc_files,
                         pkcFiles = pkc_files,
                         phenoDataFile = data_annotation,
                         phenoDataSheet=sheet,
                         phenoDataDccColName = "Sample_ID",
                         protocolDataColNames = c("aoi", "roi"),
                         experimentDataColNames = c("panel"))
pkcs <- annotation(Data)
modules <- gsub(".pkc", "", pkcs)
data.frame(PKCs = pkcs, modules = modules)

#Check data loaded
Loaded_data <- gsub( ".dcc", "", colnames(Data))
Sample_meta_data <- sheet_excel$Sample_ID
if (sum(! Sample_meta_data %in% Loaded_data) > 0){
  print("The following files in the sample sheet were not loaded.")
  Sample_meta_data[! Sample_meta_data %in% Loaded_data]
}else{
  print("All data in the sample sheet were loaded.")
}

# Shift counts to one
Data <- shiftCountsOne(Data, useDALogic = TRUE)

##################################################################
# Default QC cutoffs
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100)
       minArea = 1000)         # Minimum segment area (5000)
Data <- setSegmentQCFlags(Data, qcCutoffs = QC_params)

# Collate QC Results
QCResults <- protocolData(Data)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

#Remove WARNING data: 25. Kept the PASS QC data
Data <- Data[, QCResults$QCStatus == "PASS"]
print("Loaded data: Features and Samples")
dim(Data)
kable(QC_Summary, caption = "QC Summary by Segment")


##################################################################
#Probe QC
# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
print("Running probe QC")
Data <- setBioProbeQCFlags(Data, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)
ProbeQCResults <- fData(Data)[["QCFlags"]]
# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

####Exclude outlier probes
#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(Data, 
         fData(Data)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(Data)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
print("After probe QC: Features and Samples")
dim(ProbeQCPassed)
Data <- ProbeQCPassed 

##################################################################
####Create gene-level counts
print("Creating gene-level counts")
# Check how many unique targets the object has
print("Total target number removing negative probes:")
length(unique(featureData(Data)[["TargetName"]]))
target_Data <- aggregateCounts(Data)
dim(target_Data)

##################################################################
#Limit of Quantification
print("Estimating Limit of Quantification (LOQ)")
#Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2
# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_Data))
for(module in modules) {
  # print(module)
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),       #The meta data frame has those values already
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_Data)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_Data)[, vars[1]] * 
             pData(target_Data)[, vars[2]] ^ cutoff)
  }
}
pData(target_Data)$LOQ <- LOQ


##################################################################
print("Performing gene detection rate calculation based on LOQ")
#Another round of filtering, gene detection rate
#Prepare filtering
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_Data)$Module == module
  Mat_i <- t(esApply(target_Data[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
#LOQ_Mat is a True/False matrix. True means passed LOQ cutoff
LOQ_Mat <- LOQ_Mat[fData(target_Data)$TargetName, ]
# Gene detection rate:
pData(target_Data)$GenesDetected <- colSums(LOQ_Mat, na.rm = TRUE)
pData(target_Data)$GeneDetectionRate <- pData(target_Data)$GenesDetected / nrow(target_Data)
#feature level:
LOQ_Mat <- LOQ_Mat[, colnames(target_Data)]
fData(target_Data)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_Data)$DetectionRate <- fData(target_Data)$DetectedSegments / nrow(pData(target_Data))

#Keep samples with gene detection rate >= 0.05
print("Filtering samples with gene detection rate < 5%")
target_Data <-
  target_Data[, pData(target_Data)$GeneDetectionRate >= .05]
print("After sample filtering: Features and Samples")
dim(target_Data)


##################################
#Normalization
# Graph Q3 value vs negGeoMean of Negatives
# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
print("Normalizing data")
print("Performing Q3 normalization")
target_Data <- normalize(target_Data ,
                             norm_method = "quant", 
                             desiredQuantile = .75,
                             toElt = "q_norm")
print("Performing background normalization")
target_Data <- normalize(target_Data ,
                             norm_method = "neg", 
                             fromElt = "exprs",
                             toElt = "neg_norm")
countFile <- as.data.frame(exprs(target_Data))
#Remove NegProbe-WTX
countFile <- countFile[row.names(countFile) != "NegProbe-WTX",]
target_Data_quantile <- normalize.quantiles(as.matrix(countFile))
row.names(target_Data_quantile) <- row.names(countFile)
colnames(target_Data_quantile) <- colnames(countFile)
# boxplot(target_Data_quantile[,1:10],
#         col = "orchid", main = "Our data, Quantile Normalized, filtered",
#         log = "y", names = 1:10, xlab = "Segment",
#         ylab = "Counts, Quantile Normalized")

##################################
#DE analysis
Meta <- pData(target_Data)
Meta <- Meta[,1:ncol(Meta)-1]  #drop neg norm factors. It may cause issues when loading the data frame

#Prepare Quickomics results
QuickomicsPrefix=yaml$project_ID
#File 1: data
if (yaml$quickomics == TRUE){
  write.csv(countFile[row.names(countFile) != "NegProbe-WTX",], file = paste0(yaml$output_dir, "/", QuickomicsPrefix, "_Exp_RawData.csv"))
  write.csv(countFile[row.names(countFile) != "NegProbe-WTX",], file = paste0(yaml$output_dir, "/", QuickomicsPrefix, "_Exp_Q3NormData.csv"))
  write.csv(countFile[row.names(countFile) != "NegProbe-WTX",], file = paste0(yaml$output_dir, "/", QuickomicsPrefix, "_Exp_NegNormData.csv"))
  write.csv(target_Data_quantile, file = paste0(yaml$output_dir, "/", QuickomicsPrefix, "_Exp_QuantileNormData.csv"))
}
#File 2: sample meta
#Construct metadata: adding sampleid column, and point the group column to the region/segment/class column, if any of them exists.
samplemeta <- data.frame(row.names = row.names(Meta))
samplemeta$sampleid <- row.names(samplemeta)
if ("region" %in% colnames(Meta)){
  print("Using region as the group column")
  samplemeta$group = Meta$region
}else if ("segment" %in% colnames(Meta)){
  print("Using segment as the group column")
  samplemeta$group = Meta$segment
}else if ("class" %in% colnames(Meta)){
  print("Using class as the group column")
  samplemeta$group = Meta$class
}else{
  print("Using the first column in the meta file as the group column")
  #None of the above three exists, use the first column
  samplemeta$sampleid = Meta[,1]
}
samplemeta <- cbind(samplemeta, Meta)
samplemeta$LOQ <- unlist(samplemeta$LOQ) #This column needs to be unlisted
if (yaml$quickomics == TRUE){
  write.csv(samplemeta, file = paste0(yaml$output_dir, "/", QuickomicsPrefix, "_Sample_metadata.csv"), row.names = FALSE)
}
#File 3: DE results
de_results <- data.frame()
if (RunDE == FALSE){
  #No DE step to run. output a dummy DE file for Quickomics
  if (yaml$quickomics == TRUE){
    de_results <- data.frame(row.names = c("empty"))
    de_results$UniqueID <- row.names(de_results)
    de_results$test <- c("empty")
    de_results$Adj.P.Value <- 1
    de_results$P.Value <- 1
    de_results$logFC <- 1
    write.csv(de_results, file = paste0(yaml$output_dir, "/", QuickomicsPrefix, "_Comparison_Data.csv"), row.names = FALSE)
  }
}
#File 4: Optional gene matching file
matching <- data.frame(row.names = row.names(countFile[row.names(countFile) != "NegProbe-WTX",]))
matching$id <- 1:length(row.names(matching))
matching$UniqueID <- row.names(matching)
matching$Gene.Name <- row.names(matching)
matching$ProteinID <- ""
matching$GeneType <- "protein_coding"
if (yaml$quickomics == TRUE){
  write.csv(matching, file = paste0(yaml$output_dir, "/", QuickomicsPrefix, "_ProteinGeneName_Optional.csv"), row.names = FALSE)
}


if (RunDE == TRUE){
  print("Running DE: Comparison file is not empty")
  
  #Loop through the comparison file
  for (i in 1:dim(comparison)[1]){
    print(paste("Working on comparison group:",i))
    Current_comparison <- comparison[i,]
    print(Current_comparison)
    #Selecting Group_test
    Meta_test <- Meta[Meta[[Current_comparison$Group_name]] == Current_comparison$Group_test,]
    Meta_control <- Meta[Meta[[Current_comparison$Group_name]] == Current_comparison$Group_ctrl,]
    Meta_combine_test_control <- rbind(Meta_test, Meta_control)
    print("Samples in: Group test | Group control")
    table(Meta[[Current_comparison$Group_name]])
    print("Subsetting data, using Q3 norm values")
    Data <- assayDataElement(target_Data, elt = "q_norm")
    Data_test <- Data[,row.names(Meta_test)]
    Data_control <- Data[,row.names(Meta_control)]
    
    #Preparing data
    Data_test_control_combined <- cbind(Data_test, Data_control)
    group <- factor(c(rep(Current_comparison$Group_test, dim(Meta_test)[1]), rep(Current_comparison$Group_ctrl, dim(Meta_control)[1])),
                    levels = c(Current_comparison$Group_ctrl, Current_comparison$Group_test))
    results <- data.frame(gene = rownames(Data_test_control_combined), 
                          Log2FC = numeric(nrow(Data_test_control_combined)), 
                          p.value = numeric(nrow(Data_test_control_combined)), 
                          stringsAsFactors = FALSE)
    anno <- data.frame(row.names = c(row.names(Meta_test), row.names(Meta_control)))
    anno$group <- c(rep(Current_comparison$Group_test, dim(Meta_test)[1]), rep(Current_comparison$Group_ctrl, dim(Meta_control)[1]))
    anno_order <- anno[colnames(Data_test_control_combined),]
    
    #Running DE tests
    if (Current_comparison$Analysis_method == "Linear"){
      
    }else{
      print("Sorry, we don't support this DE method. Exit the workflow.")
    }
    for (j in 1:nrow(Data_test_control_combined)) {
      #print(j)
      # Extract expression data for the current gene
      expression_data <- Data_test_control_combined[j, ]
      
      # Fit linear model
      fit <- lm(expression_data ~ group)
      
      # Extract coefficient for treated vs control
      # coef(fit)["groupMGL_Near"]
      
      # Extract p-value for the coefficient
      results$p.value[j] <- summary(fit)$coefficients[paste0("group",Current_comparison$Group_test), "Pr(>|t|)"]
    }
    
    results$adj.p.value <- p.adjust(results$p.value, method = "BH")
    
    #Get mean expression in two groups
    Mean_test <- apply(Data_test, 1, FUN = mean)
    Mean_ctrl <- apply(Data_control, 1, FUN = mean)
    results[[paste0("Mean_",Current_comparison$Group_test)]] <- Mean_test
    results[[paste0("Mean_",Current_comparison$Group_ctrl)]] <- Mean_ctrl
    results$Log2FC <- log2(results[[paste0("Mean_",Current_comparison$Group_test)]]/results[[paste0("Mean_",Current_comparison$Group_ctrl)]])
    results_ordered <- results[order(results$adj.p.value),]
    results_ordered$test <- Current_comparison$CompareName
    write.table(results_ordered, file = paste0(yaml$output_dir,"/Comparison_result_",i,".txt"), quote = F, sep = "\t")
    #Formatting results for Quickomics:
    results_formatted <- data.frame(row.names = row.names(results))
    results_formatted$UniqueID <- results$gene
    results_formatted$test <- Current_comparison$CompareName
    results_formatted$Adj.P.Value <- results$adj.p.value
    results_formatted$P.Value <- results$p.value
    results_formatted$logFC <- results$Log2FC
    if (dim(de_results)[1] == 0){
      de_results=results_formatted
    }else{
      de_results=rbind(de_results, results_formatted)
    }
  }
}

if (yaml$quickomics == TRUE){
  write.csv(de_results, file = paste0(yaml$output_dir, "/", QuickomicsPrefix, "_Comparison_Data.csv"), row.names = FALSE)
}

print("GeoMx pipeline finished")
print("##################################################################")
#End of the GeoMx workflow
##################################
