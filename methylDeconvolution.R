## ---- Deconvolution analysis based on methylDeConv package.
# R manual: https://github.com/jysonganan/methylDeConv/blob/master/methylDeConv_0.1.1.pdf

## ------------ this code is under construction - I've faced some problems with it!



setwd("~/data")


#____Load methylation and gene expression data____#
load("TCGA_files.RData")
#Only tumor data

#beta_filtered (matrix)
#dmp_NSD1status
#dmp_tumor_ageacc
#dmp_tumor_NSD1
#dmpNSD1_f
#dna_exp_full
#dna_methyl_full
#levine_probes_info
#met_filtered
#NSD1_DMPs_betamatrix
#pd_filtered
#rna
#surv_tumor_tcga


#Call libraries
library(methylDeConv)
library(IlluminaHumanMethylation450kmanifest)
library(rafalib)

##(1)  Build the function
##WARNING: Had to remove the 'unormalized' code (RgChannelSet) because it returned error all the time!!!
# Should check it out with the function developer


MethylDeconv_pipeline <- function(input_methyl, normalized = FALSE, 
                                  array_type = "450k", 
                                  feature_selection = "oneVsAllttest",
                                  deconv_algorithm = "Houseman", 
                                  tissue = "general", 
                                  extend_reference = TRUE,
                                  custom_probes = NULL){
  if(normalized){
    # beta value input
    MethylDeconv_normalized(input_methyl = input_methyl, 
                            array_type = array_type, 
                            feature_selection = feature_selection,
                            deconv_algorithm = deconv_algorithm, 
                            tissue = tissue, 
                            extend_reference = extend_reference,
                            custom_probes = custom_probes)
  }
}




MethylDeconv_normalized <- function(input_methyl, array_type, feature_selection, deconv_algorithm,
                                    tissue, extend_reference, custom_probes){
  if (tissue == "general" && array_type == "EPIC"){
    ref <- build_reference_EPIC(extend = extend_reference)
    ref_betamatrix <- ref$ref_betamatrix
    ref_phenotype <- ref$ref_phenotype
  }
  
  if (tissue == "general" && array_type == "450k"){
    ref <- build_reference_450k(extend = extend_reference)
    ref_betamatrix <- ref$ref_betamatrix
    ref_phenotype <- ref$ref_phenotype
  }
  
  if (tissue == "brain" && array_type == "EPIC"){
    ref <- build_reference_EPIC_neuron()
    ref_betamatrix <- ref$ref_betamatrix
    ref_phenotype <- ref$ref_phenotype
  }
  
  if (tissue == "brain" && array_type == "450k"){
    ref <- build_reference_450k_neuron()
    ref_betamatrix <- ref$ref_betamatrix
    ref_phenotype <- ref$ref_phenotype
  }
  
  ## average reference profiles
  compTable <- ref_compTable(ref_betamatrix, ref_phenotype)
  ntypes <- length(unique(ref_phenotype))
  compTable <- compTable[,3:(2 + ntypes)]
  
  ## feature selection
  if (is.null(custom_probes)){
    if (feature_selection == "oneVsAllttest"){
      probes_select <- ref_probe_selection_oneVsAllttest(ref_betamatrix, ref_phenotype)
    }
    else if (feature_selection == "oneVsAllLimma"){
      probes_select <- ref_probe_selection_oneVsAllLimma(ref_betamatrix, ref_phenotype)
    }
    else if (feature_selection == "pairwiseLimma"){
      probes_select <- ref_probe_selection_pairwiseLimma(ref_betamatrix, ref_phenotype)
    }
    else if (feature_selection == "pairwiseGlmnet"){
      probes_select <- ref_probe_selection_pairwiseGlmnet_cv(ref_betamatrix, ref_phenotype)
    }
    else if (feature_selection == "multiGlmnet"){
      probes_select <- ref_probe_selection_multiclassGlmnet_cv(ref_betamatrix, ref_phenotype)
    }
    else if (feature_selection == "glmnetpreselect"){
      probes_select <- ref_probe_selection_twoStage(ref_betamatrix, ref_phenotype, ml_model = "elastic net")
    }
    else if (feature_selection == "RFpreselect"){
      probes_select <- ref_probe_selection_twoStage(ref_betamatrix, ref_phenotype, ml_model = "RF")
    }
    else if (feature_selection == "OptVariables"){
      probes_select <- ref_probe_selection_twoStage(ref_betamatrix, ref_phenotype, ml_model = "rfe")
    }
  }
  else{
    probes_select = custom_probes
  }
  
  ## deconvolution
  if (deconv_algorithm == "Houseman"){
    deconv_res = Houseman_project(input_methyl, compTable, probes_select)
  }
  else if (deconv_algorithm == "RPC"){
    deconv_res = RPC(input_methyl, compTable, probes_select)
  }
  else if (deconv_algorithm == "CBS"){
    deconv_res = CBS(input_methyl, compTable, probes_select)
  }
  else if (deconv_algorithm == "MethylResolver"){
    deconv_res = MethylResolver(input_methyl, compTable)
  }
  else if (deconv_algorithm == "MCP-counter"){
    deconv_res = enrichment_score(input_methyl, ref_betamatrix, ref_phenotype, method = "MCP-counter")
  }
  else if (deconv_algorithm == "ssGSEA"){
    deconv_res = enrichment_score(input_methyl, ref_betamatrix, ref_phenotype, method = "ssGSEA")
  }
  else if (deconv_algorithm == "ESTIMATE"){
    deconv_res = enrichment_score(input_methyl, ref_betamatrix, ref_phenotype, method = "ESTIMATE")
  }
  
  return(deconv_res)
}


ref_compTable <- function(ref_betamatrix, ref_phenotype){
  ref_phenotype <- as.factor(ref_phenotype)
  
  ffComp <- genefilter::rowFtests(ref_betamatrix, ref_phenotype)
  prof <- vapply(
    X = splitit(ref_phenotype),
    FUN = function(j) matrixStats::rowMeans2(ref_betamatrix, cols = j),
    FUN.VALUE = numeric(nrow(ref_betamatrix)))
  r <- matrixStats::rowRanges(ref_betamatrix)
  
  compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
  names(compTable)[1] <- "Fstat"
  names(compTable)[c(-2, -1, 0) + ncol(compTable)] <- c("low", "high", "range")
  return(compTable)
}


build_reference_450k <- function(extend = TRUE){
  
  library(FlowSorted.Blood.450k)
  CellLines.matrix = NULL
  cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
  ## otherwise all cell types: Bcell, CD4T, CD8T, Eos, Gran, Mono, Neu, NK, WBC, PBMC
  ref_betamatrix <- getBeta(preprocessNoob(FlowSorted.Blood.450k, dyeMethod = "single"))
  ref_phenotype <- as.data.frame(colData(FlowSorted.Blood.450k))$CellType
  keep <- which(ref_phenotype %in% cellTypes)
  ref_betamatrix <- ref_betamatrix[,keep]
  ref_phenotype <- ref_phenotype[keep]
  
  if (extend == FALSE){
    return(list(ref_betamatrix = ref_betamatrix, ref_phenotype = ref_phenotype))
  }
  else{
    library(GEOquery)
    library(minfi)
    rgSet <- read.metharray.exp("GSE40699/idat")
    ref_betamatrix_1 <- getBeta(preprocessNoob(rgSet, dyeMethod = "single"))
    geoMat <- getGEO("GSE40699")
    pD.all <- pData(geoMat[[1]])
    
    ref_betamatrix_EpiFib <- cbind(ref_betamatrix_1[,match(rownames(pD.all)[c(2,12,27,21,28,35,44,46,50,51,56)],substr(colnames(rgSet),1,9))],
                                   ref_betamatrix_1[,match(rownames(pD.all)[c(6,8,10,11,14,16,60)],substr(colnames(rgSet),1,9))])
    ref_phenotype_EpiFib <- c(rep("Epithelial", 11), rep("Fibroblast", 7))
    
    ref_betamatrix <- cbind(ref_betamatrix, ref_betamatrix_EpiFib)
    ref_phenotype <- c(ref_phenotype, ref_phenotype_EpiFib)
    
    return(list(ref_betamatrix = ref_betamatrix, ref_phenotype = ref_phenotype))
  }
}




#'One-vs-All t test feature selection
#'
#'The One-vs-All t test feature selection based on the reference matrix.
#'@param ref_betamatrix The reference matrix ref_betamatrix.
#'@param ref_phenotype The cell type information for the reference matrix.
#'@param probeSelect The selection can be "any" or "both". If "any", the function selects top probes regardless of
#'up-regulation or down-regulation. If "both", half of top probes are picked up from the up-regulated probes while
#'the other half of the top probes are picked up from the down-regulated probes. Default value is "both".
#'@param pv The p-value threshold with default value as 1e-8
#'@param MaxDMRs The number of probes selected with default value as 100.
#'@return A vector of the selected probes.
#'@export


ref_probe_selection_oneVsAllttest <- function(ref_betamatrix, ref_phenotype, probeSelect = "both", pv = 1e-8, MaxDMRs = 100){
  
  ref_phenotype <- as.factor(ref_phenotype)
  
  tIndexes <- splitit(ref_phenotype)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0,ncol(ref_betamatrix))
    x[i] <- 1
    return(genefilter::rowttests(ref_betamatrix, factor(x)))
  })
  
  if (probeSelect == "any") {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[, "p.value"] < pv, ]
      yAny <- y[order(abs(y[, "dm"]), decreasing = TRUE), ]
      c(rownames(yAny)[seq(MaxDMRs)])
    })
  } else {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[, "p.value"] < pv, ]
      yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
      yDown <- y[order(y[, "dm"], decreasing = FALSE), ]
      c(rownames(yUp)[seq_len(MaxDMRs/2)],
        rownames(yDown)[seq_len(MaxDMRs/2)])
    })
  }
  
  trainingProbes <- unique(unlist(probeList))
  return(trainingProbes)
}



#'Build the reference library for 450k arrays
#'
#'Build the reference library for 450k to include "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran",
#'And extend the reference library for 450k to include "Epithelial", "Fibroblast".
#'@param extend If TRUE, generate the extended reference library; otherwise, generate the reference library of
#'only 6 immune cell types. Default value is TRUE.
#'@return A list of beta value reference matrix (ref_betamatrix) and reference cell types (ref_phenotype).
#'@export


  

### (2) perform the analysis
deconv_res <- MethylDeconv_pipeline(beta_filtered, 
                      normalized = TRUE, 
                      array_type = "450k", 
                      feature_selection = "oneVsAllttest",
                      deconv_algorithm = "Houseman",
                      tissue = "general",
                      extend_reference = FALSE,
                      custom_probes = NULL)


## CpGs selected were enriched in gene sets and pathways 
#associated with immune response and cell differentiation, consistent with existing literatures


#Transformar em dataframe
deconv_res <- as.data.frame(deconv_res)

##############################################################################################################


# Check if sample names match
if (all(rownames(deconv_res) == pd_filtered$barcode)) {
  # Merge dataframes based on sample names
  merged_data <- cbind(deconv_res, NSD1status = pd_filtered$NSD1status, AgeAcc = pd_filtered$AgeAcc)
  
  # Reshape data for plotting
  long_data <- reshape2::melt(merged_data, 
                              id.vars = c("NSD1status", "AgeAcc"), 
                              variable.name = "Cell_Type", 
                              value.name = "Proportion")
  
  # Plotting using ggplot2
  library(ggplot2)
  
p1 <- ggplot(long_data, aes(x = Cell_Type, y = Proportion, fill = AgeAcc)) +
    geom_boxplot(position = "dodge") +
    labs(title = "Proportion of Cell Types between Epigenetic Age Accelerated and Deccelerated Groups",
         x = "Cell Type",
         y = "Proportion",
         fill = "Age Acceleration") +
    theme_minimal()
} else {
  warning("Sample names do not match between 'deconv_res' and 'pd_filtered'")
}

library(ggpubr)

ggplot(p)
#  Add p-value
p1 + stat_compare_means()
# Change method
p + stat_compare_means(method = "wilcox.test", label = "p.format")







