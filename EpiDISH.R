## Script to perform deconvolution analysis using EpiDISH package.


#setwd("~/data")
#____Load methylation and gene expression data____#
load("TCGA_files.RData")
# 1. Only tumor samples; 
# 2. There is information regarding epigenetic age acceleration and NSD1 status.

###___Files loaded into my environment:
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



#load libraries that we will use
library(EpiDISH)
library(ggplot2)
library(ggpubr)


#load the reference 
data(centEpiFibIC.m)

out.l <- epidish(beta.m = beta_filtered, 
                 ref.m = centEpiFibIC.m, 
                 method = "RPC") 
# estF is the matrix of estimated cell-type fractions. 
# ref is the reference centroid matrix used
# dataREF is the subset of the input data matrix over the probes defined in the reference matrix.




out.l$estF
dim(out.l$ref)
#[1] 489   3
dim(out.l$dataREF)
#[1] 489 114

#probes actually used to build the model and infer the cell-type fractions


## How to estimate cell-type fractions in the two-step framework.
data(centBloodSub.m)
frac.m <- hepidish(beta.m = beta_filtered, 
                   ref1.m = centEpiFibIC.m, 
                   ref2.m = centBloodSub.m, 
                   h.CT.idx = 3, 
                   method = 'RPC')
frac.m
#                                 Epi         Fib           B         NK        Mono
#TCGA-BA-4076-01A-01D-1433-05 0.8369611 0.012991022 0.000000000 0.04652298 0.103524895
#TCGA-BA-4078-01A-01D-1433-05 0.5425355 0.218326317 0.052285679 0.14465572 0.042196780
#TCGA-BA-5555-01A-01D-1511-05 0.4415564 0.194920941 0.078277156 0.18306043 0.102185072

boxplot(frac.m)

frac.m <- as.data.frame(frac.m)

merged_data$Tobacco


# Check if sample names match
if (all(rownames(frac.m) == pd_filtered$barcode)) {
  # Merge dataframes based on sample names
  merged_data <- cbind(frac.m, NSD1status = pd_filtered$NSD1status, 
                       AgeAcc = pd_filtered$AgeAcc, 
                       Tobacco = pd_filtered$cigarettes_per_day)
  

#Adjust tobacco variable
  merged_data <- merged_data %>%
    mutate(Tobacco = case_when(
      is.na(Tobacco) ~ "Non-smoker",
      is.numeric(Tobacco) ~ "Smoker",
      TRUE ~ as.character(Tobacco)  # Preserve non-numeric and non-NA values
    ))
  
# Reshape data for plotting
  long_data <- reshape2::melt(merged_data, 
                              id.vars = c("NSD1status", "AgeAcc", "Tobacco"), 
                              variable.name = "Cell_Type", 
                              value.name = "Proportion")
  
  
  
  # Plotting using ggplot2

  
  p1 <- ggplot(long_data, aes(x = Cell_Type, y = Proportion, fill = Tobacco)) +
    geom_boxplot(position = "dodge") +
    labs(title = "Proportion of Cell Types between NSD1 status of samples",
         x = "Cell Type",
         y = "Proportion",
         fill = "Tobacco") +
    theme_minimal()
} else {
  warning("Sample names do not match between 'deconv_res' and 'pd_filtered'")
}

# Change method
p1 + stat_compare_means(method = "wilcox.test", label = "p.format")

##0k
