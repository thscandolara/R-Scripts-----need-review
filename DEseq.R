# ---------------------------------
# Data Loading 
# --------------------------------
setwd("~/data")# set your workplace

# --------------------------------------------------------------------------------------------------
# Methylation data loading  --- Input files are all derived from 'ChAMP' and 'methylclock' analyses
# --------------------------------------------------------------------------------------------------
load("~/tcga_tumor_methyl_and_deseq_explore.RData") #dummy name
#Only tumoral data, and with columns describing age acceleration group and mutational status
#Não preciso rodar o data processing abaixo.


## Call libraries
library('ChAMP')
library("DESeq2") 
library("apeglm") 
library("pheatmap") 
library("RColorBrewer") 
library('reshape2')
library('Glimma')
library('edgeR')
library("DEGreport")
library('tidyverse')
library('ggplot2')
library('dplyr')
library('gplots')
library('dplyr')
library('ggtext')




# ---------------------------------------------
# Data processing: if I was to start from zero
# ---------------------------------------------

#Run champ.DMP to get DMPs between groups
dmp_tumor_NSD1 <- champ.DMP(beta = met_filtered$beta, 
                            pheno = met_filtered$pd$NSD1status, 
                            adjPVal = 0.001,
                            adjust.method = "BH")
#85767 obs -- p=0.001 (NSD1 status)

dmp_NSD1status <- as.data.frame(dmp_tumor_NSD1$WT_to_Mut)

#write.csv(dmp_NSD1status, "ALL_DMPs_NSD1.csv")


#convert rownames into a column named 'probes' 
dmp_NSD1status <- tibble::rownames_to_column(dmp_NSD1status, "probes")


DMP.GUI(DMP=dmp_tumor_NSD1[[1]],
        beta=met_filtered$beta,
        pheno=met_filtered$pd$NSD1status)

##Filter to keep probes with a absolute deltabeta value of >= 0.3:
dmpNSD1_f = with(dmp_tumor_NSD1$WT_to_Mut, which(abs(dmp_tumor_NSD1$WT_to_Mut$deltaBeta) >= 0.3))
dmpNSD1_f = dmp_tumor_NSD1$WT_to_Mut[dmpNSD1_f,]
## 5293 probes remained.

dmpNSD1_f$probes <- row.names(dmpNSD1_f)
dmpNSD1_f$probes


#write.csv(dmpNSD1_f, "DMPs_NSD1status_filtered.csv")
#write.csv(levine_probes_info, "Levine.csv")

dmp_tumor_ageacc <- as.data.frame(dmp_tumor_ageacc$DEC_to_ACC)
dmp_tumor_ageacc$probes <- row.names(dmp_tumor_ageacc)



# ---------------------------
# Plotting a heatmap of DMPs 
# ---------------------------
sig_DMPs <- beta_filtered[rownames(beta_filtered) %in% dmpNSD1_f$probes, ]


#Create a df
sample_annotations <- data.frame(
row.names = pd_filtered$barcode,
Age_Acceleration = pd_filtered$AgeAcc, 
Clinical_Stage = pd_filtered$ajcc_clinical_stage,
Nodal_Invasion = pd_filtered$ajcc_pathologic_n,
Tobacco_Exposure = pd_filtered$pack_years_smoked,
Alcohol_history = pd_filtered$alcohol_history,
NSD1_status = pd_filtered$NSD1status)



#Recoding some parameters for better plotting
sample_annotations$Clinical_Stage
sample_annotations$Clinical_Stage <- sample_annotations$Clinical_Stage %>%
  recode("Stage I" = "I",
         "Stage II" = "II",
         "Stage III" = "III",
         "Stage IV"= "IV",
         "Stage IVC" = "IV",
         "Stage IVA" = "IV",
         "Stage IVB" = "IV")

sample_annotations$Clinical_Stage <- coalesce(sample_annotations$Clinical_Stage, "NotReported")
sample_annotations$Clinical_Stage

sample_annotations$Nodal_Invasion <- sample_annotations$Nodal_Invasion %>%
  recode("NX" = "NX",
         "N1" = "N1",
         "N2a" = "N2",
         "N2b" = "N2",
         "N2c" = "N2",
         "N2" = "N2",
         "N3" = "N3",
         "N0" = "N0")

sample_annotations$Nodal_Invasion <- sample_annotations$Nodal_Invasion %>% recode("NA" = "NotReported")
sample_annotations$Nodal_Invasion <- coalesce(sample_annotations$Nodal_Invasion, "NotReported")


sample_annotations$Tobacco_Exposure <- ifelse(!is.na(sample_annotations$Tobacco_Exposure), "Yes", "No")



#Transform into z-score
z_scores <- (sig_DMPs-mean(sig_DMPs))/sd(sig_DMPs)
z_scores


#Plot
pheatmap(z_scores, 
         annotation_col = sample_annotations,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cutree_cols = 2,
         cluster_rows = F,
         cluster_cols = T,
         fontsize_row = 5,
         fontsize_col = 5,
         show_rownames=F)





#----------------
# Volcano plot
#----------------

p <- ggplot(dmp_NSD1status, aes(deltaBeta, -log10(adj.P.Val))) +
  geom_point(aes(color = ifelse(abs(deltaBeta) >= 0.3 & adj.P.Val < 0.001, ifelse(deltaBeta <= -0.3, "Hypomethylated in Mut", "Hypermethylated in Mut"), "Not Significant (Adj. P > 0.001)")), size=1) +
  scale_color_manual(values = c("Hypomethylated in Mut" = "blue", "Hypermethylated in Mut" = "red", "Not Significant (Adj. P > 0.001)" = "grey")) +
  labs(x = "DNA methylation difference (β-values)", y = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), legend.title = element_blank()) +
  geom_vline(xintercept =  c(-0.3, 0.3), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "gray")

# Exibir o gráfico
print(p)



###NOTICE: the further steps are just to save my code. If I were to re-perform this analysis, I could start from deseq files already!
#--------------------------------------------------------------------------
# Loading RNA data: If I was to start from zero again 
#---------------------------------------------------------------------------

#1) First, we need to load RNA-seq data that we've download from TCGAquery and TCGAprepare
load("/RangedSummarizedExperiment_RNA-seq.rda")
dna_exp_full <- data
rm(data)

#2) I need to change the rownames to be the gene name
rownames(dna_exp_full) <- rowData(dna_exp_full)$gene_name 


#3) Create a gene-expression value matrix. 
# exp matrix -- counts
rna <- as.data.frame(SummarizedExperiment::assay(dna_exp_full))

#####Notice: In this step I can change 'fpkm_unstrand' for other data types available 
#rna_fpkm <- as.data.frame(SummarizedExperiment::assay(dna_exp_full, "fpkm_unstrand"))


# clinical data
clinical <- data.frame(dna_exp_full@colData)


# Check how many tumor and control sample are there:
## data are stored under "definition" column in the clinical dataset.
table(clinical$definition)
#Primary solid Tumor      Solid Tissue Normal 
#115                         12




#### I will use the column “definition”, as the grouping variable for gene expression analysis. 

# replace spaces with "_" in levels of definition column
clinical$definition <-  gsub(" ", "_", clinical$definition)

# making the definition column as factor
clinical$definition <- as.factor(clinical$definition)
# releveling factor to ensure tumors would be compared to normal tissue.
#levels(clinical$definition)
#
#clinical$definition <- relevel(clinical$definition, ref = "Solid_Tissue_Normal")


##Filter to keep only tumor samples and samples that we have methylation data
clinical_exp <- subset(clinical, clinical$patient %in% met_filtered$pd$patient)
#124

clinical_exp <- clinical_exp[clinical_exp$definition != "Solid_Tissue_Normal", ]
#112

rm(clinical)

#The count matrix (rna dataset) and the rows of the column data (clinical dataset) MUST be in the same order.

# to see whether all rows of clinical are present in rna datset
all(rownames(clinical_exp) %in% colnames(rna))#TRUE
all(colnames(rna) %in% rownames(clinical_exp))#FALSE


#to keep only colnames in rna present in clinical dataset
rna_exp <- rna[, colnames(rna) %in% rownames(clinical_exp)]

all(colnames(rna_exp) %in% rownames(clinical_exp))#true

# whether they are in the same order:
all(rownames(clinical_exp) == colnames(rna_exp))#true


# if not reorder them by:
#rna <- rna[, rownames(clinical)]
#all(rownames(clinical) == colnames(rna))


#Add info regarding mutational status and ag acceleration in 'clinical_exp'
clinical_exp$AgeAcc <- met_filtered$pd$AgeAcc[match(clinical_exp$patient, met_filtered$pd$patient)]
clinical_exp$NSD1status <- met_filtered$pd$NSD1status[match(clinical_exp$patient, met_filtered$pd$patient)]

# Check if "AgeAcc" column is correctly associated with its corresponding samples
all(clinical_exp$AgeAcc == met_filtered$pd$AgeAcc[match(clinical_exp$patient, met_filtered$pd$patient)])
#TRUE

# Check if "NSD1status" column is correctly associated with its corresponding samples
all(clinical_exp$NSD1status == met_filtered$pd$NSD1status[match(clinical_exp$patient, met_filtered$pd$patient)])
#TRUE


##--- this step is important, must not be avoided!
##Check if samples are in the same ORDER in the clinical data and count matrix
all(colnames(rna_exp) == rownames(clinical_exp))
#[1] TRUE




#----------------------
# Making Deseq object
#----------------------

#The DESeq2 model internally corrects for library size, 
#so transformed or normalized values such as counts scaled by library size 
#should not be used as input: thats why I dont use FPKM as input

# Making DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData = rna_exp,
                              colData = clinical_exp,
                              design = ~ NSD1status) 
#here the algorithm will return the fold change that result from NSD1status 


# I could also use:
#1) design = ~AgeAcc same as above

#2) design = ~AgeAcc + NSD1status; 
#means that deseq2 will test the effect of the NSD1 status (the last factor), 
#controlling for the effect of Age acceleration (the first factor), 
#so that the algorithm returns the fold change result only from the effect of NSD1 status

class(dds)
#class: DESeqDataSet 
#dim: 60660 112 
#metadata(1): version
#assays(1): counts
#rownames(60660): TSPAN6 TNMD ... AL391628.1 AP006621.6
#rowData names(0):
#  colnames(112): TCGA-BA-4076-01A-01R-1436-07 TCGA-BA-4078-01A-01R-1436-07 ... TCGA-UF-A7JJ-01A-11R-A34R-07 TCGA-UF-A7JK-01A-11R-A34R-07
#colData names(80): barcode patient ... AgeAcc NSD1status



#-------------------------------------------------------
# Prefilteration: it is not necessary but recommended  
#-------------------------------------------------------

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds
#dim: 48108 112 




#--------------------------
#  set the factor level 
#-------------------------
# Who is going to be my control parameter? The mutated or the wild-type?
# In this case, my reference as control was the WILD-TYPE group!! 
# So, all data shall be related to the MUTATED samples (in comparison to the wild-type)

dds$NSD1status <- relevel(dds$NSD1status, ref = "WT")
dds$NSD1status



#----------------------
# DE analysis
#----------------------
de_dds <- DESeq(dds)


#--------------------------------------------------

#Save the normalized counts

normalized_counts <- counts(de_dds, normalized=T)
head(normalized_counts)

#--------------------------------------------------




#Recode some columns
de_dds$ajcc_pathologic_stage <- de_dds$ajcc_pathologic_stage %>%
  recode("Stage I" = "Stage I",
         "Stage II" = "Stage II",
         "Stage III" = "Stage III",
         "Stage IV" = "Stage IV",
         "Stage IVA" = "Stage IV",
         "Stage IVB" = "Stage IV",
         "Stage IVC" = "Stage IV")

de_dds$ajcc_pathologic_n <- de_dds$ajcc_pathologic_n %>%
  recode("NX" = "NX",
         "N1" = "LN Invasion",
         "N2a" = "LN Invasion",
         "N2b" = "LN Invasion",
         "N2c" = "LN Invasion",
         "N2" = "LN Invasion",
         "N3" = "LN Invasion",
         "N0" = "No LN invasion")


de_dds$smoke_history <- ifelse(de_dds$cigarettes_per_day > 0, TRUE, FALSE)


## Plot dispersion
plotDispEsts(de_dds)




#----------------------
# DEseq results
#----------------------

resultsNames(de_dds)
#[1] "Intercept"            "NSD1status_Mut_vs_WT"

deseq_res <- results(de_dds, name="NSD1status_Mut_vs_WT")
##The results function of the DESeq2 package performs independent filtering 
#by default using the mean of normalized counts as a filter statistic.



deseq_res
#log2 fold change (MLE): NSD1status Mut vs WT 
#Wald test p-value: NSD1status Mut vs WT 
#DataFrame with 48108 rows and 6 columns
#              baseMean     log2FoldChange   lfcSE         stat         pvalue      padj
#               <numeric>    <numeric>      <numeric>   <numeric>     <numeric>    <numeric>
#TSPAN6     2383.687521      0.1944444      0.1551269    1.253454     0.2100406    0.4407640
#TNMD          0.746763     -0.5842185      0.6400076   -0.912831     0.3613317    0.6047166


##log2FoldChange -- mutated vs wild-type: 
# ALL VALUES DEPICTED IN THESE COLUMN ARE RELATED TO THE MUTATED SAMPLE. 
# EXAMPLE:  a log2foldchange of 0.19 means that mutated samples has logFC 0.19 higher expression of TSPAN6 in comparison to wt (2^0.19 = 1.14, or 14%).



#----------------------
# MA plotting
#----------------------

#The MA plot shows the mean of the normalized counts versus the log2 foldchanges for all genes test
DESeq2::plotMA(deseq_res)
#genes in blue are significantly different!



#----------------------
# Exploring the results
#----------------------

summary(deseq_res)
#out of 48101 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 3491, 7.3%
#LFC < 0 (down)     : 4783, 9.9%
#outliers [1]       : 0, 0%
#low counts [2]     : 11198, 23%
#(mean count < 0)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results



#---------------------------------------------------------
# Log fold change shrinkage for visualization and ranking
#---------------------------------------------------------

#Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. 
#To shrink the LFC, we pass the dds object to the function lfcShrink. 
#Below we specify to use the apeglm method for effect size shrinkage (Zhu, Ibrahim, and Love 2018), 
#which improves on the previous estimator.

res_shrunken <- lfcShrink(de_dds, 
                          coef="NSD1status_Mut_vs_WT", 
                          res = deseq_res,
                          type="apeglm")
#The order of the names determines the direction of fold change that is reported.
#The name provided in the second element is the level that is used as baseline. 
#So for example, if we observe a log2 fold change of -2 this would mean the gene expression is lower 
#in NSD1 mutated relative to the wild-type.

DESeq2::plotMA(res_shrunken)





#-----------------
# Volcano plot
#-----------------

library(EnhancedVolcano)
EnhancedVolcano(res_shrunken,
                lab = rownames(res_shrunken),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'NSD1 Mutant vs Wild-type',
                pCutoff = 10e-16,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 3.0)




#-----------------------------------------------
## If I want to change the pAdj-value:
deseq_res_p0.05 <- results(de_dds, alpha = 0.05)
summary(deseq_res_p0.05)
#------------------------------------------------


# Checking DESEQ results:

#Transform into a tibble
deseq_res_tbl <- as.data.frame(deseq_res)
deseq_res_tbl 
deseq_res_tbl <- tibble::rownames_to_column(deseq_res_tbl, "GENE")

#Inspect the best genes with the lowest p adj
deseq_res_tbl %>%
  arrange(padj) %>%
  head(20)

#      GENE        baseMean     log2FoldChange   lfcSE       stat         pvalue         padj
#1     ERVH48.1  134.476957       5.790840      0.4464934   12.969599   1.819788e-38   6.716839e-34
#2      TUBB8P7   48.021632       3.973419      0.3096922   12.830219   1.110499e-37   2.049426e-33
#3  PROSER2.AS1  124.014482       4.422135      0.3544556   12.475850   1.011216e-35   1.244133e-31



best_genes <- deseq_res_tbl %>%
  arrange(padj) %>%
  head(20) %>% pull(GENE)



metadata(deseq_res)$filterThreshold
#23.27629%  
#0.4742693



#p-values distribution
hist(deseq_res_tbl$padj)







#----------------------------------------------------------------
# Extract most differentially expressed genes due to NSD1 status
#----------------------------------------------------------------
#Log2fold change of <1 and >1 (2 and 0.5) and adj p value of 0.05


#Step 1) filter based on p adj
deseq_res_filtered <- deseq_res_tbl %>% filter(deseq_res_tbl$padj < 0.05)
#6356 genes


#Step 2) filter based on log2fold change
deseq_res_filtered <- deseq_res_filtered %>% filter(abs(deseq_res_filtered$log2FoldChange) >1)
#3643 genes

deseq_res_filtered
#          GENE     baseMean  log2FoldChange     lfcSE        stat       pvalue         padj
#1        MYH16    60.657341       1.124078    0.3033111    3.706023   2.105394e-04   2.755677e-03
#2        WNT16    49.111698       1.566942    0.4091633    3.829624   1.283392e-04   1.834624e-03
#3        ABCB5    10.331169      -1.882327    0.4847566   -3.883036   1.031602e-04   1.531018e-03


#save deseq results -- both the original data and the filtered one
setwd("~/results_deseq/")
write.csv(deseq_res_tbl, 'deseq_results_all.csv')
write.csv(deseq_res_filtered, 'deseq_results_filtered.csv')




#-------------------
# VST transformation
#-------------------
#Perform a variance stabilizing transformation (vst function in deseq)
#The whole point of the VST is to stabilize the features so that there are all 
#on a comparable scale and you haven't inflated the noise in the data

vsd <- vst(de_dds, blind = FALSE)

vsd_corr <- assay(vsd)

##I will use the 'vsd_corr' file 




#---------------------------------
# Visualization - PCA and heatmap
#---------------------------------
#Perform PCA analysis: used to explain variance in gene expression datasets


#Plot PCA
plotPCA(vsd, intgroup = c('NSD1status', 'smoke_history'))
# I can change the variables in 'intgroup'



##Heatmap of sample-to-sample distance matrix (with clustering) based on the normalized counts
sampleDists = dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix)



#Set color scale grade
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)


# Plot
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists, 
         color = colors,
         fontsize = 4)
##what the color indicates: high similarity between samples



#____Contrasts 
# I can make several different comparisons, but it will turn into more complicated models

resultsNames(de_dds)
#[1] "Intercept"            "NSD1status_Mut_vs_WT"

#results(dds, contrast = c("NSD1status", "WT", "Mut")) 
#If I do this, my baseline will be the mutated, not the wild-type anymore



#------------------
# Volcano plot
#------------------

deseq_res_tbl %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_colour_manual(values = c("gray", "red")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  theme(legend.position = "bottom")







#---------------------------------------
# STARBUST PLOT FOR TRANS-REGULATION VISUALIZATION
#-----------------------------------------
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# adding genes to delta beta data 
tran.reg <- data.frame(ann450k)[rownames(data.frame(ann450k)) %in% rownames(dmp_NSD1status), ][, c(4,24)]
tran.reg <- tidyr::separate_rows(tran.reg, Name, UCSC_RefGene_Name) # extending collapsed cells
tran.reg$comb <- paste(tran.reg$Name,tran.reg$UCSC_RefGene_Name) # remove duplicates
tran.reg <- tran.reg[!duplicated(tran.reg$comb), ]
tran.reg <- tran.reg[, -3]
names(tran.reg)[2] <- "GENE"
dmp_nsd1 <- dmp_NSD1status


# merging with deltabeta dataframe
dmp_nsd1$Name <- rownames(dmp_nsd1)
tran.reg <- merge(tran.reg, dmp_nsd1, by = "Name")

# joining with differential expression analysis result
#setting rna deseq results row names 
rownames(deseq_res_filtered) <- deseq_res_filtered[,1]


# merging
tran.reg <- merge(tran.reg, deseq_res_filtered, by = "GENE")
# inspecting data
hist(tran.reg$logFC)
hist(tran.reg$deltaBeta) 


# defining a column for coloring
tran.reg$group <- ifelse(tran.reg$deltaBeta <= -0.2 & tran.reg$log2FoldChange <= -1.5, "hypo-down",
                         ifelse(tran.reg$deltaBeta <= -0.2 & tran.reg$log2FoldChange >= 1.5, "hypo-up",
                                ifelse(tran.reg$deltaBeta >= 0.2 & tran.reg$log2FoldChange <= -1.5, "hypr-down",
                                       ifelse(tran.reg$deltaBeta >= 0.2 & tran.reg$log2FoldChange >= 1.5, "hypr-up", "not-sig"))))

# plotting
cols <- c("hypo-down" = "#B8860B", "hypo-up" = "blue", "not-sig" = "grey", "hypr-down" = "red", "hypr-up" = "springgreen4")

ggplot(tran.reg, aes(x = deltaBeta, y = log2FoldChange, color = group)) +
  geom_point(size = 2, alpha = 1, na.rm = T) +
  scale_colour_manual(values = cols) + 
  theme_bw(base_size = 14) +
  geom_hline(yintercept = 1.5, colour="#990000", linetype="dashed") + 
  geom_hline(yintercept = -1.5, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 0.2, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -0.2, colour="#990000", linetype="dashed") +
  xlab("mean methylation differences") + 
  ylab("Log2 expression change")







#---------------------------------------------
# CORRELATION ANALYSIS BETWEEN DMPs and DEGs
#---------------------------------------------

# Aqui, estou correlacionando as matrizes de valores beta e valores de expressão normalizados
# Como input, estou utilizando apenas as sondas que foram encontradas diferencialmente metiladas
# entre as amostras tumorais NSD1 mutadas vs NSD1 selvagem e os genes encontrados diferencialmente expressos
# entre as amostras NSD1 mutadas vs NSD1 selvagem.



# For TGCA data 1-12 represents the patient and 1-15 represents the sample ID (i.e. primary solid tumor samples)

# Extracting the first 12 digits from vsd_corr and sig_DMPs column names
colnames(vsd_corr) <- substr(colnames(vsd_corr), 1, 12)
colnames(sig_DMPs) <- substr(colnames(sig_DMPs), 1, 12)

# Now sig_DMPs_subset contains only the columns that are present in vsd_corr

# Get the common column names
common <- intersect(colnames(vsd_corr), colnames(sig_DMPs))

# Create sig_DMPs_corr by subsetting sig_DMPs with common columns
sig_DMPs_corr <- sig_DMPs[, common]

rm(common)

#______________preparing methylation data for cis-regulatory analysis____________#

#Create db.genes file based on differentially methylated probes analysis according to NSD1 status;
db.genes <- data.frame(Column1= rownames(dmpNSD1_f), 
                      Column2 = dmpNSD1_f$gene)

#5293 obs. of 2 variables

#Column 1 will be my 'name' column
#Column 2 will be my 'gene' column
colnames(db.genes) <- c("name", "GENE") 

#Print
head(db.genes)
#       name     gene
#1 cg14637927        
#2 cg02973883   DCAF4


#Filter the rna expression matrix (vsd_corr, derived from vst normalization) 
#to keep only genes that were found differentially expressed by deseq analysis
de.genes <- deseq_res_filtered$GENE
cis.exp.mat <- vsd_corr[de.genes, ]
#3645 genes and 112 samples





# !!!! I M P O R T A N T:
# For the correlation analysis, we need to transform dataframe into a TIBBLE.

db.genes <- tidyr::separate_rows(db.genes, name, GENE)
str(db.genes)
#tibble [5,298 × 2] 


#Then, we combine both columns to check if there is duplicates:
db.genes$comb <- paste(db.genes$name,db.genes$GENE)

## As we have probes not associated with any gene, we remove this probes from the following analysis
db.genes <- db.genes[trimws(db.genes$GENE) != "", ]
db.genes
## Our tibble remained with 2,830 obs.

db.genes <- db.genes[!duplicated(db.genes$comb), ]
db.genes <- db.genes[, -3]



sig_DMPs_corr <- as.data.frame(sig_DMPs_corr)



#Checking if the matrices colnames dimensions are identical:
dim(sig_DMPs_corr)#[1] 5293  112
dim(cis.exp.mat)#[1] 3645    112
dim(db.genes)#[1] 2830    2

#-----------------------------------------------------------------------------------------------------------------------

#________ P E R F O R M     T H E     C O R R E L A T I O N      A N A L Y S I S  

cis.reg = data.frame( gene=character(0), cpg=character(0), pval=numeric(0), cor=numeric(0))

#Function for the correlation analysis.
for (i in 1:nrow(db.genes)){
  cpg <- db.genes[i,][1] 
  gene <- db.genes[i,][2]
  if (gene %in% rownames(cis.exp.mat)){ 
    df1 <- data.frame(exp= cis.exp.mat[as.character(gene), ])
    df2 <- t(sig_DMPs_corr[as.character(cpg), ])
    df <- merge(df1,df2, by = 0) ## This means that df rows will be combined when the first column values are identical
    res <- cor.test(df[,2], df[,3], method = "spearman")# we could use other method if wanted
    pval = round(res$p.value, 20)
    cor = round(res$estimate, 4)
    cis.reg[i,] <- c(gene, cpg, pval, cor)
  }
}

cis.reg$adj.P.Val = round(p.adjust(cis.reg$pval, "fdr"),20)
cis.reg <- cis.reg[with(cis.reg, order(cor, adj.P.Val)), ]
#2830 obs. 

# Saving this result. 
setwd("/data00/PRJ_SCRATCH/SheilaCoelho/headspace_methyl/data/results_deseq/")
write.csv(cis.reg, file="tcga_tumorlarynx_cis_reg_correlation_between_DMPs_DEgenes_NSD1status.csv")
setwd("/data00/PRJ_SCRATCH/SheilaCoelho/headspace_methyl/data/")


# We want to know WHERE this probes are located, thus:

probe.features <- ann450k@listData
probe.features <- as.data.frame(probe.features, row.names = probe.features$Name)

#Checking in the df 'probe.features', probe characteristics..
corr_probes_info <- probe.features[rownames(probe.features) %in% cis.reg$cpg, ]
#We remained with 553 Still trying to understand why we lost some probes..

#Saving:
setwd("/data00/PRJ_SCRATCH/SheilaCoelho/headspace_methyl/data/results_deseq/")
write.csv(corr_probes_info, file = "tcga_tumorlarynx_cis_reg_correlation_between_DMPs_DEgenes_NSD1status_probesinfo.csv")
setwd("/data00/PRJ_SCRATCH/SheilaCoelho/headspace_methyl/data/")



# We decided to keep only those probes with moderate to strong correlations (>= |0.5|)
# We need to filter 'cis.reg' df to find these probes 
corr_probes_f <- cis.reg[abs(cis.reg$cor) >= 0.5, ] 
#2411 obs.

# We also have to remove from this new df all rows with 'NA' values
corr_probes_f <- corr_probes_f[complete.cases(corr_probes_f), ] 
#134obs.

corr_probes_f



#----------------------------------------------------------------------
# SAVING FILES FROM GLOBAL ENVIRONMENT!
#--------------------------------------------------------------------

save(beta_filtered, 
     met_filtered, 
     dna_methyl_full, 
     levine_probes_info, 
     pd_filtered,
     dmp_tumor_NSD1,
     dmp_NSD1status,
     dmpNSD1_f,
     sig_DMPs,
     dds, 
     de_dds, 
     normalized_counts,
     deseq_res,
     res_shrunken,
     vsd,
     vsd_corr, 
     file = "tcga_tumor_larynx_methyl_and_deseq_explore_20231214.RData")




########################################################################


# analyses not necessarily included in the DESEQ pipeline 
#** lets say its an extra!!



#----------------------------------------------------------------------
# DMPs for only NSD1wt samples according to Age Acceleration phenotype.
#----------------------------------------------------------------------


#Run champ.DMP to get DMPs between groups
dmp_ageacc <- champ.DMP(beta = met_filtered$beta, 
                         pheno = met_filtered$pd$AgeAcc, 
                         adjPVal = 0.001,
                         adjust.method = "BH")


#Start to Compare : DEC, ACC
#You have found 39969 significant MVPs with a BH adjusted P-value below 0.001.



dmp_ageacc <- as.data.frame(dmp_ageacc$DEC_to_ACC)
#convert rownames into a column named 'probes' 
dmp_ageacc <- tibble::rownames_to_column(dmp_ageacc, "probes")


##Filter to keep probes with a absolute deltabeta value of >= 0.3:
dmp_ageacc_f = subset(dmp_ageacc, abs(dmp_ageacc$deltaBeta) >= 0.3)
#20 probes






#----------------
# VENN DIAGRAM
#----------------

upregulated_mut_genes <- deseq_res_filtered$GENE[deseq_res_filtered$log2FoldChange > 0]
downregulated_mut_genes <- deseq_res_filtered$GENE[deseq_res_filtered$log2FoldChange < 0]

dmpNSD1_dummy <- dmpNSD1_f
dmpNSD1_dummy <- droplevels(dmpNSD1_dummy)
dmpNSD1_dummy$gene[dmpNSD1_f$gene == ""] <- NA
dmpNSD1_dummy <- dmpNSD1_dummy %>% filter(!is.na(gene))


hyper_Mut_genes <- dmpNSD1_dummy$gene[dmpNSD1_dummy$deltaBeta > 0]
hypo_Mut_genes <- dmpNSD1_dummy$gene[dmpNSD1_dummy$deltaBeta < 0]


# Converter as listas  em conjuntos
x <- list(
  upregulated_mut_geneset <- as.character(upregulated_mut_genes),
  downregulated_mut_geneset <- as.character(downregulated_mut_genes),
  hypermethylated_mut_geneset <- as.character(hyper_Mut_genes),
  hypomethylated_mut_geneset <- as.character(hypo_Mut_genes))

# renomear
names(x) <- c(
  "Upregulated Mut",
  "Downregulated Mut",
  "Hypermethylated Mut",
  "Hypomethylated Mut"
)

#check
names(x)

library(VennDiagram)
#Criar o plot
venn.plot <- venn.diagram(x , NULL, fill=c("thistle1", "slategray2","darkolivegreen2","salmon3"), 
                          alpha=c(0.5,0.5, 0.5, 0.5), 
                          cex = 2, 
                          cat.fontface=4, 
                          category.names=c("Upregulated Mut",
                                           "Downregulated Mut",
                                           "Hypermethylated Mut",
                                           "Hypomethylated Mut"), 
                          main="DM and DE genes")


# visualizar
grid.draw(venn.plot)

# To get the list of gene present in each Venn compartment we can use the gplots package
require("gplots")
a <- venn(x, show.plot=T)
# You can inspect the contents of this object with the str() function
str(a)
# By inspecting the structure of the a object created, 
# you notice two attributes: 1) dimnames 2) intersections
# We can store the intersections in a new object named inters
inters <- attr(a,"intersections")

# We can summarize the contents of each venn compartment, as follows:
# in 1) ConditionA only, 2) ConditionB only, 3) ConditionA & ConditionB
lapply(inters, head) 

