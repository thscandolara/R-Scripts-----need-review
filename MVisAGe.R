##   'MVisAGe was designed to help users explore the effects of copy number changes and alterations in DNA methylation status on gene expression. 
##    Gene-level Pearson and Spearman correlation coefficients are computed and plotted according to the genomic position of the underlying genes.'
##    doi:10.1158/0008-5472.CAN-17-3464

# setwd('~/data')

#Loading files
#load("~/tcga_methyl_and_deseq_explore.RData")


#Call libraries
library(MVisAGe)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(biomaRt)

### Annotation data
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

### Creating a new methylation matrix (of beta values in this case) to manipulate
meth_matrix <- beta_filtered

# Transform into df for next steps (I will need to make some edits)
meth_matrix <- as.data.frame(meth_matrix)

# Create a df with probes information
probe.features <- ann450k@listData
probe.features <- as.data.frame(probe.features, row.names = probe.features$Name)

# If 'UCSC_RefGene_Name' is the column that has gene names in 'probe.features':
probe.features$gene <- sapply(strsplit(probe.features$UCSC_RefGene_Name, ";"), function(x) x[1])

# Create a new column in 'meth_matrix' - target-genes of each probe (rownames)
meth_matrix$gene <- probe.features$gene[match(rownames(meth_matrix), rownames(probe.features))]

# Concatenate each pair of gene-probe
new_row_names <- paste0(meth_matrix$gene, "_", rownames(meth_matrix))

# Change rownames in 'meth_matrix'
rownames(meth_matrix) <- new_row_names
head(meth_matrix) ##ok


# Identify which row does not have any gene associated with (they have "NA")
rows_to_exclude <- grep("^NA_", rownames(meth_matrix))

# Delete this rows --> they will not be used for the correlation
meth_matrix <- meth_matrix[-rows_to_exclude, ]
head(meth_matrix) #ok



### ----- Calculate the mean value of methylation for each gene (use the column 'gene' as a variable!!) 
gene_means <- aggregate(. ~ gene, data = meth_matrix, FUN = mean)
#In this step I am calculating the mean value using all probes methylation beta values associated with a unique gene

# Genes with a mean value (unique values) -> rownames
rownames(gene_means) <- gene_means$gene

# Remove column with mean values (the first one, in this case)
gene_means <- gene_means[, -1]

# Remove 'gene' column
gene_means$gene <- NULL  

# Resulting matrix
View(gene_means)
###18592 





### ----- Create a subset with clinical data: keep only variables that I have interest.
# The variable of my interest: epigenetic age acceleration grouping (ACC or DEC)
clinSubset <- pd_filtered[c("patient", "AgeAcc")]
rownames(clinSubset) <- NULL



### ----- Gene expression values should be NORMALIZED values: RPKM, RSEM, etc.,
#1) Loading summarized experiment data downloaded from TCGAbiolinks
load("~/RangedSummarizedExperiment_RNA-seq.rda")
rna <- data

#2) Need to change rownames --> they must be gene names
rownames(rna) <- rowData(rna)$gene_name 

#3) Extract FPKM data
exp_mat <- as.data.frame(SummarizedExperiment::assay(rna, "fpkm_unstrand"))

#4) Clinical data - identifying which sample is normal or tumoral
clinical <- data.frame(rna@colData)
table(clinical$definition)
#Primary solid Tumor      Solid Tissue Normal 
#115                         12

# Substitute spaces for "_" 
clinical$definition <-  gsub(" ", "_", clinical$definition)

# Turn 'definition' column into factor
clinical$definition <- as.factor(clinical$definition)

# Filter clinical data df: keep only samples that also have methylation data
clinical_exp <- subset(clinical, clinical$patient %in% pd_filtered$patient)
#124 samples

# Keep only tumoral samples (exclude normal)
clinical_exp <- clinical_exp[clinical_exp$definition != "Solid_Tissue_Normal", ]
#112

# Filter gene expression matrix based on clinical exp df
exp_mat <- exp_mat[, colnames(exp_mat) %in% rownames(clinical_exp)] 
#112


### ----- Adjust samples names - keep only the first 12  digits
colnames(exp_mat) <- substr(colnames(exp_mat), 1, 12)
colnames(gene_means) <- substr(colnames(gene_means), 1, 12)

# Check if all samples are equal in both matrices
commonSubjects = intersect(intersect(colnames(gene_means), colnames(exp_mat)), clinSubset[,"patient"])
length(commonSubjects)														#112 --OK!


### -----  More adjustments!
meth1 = gene_means[,which(colnames(gene_means) %in% commonSubjects)]
exp1 = exp_mat[,which(colnames(exp_mat) %in% commonSubjects)]
clinSubset1 = clinSubset[which(clinSubset[,"patient"] %in% commonSubjects),]
meth1 = meth1[,order(colnames(meth1))]
exp1 = exp1[,order(colnames(exp1))]
clinSubset1 = clinSubset1[order(clinSubset1[,"patient"]),]
all(colnames(exp1) == colnames(meth1))											#OK
all(colnames(exp1) == clinSubset1[,"patient"])										#OK
exp1 = log2(exp1 + 1) ##transformar em escala de log
hist(exp1, breaks = 50)														#OK



### -----  Create df with both gene and chromosome annotation
# I am using 'biomart' package
                              
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl")
                              
#Which options (attributes) are available
attributes <- listAttributes(ensembl)
listFilters(ensembl)
listAttributes(ensembl)


### -----  Create annotation file
#I also called a few additional columns (which I will not use for MVisage) just to check data
ann_hub <- getBM(attributes = c('hgnc_symbol',
                                'ensembl_gene_id',
                                'chromosome_name', 
                                'start_position', 
                                'band',
                                'transcript_mane_select'),
                                 mart = ensembl)



head(ann_hub)
#86767 obs. of 6 variables


# Filter out rows with "MT" in the 'chromosome_name' column (they are mitochondrial genes)
ann_hub <- subset(ann_hub, chromosome_name != "MT")
#86730

#Exclude rows that does not have HGNC info
ann_hub <- ann_hub %>%
  filter(hgnc_symbol != "")
#62809 obs

# Modify chromosome names
ann_hub <- ann_hub %>%
  filter(chromosome_name %in% c(1:22, "X", "Y"))
#56233 obs.


# Exclude genes that does not have a MANE transcript ** this was something which I've performed but it is not indicated in the package...
ann_hub <- ann_hub %>%
  group_by(hgnc_symbol) %>%
  slice(which.max(!is.na(transcript_mane_select)))
#40730 obs.

#Some additional adjustments
ann_hub <- as.data.frame(ann_hub)
rownames(ann_hub) <- ann_hub$hgnc_symbol
ann_hub <- ann_hub[, -4] # !!!! REMINDER: EXCLUDE ENSEMBL ID COLUMN AND TRANSCRIPT MANE COLUMN (they will be different column numbers)! 

# Add "chr" prefix in 'chromosome_name' column
ann_hub <- mutate(ann_hub, chromosome_name = paste0("chr", chromosome_name))

# Substitute X and Y for chr23 and chr24
ann_hub$chromosome_name <- gsub("chrX", "chr23", ann_hub$chromosome_name)
ann_hub$chromosome_name <- gsub("chrY", "chr24", ann_hub$chromosome_name)

## Rename columns -- it must be according to MVisAGE (chr, pos, cytoband) 
ann_hub <- ann_hub %>%
  rename(chr = chromosome_name, pos = start_position, cytoband = band)

# Create another df (If needed, I will have a "backup" file)
geneAnnotData = ann_hub
geneAnnotData

# Before exclude 'ann_hub', check if everything is alright if geneAnnotData.
#rm(ann_hub)

#Remove chrX and chrY from 'geneAnnotData' -- There was none already! I may have excluded them in some step above, need to check it out
geneAnnotData <- geneAnnotData[!(geneAnnotData$chr %in% c("chrX", "chrY")), ]



#Restrict this analysis for genes that are available in all datasets (methylation and gene expression matrices)
commonGenes = intersect(intersect(rownames(exp1), rownames(meth1)), rownames(geneAnnotData))
length(commonGenes)														#15086
exp1 = exp1[commonGenes,]
meth1 = meth1[commonGenes,]
geneAnnotData1 = geneAnnotData[commonGenes,]
all(rownames(exp1) == rownames(meth1))											#OK
all(rownames(exp1) == rownames(geneAnnotData1))										#OK


#Create a file to facilitate input in the formula
prepped.data = data.prep(exp.mat = exp1,
                         cn.mat = meth1,
                         gene.annot = geneAnnotData1,
                         sample.annot = clinSubset1,
                         log.exp = T,#fiz lá em cima, antes de plotar o histograma
                         gene.list = NULL)

pd.exp = prepped.data[["exp"]]
pd.cn = prepped.data[["cn"]]
pd.ga = prepped.data[["gene.annot"]]
pd.sa = prepped.data[["sample.annot"]]

#Perform the correlation analysis following MVisAGE recommendations
output.list = corr.list.compute(pd.exp,
                                pd.cn,
                                pd.ga,
                                pd.sa,
                                method = "spearman",
                                digits = 5,
                                alternative = "greater")
names(output.list)
#[1] "ACC" "DEC"

head(output.list[["DEC"]])

#Save correlations associated if each group
outputcorr_DEC <- output.list$DEC #15265
outputcorr_ACC <- output.list$ACC #15273



#Visualize
#All chromosomes plot
smooth.genome.plot(plot.list = output.list, 
                   plot.column = "R", 
                   lwd.vec = c(2), 
                   margin.vec = rep(3, 4),
                   expand.size = 50, 
                   loess.span = 100, 
                   ylim.low = -0.3, 
                   ylim.high = -0.1, 
                   legend.loc = "bottomleft", 
                   annot.colors = c("blue", "red"))

#Per region of a chromosome
smooth.region.plot(plot.list = output.list, 
                   plot.column = "R", 
                   lwd.vec = c(2),
                   expand.size = 0, 
                   loess.span = 100, 
                   ylim.low = -0.5, 
                   ylim.high = 0, 
                   legend.loc = "bottomleft", 
                   annot.colors = c("black", "red", "blue"), 
                   plot.chr = 19, 
                   plot.start = 305573, 
                   plot.stop = 58561931,
                   main.label = "Chr19", 
                   xaxis.line = 2)


#plot raw correlation coefficients: best suited to regions containing a relatively small number of genes.
#select a region with just a few genes for a better visualization
unsmooth.region.plot(plot.list = output.list,
                     plot.chr = 19,
                     plot.start = 19515736,
                     plot.stop = 20077994,
                     plot.column = "R",
                     plot.points = T,
                     plot.lines = T,
                     gene.names = T,
                     annot.colors = c("blue", "red"),
                     vert.pad = .05,
                     num.ticks = 10,
                     ylim.low = NULL,
                     ylim.high = NULL,
                     pch.vec = c(16, 16),
                     lty.vec = NULL,
                     lwd.vec = c(1, 1),
                     plot.legend = T,
                     legend.loc = "bottomleft")

prepped.data$cn
class(cn.mat)
dev.off()



cn.region.heatmap(cn.mat= pd.cn,
                  gene.annot = pd.ga,
                  plot.chr = 11,
                  plot.start = 0e6,
                  plot.stop = 135e6,
                  sample.annot = pd.sa,
                  sample.cluster = F,
                  low.thresh = -2,
                  high.thresh = 2,
                  collist = c("blue", "white", "red"),
                  annot.colors = c("blue", "red"),
                  plot.list = output.list,
                  plot.sample.annot = T,
                  cytoband.colors = c("gray90", "gray60"))





