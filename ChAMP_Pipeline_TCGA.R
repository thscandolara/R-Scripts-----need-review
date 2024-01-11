##-----------ChAMP pipeline for methylation data analysis!

###-----------------------------------------------------------------------------###
#DOC= METHYLATION AND CORRELATION ANALYSIS -- TCGA DATASET
#AUTHOR="Thalita Basso Scandolara"
###-----------------------------------------------------------------------------###
##Call libraries
library(ChAMP)
library(minfi)
library(limma)
library(ggplot2)
library(RColorBrewer)
library(edgeR)
library(survival)
library(survminer)
library(ggplot2)


#1) Set working directory
setwd("/data00/PRJ_SCRATCH/SheilaCoelho/headspace_methyl/data/")


#2) Load methylation data from TCGA (TCGAquery e TCGAprepare):
load("~/RangesSummarizedExperiment_TCGA_Methylation.rda")
dna_methyl_full <- data
rm(data)


#3) Create a beta-values dataframe
beta <- as.data.frame(SummarizedExperiment::assay(dna_methyl_full))


#4) Create a phenotype dataframe (pd)
clinical <- data.frame(dna_methyl_full@colData)


#5) Exclude samples that are not larynx (wrongly classified in TCGA - basaloid and tonsil carcinomas) from phenotype data (pd)
clinical <- clinical[!(rownames(clinical) %in% c("TCGA-CR-7404-01A-11D-2130-05", 
                                                 "TCGA-KU-A66S-01A-21D-A30F-05", 
                                                 "TCGA-QK-AA3J-01A-11D-A392-05")), ]

#6) Filtering beta matrix according to rownames included in 'clinical'
beta <- beta[, rownames(clinical)]

#6.1) Finding probes with NA values 
probe.na <- rowSums(is.na(beta))

#6.2) Excluding these NA probes
probe <- probe.na[probe.na == 0]
beta <- beta[row.names(beta) %in% names(probe), ]


# 7) Transforming my beta-vaue df into a matrix (needed for champ.filter)
beta_matrix <- as.matrix(beta)

#---------------------------------------------------------------------------------------------------------------------


#                                  S T A R T    C h A M P     P I P E L I N E


#champ.filter function: In this step, these probes can be filtered out: NoCG, SNPs start, MultiHit start, XY start..
met <- champ.filter(beta = beta_matrix,
                      pd = clinical,
                      filterXY = T,
                      filterNoCG = T,
                      filterSNPs = T,
                      filterMultiHit = T,
                      population = NULL,
                      filterDetP = TRUE)
##[ Section 2: Filtering Start >>
#Filtering NoCG Start
#Only Keep CpGs, removing 1373 probes from the analysis.
#
#[ Section 2: Filtering Start >>
#Filtering NoCG Start
#Only Keep CpGs, removing 1375 probes from the analysis.
#
#Filtering SNPs Start
#Using general 450K SNP list for filtering.
#Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
#    Removing 881 probes from the analysis.
#
#  Filtering MultiHit Start
#    Filtering probes that align to multiple locations as identified in Nordlund et al
#    Removing 10 probes from the analysis.
#
#  Filtering XY Start
#    Filtering probes located on X,Y chromosome, removing 6831 probes from the analysis.
#
#  Updating PD file
#    filterDetP parameter is FALSE, so no Sample Would be removed.
#
#  Fixing Outliers Start
#    Replacing all value smaller/equal to 0 with smallest positive value.
#    Replacing all value greater/equal to 1 with largest value below 1..
#[ Section 2: Filtering Done ]
#
# All filterings are Done, now you have 345347 probes and 130 samples.


# Saving results: no need to do more than once
#save (met, beta_matrix, clinical, file = "champ_filter_tcga.RData")



#----------------------------------------------------------------------------------------------------------------

#                   F I N D I N G    D I F F E R E N T I A L L Y    M E T H Y L A T E D     P R O B E S


#load ("champ_filter_tcga.RData")

## DMP analysis: We will extract probes that are differentially methylated between normal and tumor samples
diff.met <- champ.DMP(beta = met$beta, 
                      pheno = met$pd$shortLetterCode, 
                      adjPVal = 0.001,
                      adjust.method = "BH")
#Your pheno information contains following groups. >>
#<TP>:114 samples.
#<NT>:16 samples.


#  Final result:
#You have found 56054 significant MVPs with a BH adjusted P-value below 0.001.


# !!! ATTENTION:
#  We also want to identify probes that are at least >=|0.3| differently methylated. Thus:

#Filtering the 'diff.met' file to keep only probes with deltaBeta values of >=|0.3|
deltabeta = with(diff.met$TP_to_NT, which(abs(diff.met$TP_to_NT$deltaBeta) >= 0.3))
diff.met.db = diff.met$TP_to_NT[deltabeta,]

#This resulted in a total of 3893 probes.



##Saving the DMP dataframe with deltaBeta values >=|0.3| as a .csv file!
write.csv(diff.met.db, file = "DMPs_tcga.csv")


#Before continue, we may check if both files rownames are still the same (have I lost any sample?)
a <-rownames(diff.met.db)
b <- rownames(met$beta)
a %in% b 
#[TRUE]

#Whick values of 'a' are not in 'b'?
a[!a %in% b]
#character(0); 
#This means that all 'a' values are in 'b'! I can keep going.

#Exclude this checking files
rm(a, b)
#-----------------------------------------------------------------------------------------------------------------------


#        D M P.G U I       G R A P H I C A L     V I S U A L I Z A T I O N      O F      D M P ' s
                                                   

#DMP.GUI(DMP=diff.met.db,
#        beta=met$beta,
#        pheno=met$pd$shortLetterCode)




#######################################################################################################################################################


#                             I N T E G R A T E D    A N A L Y S I S:

#                 M E T H Y L A T I O N       A N D       G E N E    E X P R E S S I O N



## As we want to perform a correlation analysis between gene expression and methylation levels, 
# we need to have RNA-Seq data as well. 

#1) First, we need to load RNA-seq data that we've download from TCGAquery and TCGAprepare
load("~/RangesSummarizedExperiment_RNA-seq.rda")
dna_exp_full <- data
rm(data)

#2) I need to change the rownames to be the gene name - it will be necessary for the correlation analysis
rownames(dna_exp_full) <- rowData(dna_exp_full)$gene_name 

#3) Create a gene-expression value matrix. 
#In this step I can change my choice ('fpkm_unstrand') for other data types available 
rna <- as.data.frame(SummarizedExperiment::assay(dna_exp_full, "fpkm_unstrand"))


#4) Saving the file
#saveRDS (rna, file = "rna.RDS", compress = FALSE)

##Loading
#rna <- readRDS("rna.RDS")


# Create a clinical data df based on rna-seq dataset. 
clinical.exp <- data.frame(dna_exp_full@colData)

# !!!! I M P O R T A N T
# I need to perform the following steps to do the integrated analysis 

# How many normal and tumor samples do we have? 
table(clinical.exp$shortLetterCode)
#NT  TP 
#12 115


# As we are not interested in normal samples gene expression, we will exclude all 'NT' samples from the dataset

#Removing normal samples from the clinical df
clinical.exp <- clinical.exp[clinical.exp$shortLetterCode == "TP", ]


# Keep only common samples between 'clinical.exp' df and 'rna' matrix 
rna <- rna[, row.names(clinical.exp)]
all(rownames(clinical.exp) %in% colnames(rna))
str(rna)

rna <- as.matrix(rna)
str(rna)
#num [1:60660, 1:115]
#----------------------------------------------------------------------------------------------------------------



#                   F I L T E R I N G    M E T H Y L A T I O N    A N D    G E N E   E X P.   M A T R I X



# Create a filtered 'met' file - we want to keep only those DMPs with deltabeta >=|0.3|: we will use the 'diff.met.db' file

#Extract probes names from the 'diff.met.db' file
probes_diff.met <- rownames(diff.met.db)

# Filter the 'met$beta' matrix using the probes names contained on 'probes_diff.met'
met_f <- met$beta[rownames(met$beta) %in% probes_diff.met, ]
#num [1:3822, 1:130]


##Checking if rownames are identical
a <- rownames(met_f)
b <- rownames(met$beta)

a %in% b 
#[TRUE]

#Is there really no element present in 'a' that is not present in 'b'?
a[!a %in% b]
#character(0)
rm(a,b)


# !!!! IMPORTANT: in this 'met_f' file I still have NORMAL SAMPLES (not filtered yet!)

# We want to keep in 'met_f' matrix ONLY TUMOR SAMPLES. To do this, we need to:

# Remove normal samples 'TP' from the 'clinical' dataframe
clinical.met <- clinical[clinical$shortLetterCode == "TP", ]

# Filter the 'met_f' matrix to keep only samples that are present in the df 'clinical.met';
# In this df,  rownames are the samples names, but in 'met_f', samples names are the colnames!
clinical.met_rownames <- rownames(clinical.met)
met_f <- met_f[, colnames(met_f) %in% clinical.met_rownames]


# Checking if we did it right 
all(colnames(met_f) %in% rownames(clinical.met))#TRUE
all(rownames(clinical.met) %in% colnames(met_f))#TRUE
dim(met_f)
#[1] 3822  114

# yay. 
# Now, we have a matrix 'met_f' that have probes with differently methylation and deltabeta values of >=|0.3|. 
# However, 'met_f' includes only tumor samples!

# Keep going:


# We need to have the same sample's name in both matrices. These samples have the first 19 digits identical. 
colnames(met_f) <- substr(colnames(met_f),1,19) 
colnames(rna) <- substr(colnames(rna),1,19)


# We need to ensure that both 'rna' and 'met_f' df have the same set of columns (thus, samples).

#1) We first retrieve column names of each df, and find each element correspondent in 'rna' that is also present in 'met_f'.
rna <- rna[, colnames(rna) %in% colnames(met_f)] 

#2) Subset the 'met_f' df by selecting only the columns whose names match the column names of the 'rna' df. 
#  The result is stored back in the 'met_f' df.
met_f <- met_f[, colnames(rna)]

#3) Subset the 'rna' df as well. We will select only the matching column names of the 'met_f' df already filtered. 
#The result is stored back in the 'rna' df.
rna <- rna[, colnames(met_f)]


#Checking if the matrices colnames dimensions are identical:
dim(met_f)#[1] 3822  112
dim(rna)#[1] 60660   112

# Okay! I have 112 columns in both matrices. I can keep going.
# As we can see, rownames are different: we have just 3822 probes, but 60660 genes.. we will have adjustments for it





#-----------------------------------------------------------------------------------------------------------------------


#                            C R E A T I N G    T H E     DB.GENES    F I L E

# From the file ' diff.met.db', I will create another df, named 'db.genes'
# db.genes will have a column with the probes names ('name'), and their gene-target in a column named 'gene'

#First, we need to transform 'met_f' into a dataframe if it is as matrix.
met_f <- as.data.frame(met_f)

#Then, we create the db.genes df.
db.genes_f <- data.frame(Column1= rownames(diff.met.db), 
                         Column2 = diff.met.db$gene)

#Column 1 will be my 'name' column
#Column 2 will be my 'gene' column
colnames(db.genes_f) <- c("name", "gene") 


#Print
db.genes_f


# !!!! I M P O R T A N T:
# For the correlation analysis, we need to transform df into a TIBBLE.

db.genes_f <- tidyr::separate_rows(db.genes_f, name, gene)
str(db.genes_f)
#tibble [3,861 × 2] 


#Then, we combine both columns to check if there is duplicates:
db.genes_f$comb <- paste(db.genes_f$name,db.genes_f$gene)


# This keeps only the first occurrence of each unique combination --- for sure not the best option...
db.genes_f <- db.genes_f[!duplicated(db.genes_f$comb), ] 
str(db.genes_f)
#tibble [3,861 × 3] 

#Removing the 'comb' column
db.genes_f <- db.genes_f[, -3]
str(db.genes_f)
#tibble [3,861 × 2] 


## As we have probes not associated with any gene, we remove this probes from the following analysis
db.genes_f <- db.genes_f[trimws(db.genes_f$gene) != "", ]
db.genes_f
## Our tibble remained with 2,511 obs.


## Checking if the rownames of our gene expression matrix could be duplicated as well:

#First, we 'save' our rownames in a new df
df_t <- data.frame(name = row.names(rna)) 
##60660 obs.

#Then, we split every string of the saved 'df_t' 'name' column
df_t <- data.frame(do.call('rbind', strsplit(as.character(df_t$name),'|',fixed=TRUE))) 

# Search for duplicates in rowname, if any
rowName <- df_t$do.call..rbind...strsplit.as.character.df_t.name........fixed...TRUE..
#Plot
table(duplicated(rowName))
#FALSE 
# Thus, we don't have any duplicated genes 
row.names(rna) <- rowName

#We can remove datasets that we are not using again
rm(df_t, rowName) 

# Notice that we generated a table equal to the 'norm.count'
# The only difference is that we have a smaller barcode and kept only tumor samples



#-------------------------------------------------------------------------------------------------------------------------

#                         P E R F O R M     T H E     C O R R E L A T I O N      A N A L Y S I S  
#Notice: I am following the approach proposed in this: https://www.biostars.org/p/465232/ 


#Create an empty dataframe to stock generated data
cis.reg = data.frame( gene=character(0), cpg=character(0), pval=numeric(0), cor=numeric(0))

#Function for the correlation analysis.
for (i in 1:nrow(db.genes_f)){
  cpg <- db.genes_f[i,][1] 
  gene <- db.genes_f[i,][2]
  if (gene %in% rownames(rna)){ 
    df1 <- data.frame(exp= rna[as.character(gene), ])
    df2 <- t(met_f[as.character(cpg), ])
    df <- merge(df1,df2, by = 0) ## This means that df rows will be combined when the first column values are identical
    res <- cor.test(df[,2], df[,3], method = "spearman")# we could use other method if wanted
    pval = round(res$p.value, 20)
    cor = round(res$estimate, 4)
    cis.reg[i,] <- c(gene, cpg, pval, cor)
  }
}


cis.reg$adj.P.Val = round(p.adjust(cis.reg$pval, "fdr"),20)
cis.reg <- cis.reg[with(cis.reg, order(cor, adj.P.Val)), ]
#2511 obs. -- this is the number of probes contained in 'db.genes_f' df. Thus, we ran a corr analysis for all cpg's.

# Saving this result. 
#write.csv(cis.reg, file="cis.reg_tcga.csv")



# We want to know WHERE this probes are located, thus:

#Checking in the df 'probe.features', probe characteristics..
corr_probes <- probe.features[rownames(probe.features) %in% cis.reg$cpg, ]
#We remained with 2109. Still trying to understand why we lost some probes..

#Saving:
write.csv(corr_probes, file = "corr.probes_tcga.csv")


# Now, we want to subset our beta-value matrix to keep only probes that presented correlation with gene expression
# We decided to keep only those probes with moderate to strong correlations (>= |0.5|)

# We need to filter 'cis.reg' df to find these probes 
filtered_corr_probes <- cis.reg[abs(cis.reg$cor) >= 0.5, ] 


# We also have to remove from this new df all rows with 'NA' values.. probably those probes that were also excluded in 
# the filtering process with 'probe.features' df.
filtered_corr_probes <- filtered_corr_probes[complete.cases(filtered_corr_probes), ] 


#We also want to keep only significant correlations!
filtered_corr_probes <- filtered_corr_probes[filtered_corr_probes$adj.P.Val < 0.05, ]
#Remained with 112 probes.

#Saving this file
#write.csv(filtered_corr_probes, "corr.probes.significant_tcga.csv")




#-----------------------------------------------------------------------------------------------------------------------
##--------------------------------------------Additional analysis (just out of curiosity):

## Correlation analysis with gene expression normalized counts vs DMPs (samples with NSD1 mutated)



#                            C R E A T I N G    T H E     df_corr    F I L E
# From the file ' dmpNSD1_f', I will create another df, named 'db.genes'
# db.genes will have a column with the probes names ('name'), and their gene-target in a column named 'gene'

#First, we need to transform 'met_f' into a dataframe if it is as matrix.
sig_DMPs_corr <- as.data.frame(sig_DMPs_corr)

#Then, we create the db.genes df.
df_corr <- data.frame(Column1= rownames(dmpNSD1_f), 
                      Column2 = dmpNSD1_f$gene)

#Column 1 will be my 'name' column
#Column 2 will be my 'gene' column
colnames(df_corr) <- c("name", "gene") 


#Print
df_corr


# !!!! I M P O R T A N T:
# For the correlation analysis, we need to transform df into a TIBBLE.

df_corr <- tidyr::separate_rows(df_corr, name, gene)
str(df_corr)

#tibble [5,298 × 2] 


#Then, we combine both columns to check if there is duplicates:
df_corr$comb <- paste(df_corr$name,df_corr$gene)

## As we have probes not associated with any gene, we remove this probes from the following analysis
df_corr <- df_corr[trimws(df_corr$gene) != "", ]
df_corr
## Our tibble remained with 2,830 obs.
df_corr <- df_corr[!duplicated(df_corr$comb), ]
df_corr <- df_corr[, -3]

#matrices must be in dataframe form. 




#editing expression matrix rowname
df <- data.frame(name = row.names(vsd_corr)) # keeping rownames as a temporary data frame
df <- data.frame(do.call('rbind', strsplit(as.character(df$name),'|',fixed=TRUE))) 
rowName <- df$X1
# find duplicates in rowName, if any
table(duplicated(rowName))


rm(df, rowName) # removing datasets that we do not need anymore


#-------------------------------------------------------------------------------------------------------------------------
#                         P E R F O R M     T H E     C O R R E L A T I O N      A N A L Y S I S  

#Create an empty dataframe to stock generated data
corRes_DMPs_DDS = data.frame( gene=character(0), cpg=character(0), pval=numeric(0), cor=numeric(0))

#Function for the correlation analysis.
for (i in 1:nrow(df_corr)){
  cpg <- df_corr[i,][1] 
  gene <- df_corr[i,][2]
  if (gene %in% rownames(vsd_corr)){ 
    df1 <- data.frame(exp= vsd_corr[as.character(gene), ])
    df2 <- t(sig_DMPs_corr[as.character(cpg), ])
    df <- merge(df1,df2, by = 0) 
    ## This means that df rows will be combined when the first column values are identical
    res <- cor.test(df[,2], df[,3], method = "spearman")# we could use other method if wanted
    pval = round(res$p.value, 20)
    cor = round(res$estimate, 4)
    corRes_DMPs_DDS[i,] <- c(gene, cpg, pval, cor)
  }
}


corRes_DMPs_DDS$adj.P.Val = round(p.adjust(corRes_DMPs_DDS$pval, "fdr"),20)
corRes_DMPs_DDS <- corRes_DMPs_DDS[with(corRes_DMPs_DDS, order(cor, adj.P.Val)), ]
#2830 obs. -- this is the number of probes contained in 'df_corr' df. Thus, we ran a corr analysis for all cpg's.



# Now, we want to subset our beta-value matrix to keep only probes that presented correlation with gene expression
# We decided to keep only those probes with moderate to strong correlations (>= |0.5|)
# We need to filter 'cis.reg' df to find these probes 
corr_DEG_DMP <- corRes_DMPs_DDS[abs(corRes_DMPs_DDS$cor) >= 0.5, ] 
#2106 obs.

# We also have to remove from this new df all rows with 'NA' values.. probably those probes that were also excluded in 
# the filtering process with 'probe.features' df.
corr_DEG_DMP <- corr_DEG_DMP[complete.cases(corr_DEG_DMP), ] 
#288


#We also want to keep only significant correlations!
corr_DEG_DMP <- corr_DEG_DMP[corr_DEG_DMP$adj.P.Val < 0.05, ]

#Saving:
write.csv(corr_DEG_DMP, file = "correlation_DEGs_DMPs_NSD1.csv")



gen.vis <- merge(data.frame(exp= vsd_corr["NOBOX", ]), 
                 t(sig_DMPs_corr[c("cg25317872", "cg11908531", "cg21265917", "cg25317872", "cg21265917", "cg11908531", "cg01427142"), ]),
                 by = 0)

par(mfrow=c(3,2))
sapply(names(gen.vis)[3:8], function(cpg){
  plot(x= gen.vis[ ,cpg], y = gen.vis[,2], xlab = "beta value",
       xlim = c(0,1),
       ylab = "normalized expression" ,
       pch = 19,
       main = paste("NOBOX",cpg, sep = "-"),
       frame = FALSE)
  abline(lm(gen.vis[,2] ~ gen.vis[ ,cpg], data = gen.vis), col = "blue")
})




###--------------------------Starburst plot for trans-regulation visualization
#Call library
library( "IlluminaHumanMethylation450kanno.ilmn12.hg19")
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# adding genes to delta beta data 
tran.reg <- data.frame(ann450k)[rownames(data.frame(ann450k)) %in% dmp_NSD1status$probes, ][, c(4,24)]
tran.reg <- tidyr::separate_rows(tran.reg, Name, UCSC_RefGene_Name) # extending collapsed cells
tran.reg$comb <- paste(tran.reg$Name,tran.reg$UCSC_RefGene_Name) # remove duplicates
tran.reg <- tran.reg[!duplicated(tran.reg$comb), ]
tran.reg <- tran.reg[, -3]
names(tran.reg)[2] <- "GENE"

trans_sig_DMPs_corr <- dmp_NSD1status
# merging with deltabeta dataframe
trans_sig_DMPs_corr$Name <- dmp_NSD1status$probes
# Reorder the columns with "Name" as the first column
trans_sig_DMPs_corr <- trans_sig_DMPs_corr[, c("Name", names(trans_sig_DMPs_corr)[-ncol(trans_sig_DMPs_corr)])]


tran.reg <- merge(tran.reg, trans_sig_DMPs_corr, by = "Name")

# joining with differential expression analysis result
#editing expression matrix rowname

# joining with differential expression analysis result

# merging
tran.reg <- merge(tran.reg, deseq_res_tbl, by = "GENE")
# inspecting data
hist(tran.reg$log2FoldChange)
hist(tran.reg$deltaBeta) 

# defining a column for coloring
tran.reg$group <- ifelse(tran.reg$deltaBeta <= -0.3 & tran.reg$log2FoldChange <= -1.5, "hypo-down",
                         ifelse(tran.reg$deltaBeta <= -0.3 & tran.reg$log2FoldChange >= 1.5, "hypo-up",
                                ifelse(tran.reg$deltaBeta >= 0.3 & tran.reg$log2FoldChange <= -1.5, "hypr-down",
                                       ifelse(tran.reg$deltaBeta >= 0.3 & tran.reg$log2FoldChange >= 1.5, "hypr-up", "not-sig"))))

# plotting
cols <- c("hypo-down" = "#B8860B", "hypo-up" = "blue", "not-sig" = "grey", "hypr-down" = "red", "hypr-up" = "springgreen4")

ggplot(tran.reg, aes(x = deltaBeta, y = log2FoldChange, color = group)) +
  geom_point(size = 2.5, alpha = 1, na.rm = T) +
  scale_colour_manual(values = cols) + 
  theme_bw(base_size = 14) +
  geom_hline(yintercept = 1.5, colour="#990000", linetype="dashed") + 
  geom_hline(yintercept = -1.5, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 0.2, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -0.2, colour="#990000", linetype="dashed") +
  xlab("mean methylation differences") + 
  ylab("Log2 expression change")



