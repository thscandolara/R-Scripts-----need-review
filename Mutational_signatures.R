#### -----------------------   This code is under construction!   ---------------------------####



##Genomic data analysis ---> Downloaded with TCGAbiolinks code (if needed, check it out)


#Call libraries
library(dplyr)
library(maftools)


# If 'clinical_exp' is my dataframe containing clinical information regarding samples (this df is the result of several other codes, but in summary it has only
# tumor samples that have both methylation and gene expression data, n=112!), and 'maf'is the file downloaded via TCGAbiolinks:

# Extract first 16 digits from the column 'Tumor_Sample_Barcode' in  'maf' -- I want to keep only TUMOR SAMPLES that also have methylation and gene expression data
maf$Tumor_Sample <- substr(maf$Tumor_Sample_Barcode, 1, 16)

# Filter 'maf' - I want to keep matching rows, whereas 'Tumor_Sample' is also found in 'clinical_exp$sample'
maf_filtrado <- maf[maf$Tumor_Sample %in% clinical_exp$sample, ]

maf_filtrado$Tumor_Sample_Barcode <- maf_filtrado$Tumor_Sample

#create another clinical data df that I can manipulate 
clinical_data <- clinical_exp
clinical_data$Tumor_Sample_Barcode <- clinical_data$sample



##Turn 'maf' into a MAF file with maftools
maf_filtrado <- read.maf(maf_filtrado, 
                         clinical_data)
#106 samples with alterations (total n was 112)
#-Validating
#-Silent variants: 7069 
#-Summarizing
#--Possible FLAGS among top ten genes:
#TTN
#SYNE1
#MUC16
#USH2A
#-Processing clinical data
#-Finished in 5.347s elapsed (16.2s cpu) 

#print
maf_filtrado

getSampleSummary(maf_filtrado)
getGeneSummary(maf_filtrado)

#save
write.mafSummary(maf = maf_filtrado, 
                 basename = 'maf_larynx_tcga')


#Visualize results
plotmafSummary(maf = maf_filtrado, 
               rmOutlier = TRUE, 
               addStat = 'median', 
               dashboard = TRUE, 
               titvRaw = FALSE)

oncoplot(maf = maf_filtrado, top = 10)

rainfallPlot(maf = maf_filtrado, 
             detectChangePoints = TRUE, 
             pointSize = 0.4)



#exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = maf_filtrado, 
                    top = 25, 
                    pvalue = c(0.05, 0.1))

#Clinical enrichment based on epigenetic age acceleration column = checking if there is any difference between accelerated and deaccelerated groups
ageacc.ce = clinicalEnrichment(maf = maf_filtrado, 
                               clinicalFeature = 'AgeAcc')

ageacc.ce$groupwise_comparision[p_value < 0.05]

plotEnrichmentResults(enrich_res = ageacc.ce, 
                      pVal = 0.05, 
                      geneFontSize = 0.5, 
                      annoFontSize = 0.6)


#Call libraries for mutational signatures analysis
library(BSgenome.Hsapiens.UCSC.hg38)

# prep hg18 
hg18 <- BSgenome.Hsapiens.UCSC.hg38


#Step 1)
x <- trinucleotideMatrix(maf = maf_filtrado,
                         ref_genome = 'BSgenome.Hsapiens.UCSC.hg38',
                         prefix = NULL,
                         add = TRUE,
                         ignoreChr = NULL,
                         useSyn = TRUE,
                        fn = NULL)


plotApobecDiff(tnm = x, maf = maf_filtrado, pVal = 0.2)




#Step2) 
library("NMF")
x.sign <- estimateSignatures(mat = x, 
                              plotBestFitRes = FALSE, 
                              nMin = 2, 
                              nTry = 6, 
                              nrun = 5, 
                              pConstant = 0.01)

#Best possible signature is the value at which Cophenetic correlation drops significantly.
plotCophenetic(x.sign, bestFit = NULL)
##best drop was in 5 here.. 

x.ext.sign <- extractSignatures(mat = x, n = 5)



##compareSignatures returns full table of cosine similarities against COSMIC signatures, which can be further analysed. 

#Compare against original 30 signatures 
x.ext.sign.og30.cosm = compareSignatures(nmfRes = x.ext.sign, 
                                   sig_db = "legacy")

#-Comparing against COSMIC signatures
#------------------------------------
#  --Found Signature_1 most similar to COSMIC_13
#Aetiology: APOBEC Cytidine Deaminase (C>G) [cosine-similarity: 0.851]
#--Found Signature_2 most similar to COSMIC_5
#Aetiology: Unknown [cosine-similarity: 0.803]
#--Found Signature_3 most similar to COSMIC_6
#Aetiology: defective DNA mismatch repair [cosine-similarity: 0.856]
#--Found Signature_4 most similar to COSMIC_4
#Aetiology: exposure to tobacco (smoking) mutagens [cosine-similarity: 0.952]
#--Found Signature_5 most similar to COSMIC_19
#Aetiology: Unknown [cosine-similarity: 0.726]


#Compate against updated version3 60 signatures 
x.ext.sign.v3.cosm = compareSignatures(nmfRes = x.ext.sign, sig_db = "SBS")

#-Comparing against COSMIC signatures
#------------------------------------
#  --Found Signature_1 most similar to SBS13
#Aetiology: APOBEC Cytidine Deaminase (C>G) [cosine-similarity: 0.718]
#--Found Signature_2 most similar to SBS5
#Aetiology: Unknown [cosine-similarity: 0.783]
#--Found Signature_3 most similar to SBS6
#Aetiology: defective DNA mismatch repair [cosine-similarity: 0.835]
#--Found Signature_4 most similar to SBS4
#Aetiology: exposure to tobacco (smoking) mutagens [cosine-similarity: 0.936]
#--Found Signature_5 most similar to SBS5
#Aetiology: Unknown [cosine-similarity: 0.674]


##Below plot shows comparison of similarities of detected signatures against validated signatures.
library('pheatmap')
pheatmap::pheatmap(mat = x.ext.sign.og30.cosm$cosine_similarities, 
                   cluster_rows = FALSE, 
                   main = "cosine similarity against validated signatures")



maftools::plotSignatures(nmfRes = x.ext.sign,
                         contributions = F,
                         title_size = 1.2, 
                         sig_db = "SBS")




# Trying a different approach to check these signatures
library('sigminer')

#!!! This took a while:
nmf_matrix_sigest <- sig_estimate(nmf_matrix = x$nmf_matrix,
             range = 2:5,
             nrun = 10,
             use_random = FALSE,
             method = "brunet",
             seed = 123456,
             cores = 1,
             keep_nmfObj = FALSE,
             save_plots = FALSE,
             plot_basename = file.path(tempdir(), "nmf"),
             what = "all",
             verbose = FALSE)


# In this one I think the best drop in cophonetic was with 3 signatures... 
show_sig_number_survey2(nmf_matrix_sigest$survey,
  y = NULL,
  what = c("all", "cophenetic", "rss", "residuals", "dispersion", "evar", "sparseness",
           "sparseness.basis", "sparseness.coef", "silhouette", "silhouette.coef",
           "silhouette.basis", "silhouette.consensus"),
  na.rm = FALSE,
  xlab = "Total signatures",
  ylab = "",
  main = "Signature number survey using NMF package"
)


