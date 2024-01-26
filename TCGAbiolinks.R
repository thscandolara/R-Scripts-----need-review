###-------------------- THIS CODE IS TO BUILD A QUERY FOR TCGABIOLINKS AND DOWNLOAD DATA FROM GDC!
# In this, I have a query for methylation, gene expression and genomic data from TCGA larynx samples

getwd()
setwd()

#Call libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sesame")
library(sesame)


#Check list of available projects
gdc <- getGDCprojects()

# I am interested in TCGA-HNSC
getProjectSummary('TCGA-HNSC')

# --------   BUILDING A QUERY -- RNA-SEQ DATA
query_TCGArna <- GDCquery(project = 'TCGA-HNSC',
                          data.category = 'Transcriptome Profiling')

output_query_TCGArna <- getResults(query_TCGArna)
##All samples from TCGA-HNSC!

## FULL OUTPUT
#write.csv(output_query_TCGArna, file = "outputgeneexp.csv")


# Building a query to retrieve filtered data --> I selected the barcodes of samples I am interested in
query_rna <- GDCquery(project = 'TCGA-HNSC',
                      data.category = 'Transcriptome Profiling',
                      experimental.strategy = 'RNA-Seq',
                      data.type = "Gene Expression Quantification",
                      workflow.type= 'STAR - Counts',
                      access = 'open',
                      barcode = c("TCGA-BA-4076-01A-01R-1436-07",
                                  "TCGA-BA-4078-01A-01R-1436-07",
                                  "TCGA-BA-4078-01A-01T-1435-13",
                                  "TCGA-BA-5555-01A-01R-1513-13",
                                  "TCGA-BA-5555-01A-01R-1514-07",
                                  "TCGA-BA-6868-01B-12R-1914-13",
                                  "TCGA-BA-6868-01B-12R-1915-07",
                                  "TCGA-BA-6869-01A-11R-1872-13",
                                  "TCGA-BA-6869-01A-11R-1873-07",
                                  "TCGA-BA-6870-01A-11R-1872-13",
                                  "TCGA-BA-6870-01A-11R-1873-07",
                                  "TCGA-BA-A6DA-01A-31R-A31N-07",
                                  "TCGA-BA-A6DA-01A-31R-A31R-13",
                                  "TCGA-BA-A6DI-01A-11R-A30B-07",
                                  "TCGA-BA-A6DI-01A-11R-A30I-13",
                                  "TCGA-BB-4217-01A-11R-2080-13",
                                  "TCGA-BB-4217-01A-11R-2081-07",
                                  "TCGA-BB-7862-01A-21R-2231-13",
                                  "TCGA-BB-7862-01A-21R-2232-07",
                                  "TCGA-BB-7864-01A-11R-2231-13",
                                  "TCGA-BB-7864-01A-11R-2232-07",
                                  "TCGA-BB-7870-01A-11R-2231-13",
                                  "TCGA-BB-7870-01A-11R-2232-07",
                                  "TCGA-CN-4722-01A-01R-1435-13",
                                  "TCGA-CN-4722-01A-01R-1436-07",
                                  "TCGA-CN-4723-01A-01R-1435-13",
                                  "TCGA-CN-4723-01A-01R-1436-07",
                                  "TCGA-CN-4727-01A-01R-1435-13",
                                  "TCGA-CN-4727-01A-01R-1436-07",
                                  "TCGA-CN-4735-01A-01R-1435-13",
                                  "TCGA-CN-4735-01A-01R-1436-07",
                                  "TCGA-CN-4738-01A-02R-1513-13",
                                  "TCGA-CN-4738-01A-02R-1514-07",
                                  "TCGA-CN-4739-01A-02R-1513-13",
                                  "TCGA-CN-4739-01A-02R-1514-07",
                                  "TCGA-CN-5355-01A-01R-1435-13",
                                  "TCGA-CN-5355-01A-01R-1436-07",
                                  "TCGA-CN-5356-01A-01R-1435-13",
                                  "TCGA-CN-5356-01A-01R-1436-07",
                                  "TCGA-CN-5360-01A-01R-1435-13",
                                  "TCGA-CN-5360-01A-01R-1436-07",
                                  "TCGA-CN-5361-01A-01R-1435-13",
                                  "TCGA-CN-5361-01A-01R-1436-07",
                                  "TCGA-CN-5363-01A-01R-1435-13",
                                  "TCGA-CN-5363-01A-01R-1436-07",
                                  "TCGA-CN-6010-01A-11R-1685-13",
                                  "TCGA-CN-6010-01A-11R-1686-07",
                                  "TCGA-CN-6012-01A-11R-1685-13",
                                  "TCGA-CN-6012-01A-11R-1686-07",
                                  "TCGA-CN-6021-01A-11R-1685-13",
                                  "TCGA-CN-6021-01A-11R-1686-07",
                                  "TCGA-CN-6022-01A-21R-1685-13",
                                  "TCGA-CN-6022-01A-21R-1686-07",
                                  "TCGA-CN-6023-01A-11R-1685-13",
                                  "TCGA-CN-6023-01A-11R-1686-07",
                                  "TCGA-CN-6988-01A-11R-1914-13",
                                  "TCGA-CN-6988-01A-11R-1915-07",
                                  "TCGA-CN-6989-01A-11R-1914-13",
                                  "TCGA-CN-6989-01A-11R-1915-07",
                                  "TCGA-CN-6992-01A-11R-1914-13",
                                  "TCGA-CN-6992-01A-11R-1915-07",
                                  "TCGA-CN-6997-01A-11R-2015-13",
                                  "TCGA-CN-6997-01A-11R-2016-07",
                                  "TCGA-CN-A497-01A-11R-A24H-07",
                                  "TCGA-CN-A497-01A-11R-A24J-13",
                                  "TCGA-CN-A49B-01A-31R-A24H-07",
                                  "TCGA-CN-A49B-01A-31R-A24J-13",
                                  "TCGA-CN-A63T-01A-11R-A28V-07",
                                  "TCGA-CN-A63T-01A-11R-A28Z-13",
                                  "TCGA-CN-A63U-01A-11R-A30B-07",
                                  "TCGA-CN-A63U-01A-11R-A30I-13",
                                  "TCGA-CN-A63W-01A-11R-A30B-07",
                                  "TCGA-CN-A63W-01A-11R-A30I-13",
                                  "TCGA-CN-A641-01A-11R-A30B-07",
                                  "TCGA-CN-A641-01A-11R-A30I-13",
                                  "TCGA-CN-A6V3-01A-12R-A34O-13",
                                  "TCGA-CN-A6V3-01A-12R-A34R-07",
                                  "TCGA-CR-6474-01A-11R-1872-13",
                                  "TCGA-CR-6474-01A-11R-1873-07",
                                  "TCGA-CR-7364-01A-11R-2015-13",
                                  "TCGA-CR-7364-01A-11R-2016-07",
                                  "TCGA-CR-7370-01A-11R-2131-13",
                                  "TCGA-CR-7370-01A-11R-2132-07",
                                  "TCGA-CR-7371-01A-11R-2015-13",
                                  "TCGA-CR-7371-01A-11R-2016-07",
                                  "TCGA-CR-7374-01A-11R-2015-13",
                                  "TCGA-CR-7374-01A-11R-2016-07",
                                  "TCGA-CR-7388-01A-11R-2015-13",
                                  "TCGA-CR-7388-01A-11R-2016-07",
                                  "TCGA-CR-7389-01A-11R-2015-13",
                                  "TCGA-CR-7389-01A-11R-2016-07",
                                  "TCGA-CR-7398-01A-11R-2015-13",
                                  "TCGA-CR-7398-01A-11R-2016-07",
                                  "TCGA-CR-7399-01A-11R-2015-13",
                                  "TCGA-CR-7399-01A-11R-2016-07",
                                  "TCGA-CR-7402-01A-11R-2015-13",
                                  "TCGA-CR-7402-01A-11R-2016-07",
                                  "TCGA-CR-7404-01A-11R-2131-13",
                                  "TCGA-CR-7404-01A-11R-2132-07",
                                  "TCGA-CV-5430-01A-02R-1685-13",
                                  "TCGA-CV-5430-01A-02R-1686-07",
                                  "TCGA-CV-5431-01A-01R-1513-13",
                                  "TCGA-CV-5431-01A-01R-1514-07",
                                  "TCGA-CV-5432-01A-02R-1685-13",
                                  "TCGA-CV-5432-01A-02R-1686-07",
                                  "TCGA-CV-5434-01A-01R-1685-13",
                                  "TCGA-CV-5434-01A-01R-1686-07",
                                  "TCGA-CV-5435-01A-01R-1685-13",
                                  "TCGA-CV-5435-01A-01R-1686-07",
                                  "TCGA-CV-5440-01A-01R-1513-13",
                                  "TCGA-CV-5440-01A-01R-1514-07",
                                  "TCGA-CV-5441-01A-01R-1513-13",
                                  "TCGA-CV-5441-01A-01R-1514-07",
                                  "TCGA-CV-5443-01A-01R-1513-13",
                                  "TCGA-CV-5443-01A-01R-1514-07",
                                  "TCGA-CV-5444-01A-02R-1513-13",
                                  "TCGA-CV-5444-01A-02R-1514-07",
                                  "TCGA-CV-5978-01A-11R-1685-13",
                                  "TCGA-CV-5978-01A-11R-1686-07",
                                  "TCGA-CV-6935-01A-11R-1914-13",
                                  "TCGA-CV-6935-01A-11R-1915-07",
                                  "TCGA-CV-6935-11A-01R-1914-13",
                                  "TCGA-CV-6935-11A-01R-1915-07",
                                  "TCGA-CV-6962-01A-11R-1914-13",
                                  "TCGA-CV-6962-01A-11R-1915-07",
                                  "TCGA-CV-6962-11A-01R-1914-13",
                                  "TCGA-CV-6962-11A-01R-1915-07",
                                  "TCGA-CV-7089-01A-11R-2015-13",
                                  "TCGA-CV-7089-01A-11R-2016-07",
                                  "TCGA-CV-7101-01A-11R-2015-13",
                                  "TCGA-CV-7101-01A-11R-2016-07",
                                  "TCGA-CV-7101-11A-01R-2015-13",
                                  "TCGA-CV-7101-11A-01R-2016-07",
                                  "TCGA-CV-7177-01A-11R-2015-13",
                                  "TCGA-CV-7177-01A-11R-2016-07",
                                  "TCGA-CV-7177-11A-01R-2015-13",
                                  "TCGA-CV-7177-11A-01R-2016-07",
                                  "TCGA-CV-7242-01A-11R-2015-13",
                                  "TCGA-CV-7242-01A-11R-2016-07",
                                  "TCGA-CV-7242-11A-01R-2015-13",
                                  "TCGA-CV-7242-11A-01R-2016-07",
                                  "TCGA-CV-7245-01A-11R-2015-13",
                                  "TCGA-CV-7245-01A-11R-2016-07",
                                  "TCGA-CV-7245-11A-01R-2015-13",
                                  "TCGA-CV-7245-11A-01R-2016-07",
                                  "TCGA-CV-7247-01A-11R-2015-13",
                                  "TCGA-CV-7247-01A-11R-2016-07",
                                  "TCGA-CV-7248-01A-11R-2015-13",
                                  "TCGA-CV-7248-01A-11R-2016-07",
                                  "TCGA-CV-7250-01A-11R-2015-13",
                                  "TCGA-CV-7250-01A-11R-2016-07",
                                  "TCGA-CV-7250-11A-01R-2015-13",
                                  "TCGA-CV-7250-11A-01R-2016-07",
                                  "TCGA-CV-7261-01A-11R-2015-13",
                                  "TCGA-CV-7261-01A-11R-2016-07",
                                  "TCGA-CV-7261-11A-01R-2015-13",
                                  "TCGA-CV-7261-11A-01R-2016-07",
                                  "TCGA-CV-7410-01A-21R-2080-13",
                                  "TCGA-CV-7410-01A-21R-2081-07",
                                  "TCGA-CV-7415-01A-11R-2080-13",
                                  "TCGA-CV-7415-01A-11R-2081-07",
                                  "TCGA-CV-7418-01A-11R-2080-13",
                                  "TCGA-CV-7418-01A-11R-2081-07",
                                  "TCGA-CV-7421-01A-11R-2080-13",
                                  "TCGA-CV-7421-01A-11R-2081-07",
                                  "TCGA-CV-7422-01A-21R-2080-13",
                                  "TCGA-CV-7422-01A-21R-2081-07",
                                  "TCGA-CV-7424-01A-11R-2080-13",
                                  "TCGA-CV-7424-01A-11R-2081-07",
                                  "TCGA-CV-7424-11A-01R-2080-13",
                                  "TCGA-CV-7424-11A-01R-2081-07",
                                  "TCGA-CV-7430-01A-11R-2131-13",
                                  "TCGA-CV-7430-01A-11R-2132-07",
                                  "TCGA-CV-7433-01A-11R-2131-13",
                                  "TCGA-CV-7433-01A-11R-2132-07",
                                  "TCGA-CV-7437-01A-21R-2131-13",
                                  "TCGA-CV-7437-01A-21R-2132-07",
                                  "TCGA-CV-7437-11A-01R-2131-13",
                                  "TCGA-CV-7437-11A-01R-2132-07",
                                  "TCGA-CV-7440-01A-11R-2131-13",
                                  "TCGA-CV-7440-01A-11R-2132-07",
                                  "TCGA-CV-7440-11A-01R-2186-13",
                                  "TCGA-CV-7440-11A-01R-2187-07",
                                  "TCGA-CV-A45W-01A-11R-A24Z-07",
                                  "TCGA-CV-A45W-01A-11R-A25B-13",
                                  "TCGA-CV-A45Z-01A-21R-A24Z-07",
                                  "TCGA-CV-A45Z-01A-21R-A25B-13",
                                  "TCGA-CV-A460-01A-21R-A24Z-07",
                                  "TCGA-CV-A460-01A-21R-A25B-13",
                                  "TCGA-CV-A461-01A-41R-A25Z-13",
                                  "TCGA-CV-A461-01A-41R-A266-07",
                                  "TCGA-CV-A6K1-01A-11R-A31N-07",
                                  "TCGA-CV-A6K1-01A-11R-A31R-13",
                                  "TCGA-D6-6517-01A-11R-1872-13",
                                  "TCGA-D6-6517-01A-11R-1873-07",
                                  "TCGA-D6-6824-01A-11R-1914-13",
                                  "TCGA-D6-6824-01A-11R-1915-07",
                                  "TCGA-D6-6826-01A-11R-1914-13",
                                  "TCGA-D6-6826-01A-11R-1915-07",
                                  "TCGA-D6-8568-01A-11R-2401-13",
                                  "TCGA-D6-8568-01A-11R-2403-07",
                                  "TCGA-D6-A6EK-01A-11R-A31N-07",
                                  "TCGA-D6-A6EK-01A-11R-A31R-13",
                                  "TCGA-D6-A6EQ-01A-11R-A31N-07",
                                  "TCGA-D6-A6EQ-01A-11R-A31R-13",
                                  "TCGA-D6-A6ES-01A-12R-A31N-07",
                                  "TCGA-D6-A6ES-01A-12R-A31R-13",
                                  "TCGA-D6-A74Q-01A-11R-A34O-13",
                                  "TCGA-D6-A74Q-01A-11R-A34R-07",
                                  "TCGA-DQ-5629-01A-01R-1872-13",
                                  "TCGA-DQ-5629-01A-01R-1873-07",
                                  "TCGA-DQ-7589-01A-11R-2231-13",
                                  "TCGA-DQ-7589-01A-11R-2232-07",
                                  "TCGA-DQ-7595-01A-11R-2231-13",
                                  "TCGA-DQ-7595-01A-11R-2232-07",
                                  "TCGA-F7-7848-01A-11R-2131-13",
                                  "TCGA-F7-7848-01A-11R-2132-07",
                                  "TCGA-F7-8298-01A-11R-2401-13",
                                  "TCGA-F7-8298-01A-11R-2403-07",
                                  "TCGA-F7-A50I-01A-11R-A28V-07",
                                  "TCGA-F7-A50I-01A-11R-A28Z-13",
                                  "TCGA-F7-A622-01A-11R-A28V-07",
                                  "TCGA-F7-A622-01A-11R-A28Z-13",
                                  "TCGA-F7-A623-01A-11R-A28V-07",
                                  "TCGA-F7-A623-01A-11R-A28Z-13",
                                  "TCGA-H7-A6C5-01A-11R-A30I-13",
                                  "TCGA-H7-A6C5-11A-11R-A30B-07",
                                  "TCGA-H7-A6C5-11A-11R-A30I-13",
                                  "TCGA-HD-7229-01A-11R-2015-13",
                                  "TCGA-HD-7229-01A-11R-2016-07",
                                  "TCGA-KU-A66S-01A-21R-A30B-07",
                                  "TCGA-KU-A66S-01A-21R-A30I-13",
                                  "TCGA-QK-A8Z8-01A-11R-A39B-13",
                                  "TCGA-QK-A8Z8-01A-11R-A39I-07",
                                  "TCGA-QK-AA3J-01A-11R-A39B-13",
                                  "TCGA-QK-AA3J-01A-11R-A39I-07",
                                  "TCGA-T3-A92M-01A-31R-A39B-13",
                                  "TCGA-T3-A92M-01A-31R-A39I-07",
                                  "TCGA-TN-A7HJ-01A-12R-A34O-13",
                                  "TCGA-TN-A7HJ-01A-12R-A34R-07",
                                  "TCGA-UF-A718-01A-22R-A34O-13",
                                  "TCGA-UF-A718-01A-22R-A34R-07",
                                  "TCGA-UF-A71D-01A-12R-A34O-13",
                                  "TCGA-UF-A71D-01A-12R-A34R-07",
                                  "TCGA-UF-A7J9-01A-12R-A34O-13",
                                  "TCGA-UF-A7J9-01A-12R-A34R-07",
                                  "TCGA-UF-A7JF-01A-11R-A34O-13",
                                  "TCGA-UF-A7JF-01A-11R-A34R-07",
                                  "TCGA-UF-A7JH-01A-21R-A34O-13",
                                  "TCGA-UF-A7JH-01A-21R-A34R-07",
                                  "TCGA-UF-A7JJ-01A-11R-A34O-13",
                                  "TCGA-UF-A7JJ-01A-11R-A34R-07",
                                  "TCGA-UF-A7JK-01A-11R-A34O-13",
                                  "TCGA-UF-A7JK-01A-11R-A34R-07"))



getResults(query_rna)

#download files -- remember, this step should be done just once =)
#GDCdownload(query_rna)
##GDCdownload will download 127 files.

getwd()
#prepare data from these files (preciso estar no wd, senão dá erro!)

tcga_hnsc_rnadata <- GDCprepare(query_rna, 
                                summarizedExperiment = TRUE, 
                                save=T, 
                                save.filename="rna_larynx.rda")
#Available assays in SummarizedExperiment : 
#=> unstranded
#=> stranded_first
#=> stranded_second
#=> tpm_unstrand
#=> fpkm_unstrand
#=> fpkm_uq_unstrand

## If I already want to save a rna matrix, I can just do this:
hnsc_matrix <- assay(tcga_hnsc_rnadata, 'fpkm_unstrand')



###########################################################################################


# --------   BUILDING A QUERY -- METHYLATION DATA

query_TCGA_methyl <- GDCquery(project = 'TCGA-HNSC',
                              data.category = 'DNA Methylation',
                              platform = 'Illumina Human Methylation 450',
                              access = 'open',
                              data.type = 'Methylation Beta Value',
                              barcode = c("TCGA-BA-4076-01A-01D-1433-05",
                                          "TCGA-BA-4078-01A-01D-1433-05",
                                          "TCGA-BA-5555-01A-01D-1511-05",
                                          "TCGA-BA-6868-01B-12D-1913-05",
                                          "TCGA-BA-6869-01A-11D-1871-05",
                                          "TCGA-BA-6870-01A-11D-1871-05",
                                          "TCGA-BA-A6DA-01A-31D-A31M-05",
                                          "TCGA-BA-A6DI-01A-11D-A30F-05",
                                          "TCGA-BB-4217-01A-11D-2079-05",
                                          "TCGA-BB-7862-01A-21D-2230-05",
                                          "TCGA-BB-7864-01A-11D-2230-05",
                                          "TCGA-BB-7870-01A-11D-2230-05",
                                          "TCGA-CN-4722-01A-01D-1433-05",
                                          "TCGA-CN-4723-01A-01D-1433-05",
                                          "TCGA-CN-4727-01A-01D-1433-05",
                                          "TCGA-CN-4735-01A-01D-1433-05",
                                          "TCGA-CN-4738-01A-02D-1511-05",
                                          "TCGA-CN-4739-01A-02D-1511-05",
                                          "TCGA-CN-5355-01A-01D-1433-05",
                                          "TCGA-CN-5356-01A-01D-1433-05",
                                          "TCGA-CN-5360-01A-01D-1433-05",
                                          "TCGA-CN-5361-01A-01D-1433-05",
                                          "TCGA-CN-5363-01A-01D-1433-05",
                                          "TCGA-CN-6010-01A-11D-1684-05",
                                          "TCGA-CN-6012-01A-11D-1684-05",
                                          "TCGA-CN-6021-01A-11D-1684-05",
                                          "TCGA-CN-6022-01A-21D-1684-05",
                                          "TCGA-CN-6023-01A-11D-1684-05",
                                          "TCGA-CN-6988-01A-11D-1913-05",
                                          "TCGA-CN-6989-01A-11D-1913-05",
                                          "TCGA-CN-6992-01A-11D-1913-05",
                                          "TCGA-CN-6997-01A-11D-2014-05",
                                          "TCGA-CN-A497-01A-11D-A24I-05",
                                          "TCGA-CN-A49B-01A-31D-A24I-05",
                                          "TCGA-CN-A63T-01A-11D-A28S-05",
                                          "TCGA-CN-A63U-01A-11D-A30F-05",
                                          "TCGA-CN-A63W-01A-11D-A30F-05",
                                          "TCGA-CN-A641-01A-11D-A30F-05",
                                          "TCGA-CN-A6V3-01A-12D-A34K-05",
                                          "TCGA-CR-6474-01A-11D-1871-05",
                                          "TCGA-CR-7364-01A-11D-2014-05",
                                          "TCGA-CR-7370-01A-11D-2130-05",
                                          "TCGA-CR-7371-01A-11D-2014-05",
                                          "TCGA-CR-7374-01A-11D-2014-05",
                                          "TCGA-CR-7388-01A-11D-2014-05",
                                          "TCGA-CR-7389-01A-11D-2014-05",
                                          "TCGA-CR-7398-01A-11D-2014-05",
                                          "TCGA-CR-7399-01A-11D-2014-05",
                                          "TCGA-CR-7402-01A-11D-2014-05",
                                          "TCGA-CR-7404-01A-11D-2130-05",
                                          "TCGA-CV-5430-01A-02D-1684-05",
                                          "TCGA-CV-5430-11B-01D-1684-05",
                                          "TCGA-CV-5431-01A-01D-1511-05",
                                          "TCGA-CV-5431-11A-01D-1511-05",
                                          "TCGA-CV-5432-01A-02D-1684-05",
                                          "TCGA-CV-5432-11B-01D-1684-05",
                                          "TCGA-CV-5434-01A-01D-1684-05",
                                          "TCGA-CV-5434-11B-01D-1684-05",
                                          "TCGA-CV-5435-01A-01D-1684-05",
                                          "TCGA-CV-5435-11B-01D-1684-05",
                                          "TCGA-CV-5440-01A-01D-1511-05",
                                          "TCGA-CV-5440-11A-01D-1511-05",
                                          "TCGA-CV-5441-01A-01D-1511-05",
                                          "TCGA-CV-5441-11A-01D-1511-05",
                                          "TCGA-CV-5443-01A-01D-1511-05",
                                          "TCGA-CV-5443-11A-01D-1511-05",
                                          "TCGA-CV-5444-01A-02D-1511-05",
                                          "TCGA-CV-5444-11A-01D-1511-05",
                                          "TCGA-CV-5978-01A-11D-1684-05",
                                          "TCGA-CV-5978-11A-01D-1684-05",
                                          "TCGA-CV-6935-01A-11D-1913-05",
                                          "TCGA-CV-6935-11A-01D-1913-05",
                                          "TCGA-CV-6962-01A-11D-1913-05",
                                          "TCGA-CV-6962-11A-01D-1913-05",
                                          "TCGA-CV-7089-01A-11D-2014-05",
                                          "TCGA-CV-7089-11A-01D-2014-05",
                                          "TCGA-CV-7101-01A-11D-2014-05",
                                          "TCGA-CV-7101-11A-01D-2014-05",
                                          "TCGA-CV-7177-01A-11D-2014-05",
                                          "TCGA-CV-7242-01A-11D-2014-05",
                                          "TCGA-CV-7245-01A-11D-2014-05",
                                          "TCGA-CV-7245-11A-01D-2014-05",
                                          "TCGA-CV-7247-01A-11D-2014-05",
                                          "TCGA-CV-7248-01A-11D-2014-05",
                                          "TCGA-CV-7250-01A-11D-2014-05",
                                          "TCGA-CV-7250-11A-01D-2014-05",
                                          "TCGA-CV-7261-01A-11D-2014-05",
                                          "TCGA-CV-7410-01A-21D-2079-05",
                                          "TCGA-CV-7415-01A-11D-2079-05",
                                          "TCGA-CV-7418-01A-11D-2079-05",
                                          "TCGA-CV-7421-01A-11D-2079-05",
                                          "TCGA-CV-7422-01A-21D-2079-05",
                                          "TCGA-CV-7424-01A-11D-2079-05",
                                          "TCGA-CV-7430-01A-11D-2130-05",
                                          "TCGA-CV-7433-01A-11D-2130-05",
                                          "TCGA-CV-7437-01A-21D-2130-05",
                                          "TCGA-CV-7440-01A-11D-2130-05",
                                          "TCGA-CV-A45W-01A-11D-A254-05",
                                          "TCGA-CV-A45Z-01A-21D-A254-05",
                                          "TCGA-CV-A460-01A-21D-A254-05",
                                          "TCGA-CV-A461-01A-41D-A265-05",
                                          "TCGA-CV-A6K1-01A-11D-A31M-05",
                                          "TCGA-D6-6517-01A-11D-1871-05",
                                          "TCGA-D6-6824-01A-11D-1913-05",
                                          "TCGA-D6-6826-01A-11D-1913-05",
                                          "TCGA-D6-8568-01A-11D-2398-05",
                                          "TCGA-D6-A6EK-01A-11D-A31M-05",
                                          "TCGA-D6-A6EQ-01A-11D-A31M-05",
                                          "TCGA-D6-A6ES-01A-12D-A31M-05",
                                          "TCGA-D6-A74Q-01A-11D-A34K-05",
                                          "TCGA-DQ-5629-01A-01D-1871-05",
                                          "TCGA-DQ-7589-01A-11D-2230-05",
                                          "TCGA-DQ-7595-01A-11D-2230-05",
                                          "TCGA-F7-7848-01A-11D-2130-05",
                                          "TCGA-F7-8298-01A-11D-2398-05",
                                          "TCGA-F7-A50I-01A-11D-A28S-05",
                                          "TCGA-F7-A622-01A-11D-A28S-05",
                                          "TCGA-F7-A623-01A-11D-A28S-05",
                                          "TCGA-H7-A6C5-01A-11D-A30F-05",
                                          "TCGA-HD-7229-01A-11D-2014-05",
                                          "TCGA-KU-A66S-01A-21D-A30F-05",
                                          "TCGA-QK-A8Z8-01A-11D-A392-05",
                                          "TCGA-QK-A8ZB-01A-11D-A392-05",
                                          "TCGA-QK-AA3J-01A-11D-A392-05",
                                          "TCGA-T3-A92M-01A-31D-A392-05",
                                          "TCGA-TN-A7HJ-01A-12D-A34K-05",
                                          "TCGA-UF-A718-01A-22D-A34K-05",
                                          "TCGA-UF-A71D-01A-12D-A34K-05",
                                          "TCGA-UF-A7J9-01A-12D-A34K-05",
                                          "TCGA-UF-A7JF-01A-11D-A34K-05",
                                          "TCGA-UF-A7JH-01A-21D-A34K-05",
                                          "TCGA-UF-A7JJ-01A-11D-A34K-05",
                                          "TCGA-UF-A7JK-01A-11D-A34K-05"))

output_query_methyl <- getResults(query_TCGA_methyl)

#download methylation data
#GDCdownload(query_TCGA_methyl)

tcga_hnsc_methyldata <- GDCprepare(query_TCGA_methyl, 
                                   summarizedExperiment = T, 
                                   save=T, 
                                   save.filename="methyl_larynx.rda")
# MAtrix of BETA VALUES 
beta <- assay(tcga_hnsc_methyldata)

idx <- tcga_hnsc_methyldata %>%
  assay %>%
  rowVars() %>%
  order(decreasing = T) %>%
  head(10)

##plot
install.packages("pheatmap")
library(pheatmap)

pheatmap(assay(tcga_hnsc_methyldata)[idx,])



# --------   BUILDING A QUERY -- GENOMIC DATA
getProjectSummary('TCGA-HNSC')

query <- GDCquery(
  project = "TCGA-HNSC", 
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")

GDCdownload(query)

#Will be used for mutational signatures code!
maf <- GDCprepare(query)
