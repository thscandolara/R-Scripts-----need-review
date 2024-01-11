## Prediction Analysis of Microarrays for R (pamr) 

setwd()

##Call packages
library(ChAMP)
library(minfi)
library(limma)
library(ggplot2)
library(RColorBrewer)
library(edgeR)
library(survival)
library(survminer)
library(ggplot2)
library(pamr)



###Load data from TCGAquery e TCGAprepare:
load("~/DNAmethylation_larynx.rda")
dna_methyl_full <- data
rm(data)



#____Criando a matrix com os valores de metilação (beta-values)
beta <- as.data.frame(SummarizedExperiment::assay(dna_methyl_full))


#____Criando um df com os dados clinicos
clinical <- data.frame(dna_methyl_full@colData)


##Excluir amostras de pacientes classificados como basalóide e/ou tinham outro sítio que não laringe (tonsil)

clinical <- clinical[!(rownames(clinical) %in% c("TCGA-CR-7404-01A-11D-2130-05", "TCGA-KU-A66S-01A-21D-A30F-05", "TCGA-QK-AA3J-01A-11D-A392-05")), ]



# Filtrar o dataset "beta" com base nas rownames do dataset "clinical" filtrado:
beta <- beta[, rownames(clinical)]



#__Selecionando todas as sondas as quais tiveram valores nulos/zeros/NAs na lista 'probe.na'
probe.na <- rowSums(is.na(beta))

#__Manter na minha tabela 'beta' apenas as sondas que não tinham valores NA na lista 'probe.na'
probe <- probe.na[probe.na == 0]
beta <- beta[row.names(beta) %in% names(probe), ]

#__Transformando minha df em uma matriz (input necessário no champ.filter)
beta_matrix <- as.matrix(beta)



##Pacote ChAMP: função champ.filter
##_________Filtrar as CpGs -- filtrar para NoCG, SNPs start, MultiHit start, XY start..
met <- champ.filter(beta = beta_matrix,
                    pd = clinical,
                    filterXY = T,
                    filterNoCG = T,
                    filterSNPs = T,
                    filterMultiHit = T,
                    population = NULL,
                    filterDetP = TRUE)



###Transformar estas colunas em fatores
as.factor(met$pd$ajcc_pathologic_stage)
as.factor(met$pd$ajcc_pathologic_n)
as.factor(met$pd$shortLetterCode)


## Dicotomizar a coluna de classificação de invasão linfonodal
no_invasion <- c("N0")
nodal_invasion <- c("N1", "N2", "N2a", "N2b", "N2c", "N3")

met$pd$ajcc_pathologic_n <- ifelse(met$pd$ajcc_pathologic_n %in% no_invasion, "No Invasion",
                                   ifelse(met$pd$ajcc_pathologic_n %in% nodal_invasion, "Nodal Invasion", NA))

as.factor(met$pd$ajcc_pathologic_n)



##_____________REMOVER AMOSTRAS NORMAIS DA ANÁLISE (CONFUSÃO)____________________##

# Obter índices das linhas a serem removidas (amostras normais)
indices_normais <- which(met$pd$shortLetterCode == "NT")

# Remover as linhas correspondentes no met$beta
met$beta <- met$beta[, -indices_normais]

# Remover as linhas correspondentes no met$pd$ajcc_pathologic_n
met$pd <- subset(met$pd, !met$pd$shortLetterCode %in% "NT")

#Se eu quiser refazer a análise contendo amostras normais tbm, preciso criar um novo arquivo met, sem esse filtro






###_______________ELIMINAR VALORES 'NA' DA COLUNA 'ajcc_pathologic_n' E MANTER MESMAS AMOSTRAS NA MATRIZ BETA______


# Filtrar linhas com valores não NA em ajcc_pathologic_n
non_na_rows <- complete.cases(met$pd$ajcc_pathologic_n)
met$pd <- met$pd[non_na_rows, ]

# Obter os nomes das linhas não NA em ajcc_pathologic_n
rows_to_keep <- rownames(met$pd)

# Filtrar as colunas de met$beta usando as linhas selecionadas
met$beta <- met$beta[, intersect(colnames(met$beta), rows_to_keep)]


# Contar o número de ocorrências de "Nodal Invasion" e "No Invasion" em ajcc_pathologic_n
sum(met$pd$ajcc_pathologic_n == "Nodal Invasion", na.rm = TRUE)#54
sum(met$pd$ajcc_pathologic_n == "No Invasion", na.rm = TRUE)#39




###_____________________PAM - ANÁLISE COM DADOS DE INVASÃO LINFONODAL _____________________________


pamr_data <- list(x=met$beta,y=factor(met$pd$ajcc_pathologic_n), 
                  geneid=rownames(met$beta),
                  genenames=paste("cg",as.character(1:nrow(met$beta)),sep=""))


# Realizar a call com a função pamr.train()
train <- pamr.train(data = pamr_data)
train


## Cross-validate the classifier
results <- pamr.cv(train,pamr_data) 
results

## Plot the cross-validated error curves
pamr.plotcv(results)

## Compute the confusion matrix for a particular model (threshold=2.1) 
pamr.confusion(results, threshold=2.1)

#                     No Invasion     Nodal Invasion      Class Error rate
#No Invasion             21                  18              0.4615385
#Nodal Invasion          17                  37              0.3148148
#Overall error rate= 0.373



## Plot the cross-validated class probabilities by class(threshold=2.1) 
pamr.plotcvprob(results, pamr_data, threshold=2.1)


## Plot the class centroids
pamr.plotcen(train, pamr_data, threshold=2.1)

## Make a gene plot of the most significant genes
pamr.geneplot(train, pamr_data, threshold= 2.1)
##Muitos genes, não vai dar pra plotar


# Estimate false discovery rates and plot them
fdr.obj<- pamr.fdr(train, pamr_data)
str(fdr.obj)
pamr.plotfdr(fdr.obj)#NULL

##Gerar um gráfico
plot(fdr.obj$results[, "Threshold"], fdr.obj$results[, "Median FDR"], type = "l", 
     xlab = "Threshold", ylab = "Median FDR", main = "False Discovery Rates")


## List the significant genes
pamrlist_nodal <- pamr.listgenes(train, pamr_data, threshold=2.1)

write.csv(pamrlist_nodal, "pamr.list_sondas_genes.csv")






###_____________________PAM - ANÁLISE COM DADOS DE ESTADIAMENTO _____________________________

#Chamar novamente toda a matriz beta e dados clinicos para não dar confusão...
#Rodar da filtragem de amostras normais

met <- champ.filter(beta = beta_matrix,
                    pd = clinical,
                    filterXY = T,
                    filterNoCG = T,
                    filterSNPs = T,
                    filterMultiHit = T,
                    population = NULL,
                    filterDetP = TRUE)



###Transformar estas colunas em fatores
as.factor(met$pd$ajcc_pathologic_stage)
as.factor(met$pd$ajcc_pathologic_n)
as.factor(met$pd$shortLetterCode)


##Dicotomizar coluna de estadiamento patológico
early_stage <- c("Stage I", "Stage II")
advanced_stage <- c("Stage III", "Stage IVA", "Stage IVB", "Stage IVC")

met$pd$ajcc_pathologic_stage <- ifelse(met$pd$ajcc_pathologic_stage %in% early_stage, "Early Stage",
                                       ifelse(met$pd$ajcc_pathologic_stage %in% advanced_stage, "Advanced Stage", NA))

as.factor(met$pd$ajcc_pathologic_stage)
str(met$pd$ajcc_pathologic_stage) #113


###Remover normais

# Obter índices das linhas a serem removidas (amostras normais)
indices_normais <- which(met$pd$shortLetterCode == "NT")

# Remover as linhas correspondentes no met$beta
met$beta <- met$beta[, -indices_normais]

# Remover as linhas correspondentes no met$pd$ajcc_pathologic_n
met$pd <- subset(met$pd, !met$pd$shortLetterCode %in% "NT")

#Se eu quiser refazer a análise contendo amostras normais tbm, preciso criar um novo arquivo met, sem esse filtro


###_______________ELIMINAR VALORES 'NA' E MANTER MESMAS AMOSTRAS NA MATRIZ BETA______

# Filtrar linhas com valores não NA em ajcc_pathologic_stage
non_na_rows2 <- complete.cases(met$pd$ajcc_pathologic_stage)
met$pd <- met$pd[non_na_rows2, ]

# Obter os nomes das linhas não NA em ajcc_pathologic_stage
rows_to_keep2 <- rownames(met$pd)

# Filtrar as colunas de met$beta usando as linhas selecionadas
met$beta <- met$beta[, intersect(colnames(met$beta), rows_to_keep2)]


# Contar o número de ocorrências de "Early Stage" e "Advanced Stage" em ajcc_pathologic_stage
sum(met$pd$ajcc_pathologic_stage == "Early Stage", na.rm = TRUE)#11
sum(met$pd$ajcc_pathologic_stage == "Advanced Stage", na.rm = TRUE)#86



str(met$beta) 
str(met$pd$ajcc_pathologic_stage)


##Criar o pam_data
pamr_data2 <- list(x=met$beta,y=met$pd$ajcc_pathologic_stage, geneid=as.character(1:nrow(met$beta)),
                   genenames=paste("g",as.character(1:nrow(met$beta)),sep=""))


# Realizar a call com a função pamr.train()
train2 <- pamr.train(data = pamr_data2)
train2


## Cross-validate the classifier
results2 <- pamr.cv(train2,pamr_data2) 
results2

## Plot the cross-validated error curves
pamr.plotcv(results2)

## Compute the confusion matrix for a particular model (threshold=2.5) 
pamr.confusion(results2, threshold=1.55)
#                 Advanced Stage      Early Stage  Class Error rate
#Advanced Stage             57          29            0.3372093
#Early Stage                 7           4            0.6363636
#Overall error rate= 0.367



## Plot the cross-validated class probabilities by class
pamr.plotcvprob(results2, pamr_data2, threshold=1.55)


## Plot the class centroids
pamr.plotcen(train2, pamr_data2, threshold=1.55)


## Make a gene plot of the most significant genes
pamr.geneplot(train2, pamr_data2, threshold= 1.55)
##Muita coisa pra plotar, n dá certo

# Estimate false discovery rates and plot them
fdr.obj2 <- pamr.fdr(train2, pamr_data2)

str(fdr.obj2)
pamr.plotfdr(fdr.obj2)#NULL

##Gerar um gráfico
plot(fdr.obj2$results[, "Threshold"], fdr.obj2$results[, "Median FDR"], type = "l", 
     xlab = "Threshold", ylab = "Median FDR", main = "False Discovery Rates")



## List the significant genes
pamlist_stage <- pamr.listgenes(train2, pamr_data2, threshold=1.55)
#Error in dimnames(res) <- list(NULL, c("id", gnhdr, schdr)) : 
#length of 'dimnames' [2] not equal to array extent

write.csv(pamlist_stage, "pamr.list_sondas_stage.csv")





###################################################################################################################################################


##___________________________________________________________Análise PAM apenas com sondas que foram correlacionadas com expressão genica
##Preciso ter o dataframe "clinical.met"

clin.surv <- data.frame(
  vital_status = clinical.met$vital_status,
  days_to_death = clinical.met$days_to_death,
  year_of_birth = clinical.met$year_of_birth,
  year_of_death = clinical.met$year_of_death,
  cigarettes_per_day = clinical.met$cigarettes_per_day,
  alcohol_history = clinical.met$alcohol_history,
  ajcc_clinical_stage = clinical.met$ajcc_clinical_stage,
  year_of_diagnosis = clinical.met$year_of_diagnosis,
  ajcc_clinical_n = clinical.met$ajcc_clinical_n,
  ajcc_pathologic_n = clinical.met$ajcc_pathologic_n,
  years_smoked = clinical.met$years_smoked,
  barcode = clinical.met$barcode,
  race = clinical.met$race,
  gender = clinical.met$gender,
  ethnicity = clinical.met$ethnicity
)

rownames(clin.surv) <- clin.surv$barcode
rownames(clin.surv) <- substr(rownames(clin.surv),1,19)
clin.surv.rows <- rownames(clin.surv)

#Manter apenas os valores de beta de amostras que estão contidas no clin.surv.rows.
beta_m_surv <- beta_m_surv[, colnames(beta_m_surv) %in% clin.surv.rows]

#Manter as amostras em clin.surv iguais as amostras em beta_m_surv.
clin.surv <- clin.surv[rownames(clin.surv) %in% colnames(beta_m_surv), ]


# Checando se o resultado deu certo 
all(colnames(beta_m_surv) %in% rownames(clin.surv))#TRUE
all(rownames(clin.surv) %in% colnames(beta_m_surv))#TRUE

###Transformar estas colunas em fatores
as.factor(clin.surv$ajcc_pathologic_n)


## Dicotomizar a coluna de classificação de invasão linfonodal
no_invasion1 <- c("N0")
nodal_invasion1 <- c("N1", "N2", "N2a", "N2b", "N2c", "N3")

clin.surv$ajcc_pathologic_n <- ifelse(clin.surv$ajcc_pathologic_n %in% no_invasion1, "No Invasion",
                                      ifelse(clin.surv$ajcc_pathologic_n %in% nodal_invasion1, "Nodal Invasion", NA))

as.factor(clin.surv$ajcc_pathologic_n)
clin.surv<- as.data.frame(clin.surv)


##Dicotomizar coluna de estadiamento patológico
early_stage1 <- c("Stage I", "Stage II")
advanced_stage1 <- c("Stage III", "Stage IVA", "Stage IVB", "Stage IVC")

clin.surv$ajcc_clinical_stage <- ifelse(clin.surv$ajcc_clinical_stage  %in% early_stage1, "Early Stage",
                                        ifelse(clin.surv$ajcc_clinical_stage  %in% advanced_stage1, "Advanced Stage", NA))

as.factor(clin.surv$ajcc_clinical_stage)
str(clin.surv$ajcc_clinical_stage )


###_______________ELIMINAR VALORES 'NA' DA COLUNA 'ajcc_pathologic_n' E MANTER MESMAS AMOSTRAS NA MATRIZ BETA______

# Filtrar linhas com valores não NA em ajcc_pathologic_n
non_na_rows1 <- complete.cases(clin.surv$ajcc_pathologic_n)
clin.surv <- clin.surv[non_na_rows1, ]

# Obter os nomes das linhas não NA em ajcc_pathologic_n
rows_to_keep1 <- rownames(clin.surv)

# Filtrar as colunas de beta_m_surv usando as linhas selecionadas
beta_m_surv <- beta_m_surv[, intersect(colnames(beta_m_surv), rows_to_keep1)] 

as.matrix(beta_m_surv)

# Contar o número de ocorrências de "Nodal Invasion" e "No Invasion" em ajcc_pathologic_n
sum(clin.surv$ajcc_pathologic_n == "Nodal Invasion", na.rm = TRUE)#54
sum(clin.surv$ajcc_pathologic_n == "No Invasion", na.rm = TRUE)#37




###_____________________PAM - ANÁLISE COM DADOS DE INVASÃO LINFONODAL _____________________________

pamr_data2 <- list(x= as.matrix(beta_m_surv),
                   y=factor(clin.surv$ajcc_pathologic_n), 
                   geneid=rownames(beta_m_surv),
                   genenames=paste("cg",as.character(1:nrow(beta_m_surv)),sep=""))



# Realizar a call com a função pamr.train()
train2 <- pamr.train(data = pamr_data2)
train2


## Cross-validate the classifier
results2 <- pamr.cv(train2,pamr_data2) 
results2

## Plot the cross-validated error curves
pamr.plotcv(results2)

## Compute the confusion matrix for a particular model (threshold=0.2) 
pamr.confusion(results2, threshold=0.2)

#                     No Invasion   Nodal Invasion     Class Error rate
#No Invasion             20             17              0.4594595
#Nodal Invasion          18             36              0.3333333
#Overall error rate= 0.381



## Plot the cross-validated class probabilities by class(threshold=0.2) 
pamr.plotcvprob(results2, pamr_data2, threshold=0.2)


## Plot the class centroids
pamr.plotcen(train2, pamr_data2, threshold=0.2)

## Make a gene plot of the most significant genes
pamr.geneplot(train2, pamr_data2, threshold= 0.2)
##Muitos genes, não vai dar pra plotar


# Estimate false discovery rates and plot them
fdr.obj2<- pamr.fdr(train2, pamr_data2)
str(fdr.obj2)
pamr.plotfdr(fdr.obj2)#NULL

##Gerar um gráfico
plot(fdr.obj2$results[, "Threshold"], fdr.obj2$results[, "Median FDR"], type = "l", 
     xlab = "Threshold", ylab = "Median FDR", main = "False Discovery Rates")


## List the significant genes
pamrlist_nodal.ge <- pamr.listgenes(train2, pamr_data2, threshold=0.2)

pamrlist_nodal.ge <- as.data.frame(pamrlist_nodal.ge)

write.csv(pamrlist_nodal.ge, "pamrlist_nodal_ge.csv")




# Obter as sondas encontradas no pam
pam_probes <- pamrlist_nodal.ge$id
print(pam_probes)



###_______________ELIMINAR VALORES 'NA' DA COLUNA 'clin.surv$ajcc_clinical_stage' E MANTER MESMAS AMOSTRAS NA MATRIZ BETA______

# Filtrar linhas com valores não NA em clin.surv$ajcc_clinical_stage
non_na_rows2 <- complete.cases(clin.surv$ajcc_clinical_stage)
clin.surv <- clin.surv[non_na_rows2, ]

# Obter os nomes das linhas não NA em ajcc_pathologic_n
rows_to_keep2 <- rownames(clin.surv)

# Filtrar as colunas de beta_m_surv usando as linhas selecionadas
beta_m_surv <- beta_m_surv[, intersect(colnames(beta_m_surv), rows_to_keep2)] 

as.matrix(beta_m_surv)

# Contar o número de ocorrências de "Early Stage" e "Advanced Stage" em clin.surv$ajcc_clinical_stage
sum(clin.surv$ajcc_clinical_stage == "Early Stage", na.rm = TRUE)#7
sum(clin.surv$ajcc_clinical_stage == "Advanced Stage", na.rm = TRUE)#81




###_____________________PAM - ANÁLISE COM DADOS DE ESTADIAMENTO _____________________________

pamr_data3 <- list(x= as.matrix(beta_m_surv),
                   y=factor(clin.surv$ajcc_clinical_stage), 
                   geneid=rownames(beta_m_surv),
                   genenames=paste("cg",as.character(1:nrow(beta_m_surv)),sep=""))



# Realizar a call com a função pamr.train()
train3 <- pamr.train(data = pamr_data3)
train3


## Cross-validate the classifier
results3 <- pamr.cv(train3,pamr_data3) 
results3

## Plot the cross-validated error curves
pamr.plotcv(results3)

## Compute the confusion matrix for a particular model (threshold=0.2) 
pamr.confusion(results3, threshold=0.1)




## Plot the cross-validated class probabilities by class(threshold=0.2) 
pamr.plotcvprob(results3, pamr_data3, threshold=0.2)


## Plot the class centroids
pamr.plotcen(train3, pamr_data3, threshold=0.2)

## Make a gene plot of the most significant genes
pamr.geneplot(train3, pamr_data3, threshold= 0.2)
##Muitos genes, não vai dar pra plotar


# Estimate false discovery rates and plot them
fdr.obj3<- pamr.fdr(train3, pamr_data3)
str(fdr.obj3)
pamr.plotfdr(fdr.obj3)#NULL

##Gerar um gráfico
plot(fdr.obj3$results[, "Threshold"], fdr.obj3$results[, "Median FDR"], type = "l", 
     xlab = "Threshold", ylab = "Median FDR", main = "False Discovery Rates")


## List the significant genes
pamrlist_stage.ge <- pamr.listgenes(train3, pamr_data3, threshold=0.2)

pamrlist_stage.ge <- as.data.frame(pamrlist_stage.ge)

write.csv(pamrlist_stage.ge, "pamrlist_stage_ge.csv")
