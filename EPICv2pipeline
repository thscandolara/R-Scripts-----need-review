# ============================================================
# PIPELINE: Análise de Metilação EPICv2 | Câncer de Colo de Útero
# Referencia: https://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
# Autora: TBS
# Data: Abril de 2026
# ============================================================

# Diretórios
baseDir  <- '/scr/tscandolara/genomas_cervix/data/Metiloma'
idatsDir <- '/scr/tscandolara/genomas_cervix/data/Metiloma/idats/'
setwd(baseDir)

# Pacotes
library(limma)
library(minfi)
library(DMRcate)
library(ChAMP)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(enrichplot)
library(openxlsx)

pal <- brewer.pal(8, "Dark2")

# ============================================================
# BLOCO 1 — LEITURA DOS DADOS
# ============================================================

targets <- read.metharray.sheet(idatsDir, pattern = "samplesheet_complete.csv")
targets <- targets[, -1]   # remover coluna extra gerada pelo read.metharray.sheet

# Verificar colunas disponíveis-importante para o SVD depois
colnames(targets)

# Nomes descritivos
targets$ID <- paste(targets$Sample_type, targets$Sample_Name, sep = "_")

rgSet <- read.metharray.exp(targets = targets)
sampleNames(rgSet) <- targets$ID

cat("rgSet carregado:", nrow(rgSet), "probes,", ncol(rgSet), "amostras\n")
#rgSet carregado: 1105209 probes, 152 amostras

# ============================================================
# BLOCO 2 — QC E FILTRAGEM DE AMOSTRAS
# ============================================================

qcReport(rgSet, sampNames = targets$ID,
         sampGroups = targets$Sample_type,
         pdf = "qcReport_Raw.pdf")

detP <- detectionP(rgSet)

# Plot detection p-values
detP_df <- data.frame(
  Sample  = colnames(detP),
  Mean_P  = colMeans(detP),
  Group   = factor(targets$Sample_type)
) 

detP_plot <- ggplot(detP_df, aes(x = Sample, y = Mean_P, fill = Group)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.01, color = "red", linetype = "dashed") +
  labs(title = "QC: Mean Detection P-value", x = NULL,
       y = "Mean detection p-values") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        legend.position = "bottom")
detP_plot

ggsave("qcReport_detP_raw.png", plot = detP_plot, width = 15, height = 7)

# Remover amostras com >10% probes falhando (detP >= 0.01)
pct_fail <- colMeans(detP >= 0.01)

keepSamples  <- colMeans(detP >= 0.01) < 0.1
removedSamples <- colnames(detP)[!keepSamples]

rgSetFlt  <- rgSet[, keepSamples]
detPFlt   <- detP[, keepSamples]
targetsFlt <- targets[keepSamples, ]

message("Amostras removidas: ", length(removedSamples),
        " (", paste(removedSamples, collapse = ", "), ")",
        "\nAmostras restantes: ", ncol(rgSetFlt))

cat("\n% de probes com detP > 0.01 por ammostra (ordenado da maior para a menor):\n")
print(round(sort(pct_fail, decreasing = TRUE)*100, 3))
# Amostras removidas: 1 (TUMOR_M022)
# Amostras restantes: 151


# Identificar probes ruins ANTES de normalizar (lista negra)
# Critério: falha em qualquer amostra restante (mais conservador)
probes_detP_bad <- rownames(detPFlt)[rowMeans(detPFlt > 0.01) > 0]
cat("Probes na lista negra:", length(probes_detP_bad), "\n")
#Probes na lista negra: 127985 


# ============================================================
# BLOCO 3 — SVD PRÉ-NORMALIZAÇÃO (dados brutos)
# ============================================================
# Objetivo: verificar se há batch effect nos dados crus,
# antes de qualquer correção. Usa preprocessRaw para obter
# betas sem normalização.

mSetRaw <- preprocessRaw(rgSetFlt)
beta_raw <- getBeta(mSetRaw)

# Substituir NAs por média da probe (champ.SVD não aceita NAs)
beta_raw_noNA <- t(apply(beta_raw, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

# Montar pd — ajustar os nomes das colunas conforme colnames(targetsFlt)
pd_svd <- data.frame(
  Sample_type    = targetsFlt$Sample_type,
  Relapse_status = targetsFlt$Relapse_status,
  Slide          = targetsFlt$Slide,        # ajuste se necessário
  Array          = targetsFlt$Array   # ajuste se necessário
)
rownames(pd_svd) <- targetsFlt$ID

stopifnot(all(colnames(beta_raw_noNA) == rownames(pd_svd)))

dir.create("SVD_preNorm", showWarnings = FALSE)


champ.SVD(
  beta       = as.data.frame(beta_raw_noNA),
  pd         = pd_svd,
  RGEffect   = FALSE,
  PDFplot    = TRUE,
  Rplot      = TRUE,
  resultsDir = "./SVD_preNorm/"
)
#[SVD analysis will be proceed with 936990 probes and 151 samples.]



# ============================================================
# BLOCO 4 — NORMALIZAÇÃO
# ============================================================
# preprocessFunnorm: recomendado para datasets com diferenças
# globais de metilação (tumor vs normal/lesão).
# Inclui correção de background (noob) + normalização funcional.

mSetNorm <- preprocessFunnorm(rgSetFlt)

cat("Após normalização:", nrow(mSetNorm), "probes,",
    ncol(mSetNorm), "amostras\n")
#Após normalização: 930075 probes, 151 amostras



# ============================================================
# BLOCO 5 — SVD PÓS-NORMALIZAÇÃO (antes de filtrar probes)
# ============================================================
# Objetivo: verificar se a normalização corrigiu batch effects
# técnicos. Usa todas as probes normalizadas para máxima sensibilidade.

beta_norm <- getBeta(mSetNorm)

beta_norm_noNA <- t(apply(beta_norm, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

stopifnot(all(colnames(beta_norm_noNA) == rownames(pd_svd)))

dir.create("SVD_posNorm", showWarnings = FALSE)

champ.SVD(
  beta       = as.data.frame(beta_norm_noNA),
  pd         = pd_svd,
  RGEffect   = FALSE,
  PDFplot    = TRUE,
  Rplot      = TRUE,
  resultsDir = "./SVD_posNorm/"
)

# >> Interpretar os resultados SVD antes de continuar:
# - Se Slide/Array correlacionam com PC1/PC2: batch effect presente.
#   Opção 1: incluir Slide como covariável no modelo (preferencial se
#            a distribuição de grupos biológicos for balanceada por slide).
#   Opção 2: corrigir com champ.runCombat() abaixo.
# - Se Sample_type separa bem: variação biológica está sendo capturada.
# - Anotar qual PC captura Relapse_status (informa poder estatístico).


# [OPCIONAL] Correção de batch com Combat — só rodar se SVD indicar
# batch effect não corrigível pelo modelo:
# mSetNorm_corrected <- champ.runCombat(
#   beta      = beta_norm_noNA,
#   pd        = pd_svd,
#   batchname = "Slide"
# )
# Nesse caso, usar mSetNorm_corrected no lugar de beta_norm para
# extrair M-values: mVals <- minfi::logit2(mSetNorm_corrected)


# ============================================================
# BLOCO 6 — FILTRAGEM DE PROBES
# ============================================================

# 6a. Remover probes com detP alto (lista negra do Bloco 2)
probes_to_remove <- intersect(rownames(mSetNorm), probes_detP_bad)
mSetNormFlt <- mSetNorm[!(rownames(mSetNorm) %in% probes_to_remove), ]
cat("Após remoção por detP:", nrow(mSetNormFlt), "probes\n")
#Após remoção por detP: 802307 probes


# 6b. Remover probes com SNPs no sítio CpG
mSetNormFlt <- dropLociWithSnps(mSetNormFlt)
cat("Após remoção de SNPs:", nrow(mSetNormFlt), "probes\n")
#Após remoção de SNPs: 791457 probes


# 6c. Remover probes cross-reativas (EPICv2 — Pidsley et al. 2024)
manifest2 <- read.csv(
  file.path('/scr/tscandolara/genomas_cervix/data/Metiloma/Augmented_EPICV2_manifest',
            'Augmented EPICv2 manifest_12864_2024_10027_MOESM4_ESM.csv'),
  stringsAsFactors = FALSE
)

xReactive2 <- manifest2$IlmnID[manifest2$CH_WGBS_evidence == "Y"]
keep_xr    <- !(featureNames(mSetNormFlt) %in% xReactive2)
mSetNormFlt <- mSetNormFlt[keep_xr, ]
cat("Após remoção de cross-reativas:", nrow(mSetNormFlt), "probes\n")
#Após remoção de cross-reativas: 780174 probes


# ============================================================
# BLOCO 7 — EXPLORAÇÃO: MDS E DENSITY PLOTS
# ============================================================

# Density: antes vs depois da normalização e filtragem
dev.new()
par(mfrow = c(1, 2))
densityPlot(rgSetFlt, sampGroups = targetsFlt$Sample_type,
            main = "Raw", legend = FALSE)
legend("top", legend = levels(factor(targetsFlt$Sample_type)),
       text.col = pal[1:2])
densityPlot(getBeta(mSetNormFlt), sampGroups = targetsFlt$Sample_type,
            main = "Normalized and Filtered", legend = FALSE)
legend("top", legend = levels(factor(targetsFlt$Sample_type)),
       text.col = pal[1:2])


# MDS — Sample type e Relapse status
dev.new()
par(mfrow = c(1, 2))
plotMDS(getM(mSetNormFlt), top = 1000, gene.selection = "common",
        col = pal[factor(targetsFlt$Sample_type)], cex = 0.8,
        main = "MDS — Sample type")
legend("right", legend = levels(factor(targetsFlt$Sample_type)),
       text.col = pal, cex = 0.7, bg = "white")
plotMDS(getM(mSetNormFlt), top = 1000, gene.selection = "common",
        col = pal[factor(targetsFlt$Relapse_status)], cex = 0.8,
        main = "MDS — Relapse status")
legend("right", legend = levels(factor(targetsFlt$Relapse_status)),
       text.col = pal, cex = 0.7, bg = "white")

# ============================================================
# BLOCO 8 — EXTRAIR M-VALUES E BETA VALUES
# ============================================================

mVals <- getM(mSetNormFlt)
bVals <- getBeta(mSetNormFlt)

cat("Dimensões finais — mVals:", nrow(mVals), "probes x", ncol(mVals), "amostras\n")
#Dimensões finais — mVals: 780174 probes x 151 amostras


# Density plots das matrizes finais
dev.new()
par(mfrow = c(1, 2))
densityPlot(bVals, sampGroups = targetsFlt$Sample_type,
            main = "Beta values", legend = FALSE, xlab = "Beta")
legend("top", legend = levels(factor(targetsFlt$Sample_type)),
       text.col = pal[1:2])
densityPlot(mVals, sampGroups = targetsFlt$Sample_type,
            main = "M-values", legend = FALSE, xlab = "M")
legend("topleft", legend = levels(factor(targetsFlt$Sample_type)),
       text.col = pal[1:2])

# ============================================================
# BLOCO 9 — DMPs: TUMOR vs LESÃO
# ============================================================

sampletype <- factor(targetsFlt$Sample_type, levels = c("LESAO", "TUMOR"))
design_dmp  <- model.matrix(~0 + sampletype)
colnames(design_dmp) <- c("LESAO", "TUMOR")

# >> Se SVD indicou batch effect em Slide: adicionar como covariável
# design_tv <- model.matrix(~0 + sampletype + targetsFlt$Slide)
# colnames(design_tv)[1:2] <- c("LESAO", "TUMOR")

fit_dmp  <- lmFit(mVals, design_dmp)
cont_dmp <- makeContrasts(TUMORvsLESAO = TUMOR - LESAO, levels = design_dmp)
fit2_dmp <- contrasts.fit(fit_dmp, cont_dmp)
fit2_dmp <- eBayes(fit2_dmp)

cat("\nDMPs TUMOR vs LESÃO (FDR < 0.05):\n")
print(summary(decideTests(fit2_dmp)))
# TUMORvsLESAO
# Down         173248
# NotSig       435510
# Up           171416


dmps_TvsL <- topTable(fit2_dmp, number = Inf,
                    coef = "TUMORvsLESAO", adjust.method = "BH")

# Delta-beta
lesao_idx  <- which(sampletype == "LESAO")
tumor_idx  <- which(sampletype == "TUMOR")
beta_lesao <- rowMeans(bVals[, lesao_idx], na.rm = TRUE)
beta_tumor <- rowMeans(bVals[, tumor_idx], na.rm = TRUE)

dmps_TvsL <- cbind(
  dmps_TvsL,
  data.frame(
    Beta_LESAO = beta_lesao[rownames(dmps_TvsL)],
    Beta_TUMOR = beta_tumor[rownames(dmps_TvsL)],
    Delta_Beta = (beta_tumor - beta_lesao)[rownames(dmps_TvsL)]
  )
)

# Anotação com manifest
rownames(manifest2) <- manifest2$IlmnID
manifest_cols <- c("Name", "CHR", "MAPINFO", "UCSC_RefGene_Name",
                   "UCSC_RefGene_Group", "UCSC_CpG_Islands_Name",
                   "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Group",
                   "DNase_Hypersensitivity_NAME")
manifest_flt <- manifest2[rownames(dmps_TvsL), manifest_cols]

clean_genes <- function(x) {
  sapply(strsplit(as.character(x), ";"), function(g) {
    g <- unique(g[g != "" & !is.na(g)])
    if (length(g) == 0) return(NA)
    paste(g, collapse = ";")
  })
}

manifest_flt$UCSC_RefGene_Name <- clean_genes(manifest_flt$UCSC_RefGene_Name)

dmps_TvsL_anotadas <- cbind(dmps_TvsL, manifest_flt)

# Filtragem final
dmps_TvsL_sig <- dmps_TvsL_anotadas[
  dmps_TvsL_anotadas$adj.P.Val < 0.05 & abs(dmps_TvsL_anotadas$Delta_Beta) > 0.2, ]
cat("DMPs significativas (FDR<0.05, |Δβ|>0.2):", nrow(dmps_TvsL_sig), "\n")
#DMPs significativas (FDR<0.05, |Δβ|>0.2): 29642 

write.xlsx(dmps_TvsL_sig, "DMPs_significativas_TUMORvsLESAO.xlsx", asTable = TRUE)


# ============================================================
# BLOCO 10 — DMRs: TUMOR vs LESÃO
# ============================================================

ALLMs_TvsL <- rmSNPandCH(mVals, rmcrosshyb = FALSE)

# Usar targetsFlt$Sample_type diretamente (não gsub nos colnames)
type_TvsL     <- factor(targetsFlt$Sample_type, levels = c("LESAO", "TUMOR"))
design_dmr  <- model.matrix(~type_TvsL)
colnames(design_dmr)

myannotation_tv <- cpg.annotate(
  "array",
  object        = ALLMs_TvsL,
  what          = "M",
  arraytype     = "EPICv2",
  epicv2Filter  = "mean",
  epicv2Remap   = TRUE,
  analysis.type = "differential",
  design        = design_dmr,
  coef          = 2,
  fdr           = 0.05
)
#Your contrast returned 341688 individually significant probes.


dmrcoutput_TvsL  <- dmrcate(myannotation_tv, lambda = 1000, C = 2)
dmrs_TvsL       <- extractRanges(dmrcoutput_TvsL, genome = "hg38")
dmrdf_TvsL       <- as.data.frame(dmrs_TvsL)

cat("DMRs encontradas (TUMOR vs LESÃO):", nrow(dmrdf_TvsL), "\n")
#DMRs encontradas (TUMOR vs LESÃO): 56128 
write.xlsx(dmrdf_TvsL, "DMRs_TUMORvsLESAO.xlsx", asTable = TRUE)

head(dmrdf_TvsL$overlapping.genes)

# ============================================================
# BLOCO DE ENRIQUECIMENTO
# ============================================================

# ---------------------------------------------------------------
# Background — todos os genes cobertos pelas probes testadas
# ---------------------------------------------------------------
background_genes <- unique(trimws(unlist(
  strsplit(na.omit(manifest2$UCSC_RefGene_Name[
    manifest2$IlmnID %in% rownames(mVals)
  ]), ";")
)))
background_genes <- background_genes[background_genes != ""]

background_entrez <- bitr(
  background_genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
background_entrez <- background_entrez[!duplicated(background_entrez$SYMBOL), ]
cat("Background final:", nrow(background_entrez), "genes\n")
#Background final: 25569 genes

# ---------------------------------------------------------------
# Genes das DMPs significativas
# ---------------------------------------------------------------
genes_dmps_TvsL <- unique(trimws(unlist(
  strsplit(na.omit(dmps_TvsL_sig$UCSC_RefGene_Name), ";")
)))
genes_dmps_TvsL <- genes_dmps_TvsL[genes_dmps_TvsL != ""]
cat("Genes DMPs:", length(genes_dmps_TvsL), "\n")
#Genes DMPs: 3996 

genes_dmps_entrez <- bitr(
  genes_dmps_TvsL,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
genes_dmps_entrez <- genes_dmps_entrez[!duplicated(genes_dmps_entrez$SYMBOL), ]
cat("Genes DMPs com mapeamento Entrez:", nrow(genes_dmps_entrez), "\n")
#Genes DMPs com mapeamento Entrez: 3936 

# ---------------------------------------------------------------
# Genes das DMRs
# ---------------------------------------------------------------
genes_dmrs_TvsL <- unique(trimws(unlist(
  strsplit(na.omit(dmrdf_TvsL$overlapping.genes), ", ")
)))
genes_dmrs_TvsL <- genes_dmrs_TvsL[genes_dmrs_TvsL != ""]
cat("Genes DMRs:", length(genes_dmrs_TvsL), "\n")
#Genes DMRs: 23331

genes_dmrs_entrez <- bitr(
  genes_dmrs_TvsL,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
genes_dmrs_entrez <- genes_dmrs_entrez[!duplicated(genes_dmrs_entrez$SYMBOL), ]
cat("Genes DMRs com mapeamento Entrez:", nrow(genes_dmrs_entrez), "\n")
#Genes DMRs com mapeamento Entrez: 17541 

# ============================================================
# ENRIQUECIMENTO — DMPs TUMOR vs LESÃO
# ============================================================

# --- GO BP — DMPs ----------------------------------------------
ego_dmps_TvsL <- enrichGO(
  gene          = genes_dmps_entrez$ENTREZID,
  universe      = background_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
ego_dmps_TvsL <- pairwise_termsim(ego_dmps_TvsL)
cat("GO BP DMPs — termos significativos:", nrow(as.data.frame(ego_dmps_TvsL)), "\n")
#GO BP DMPs — termos significativos: 909

dev.new(); barplot(ego_dmps_TvsL, showCategory = 20, title = "GO BP - DMPs TUMOR vs LESÃO")
dev.new(); dotplot(ego_dmps_TvsL, showCategory = 20, title = "GO BP - DMPs TUMOR vs LESÃO")
dev.new(); emapplot(ego_dmps_TvsL, showCategory = 20)
dev.new(); cnetplot(ego_dmps_TvsL, showCategory = 5)

write.xlsx(as.data.frame(ego_dmps_TvsL), "GO_BP_DMPs_TUMORvsLESAO.xlsx")

# --- KEGG — DMPs -----------------------------------------------
ekegg_dmps_TvsL <- enrichKEGG(
  gene          = genes_dmps_entrez$ENTREZID,
  universe      = background_entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
ekegg_dmps_TvsL <- setReadable(ekegg_dmps_TvsL, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
cat("KEGG DMPs — pathways significativos:", nrow(as.data.frame(ekegg_dmps_TvsL)), "\n")
#KEGG DMPs — pathways significativos: 38 

dev.new(); barplot(ekegg_dmps_TvsL, showCategory = 20, title = "KEGG - DMPs TUMOR vs LESÃO")
dev.new(); dotplot(ekegg_dmps_TvsL, showCategory = 20, title = "KEGG - DMPs TUMOR vs LESÃO")
dev.new(); cnetplot(ekegg_dmps_TvsL, showCategory = 5)

write.xlsx(as.data.frame(ekegg_dmps_TvsL), "KEGG_DMPs_TUMORvsLESAO.xlsx")

# --- Reactome — DMPs -------------------------------------------
ereactome_dmps_TvsL <- enrichPathway(
  gene          = genes_dmps_entrez$ENTREZID,
  universe      = background_entrez$ENTREZID,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
ereactome_dmps_TvsL <- pairwise_termsim(ereactome_dmps_TvsL)
cat("Reactome DMPs — pathways significativos:", nrow(as.data.frame(ereactome_dmps_TvsL)), "\n")

dev.new(); barplot(ereactome_dmps_TvsL, showCategory = 20, title = "Reactome — DMPs TUMOR vs LESÃO")
dev.new(); dotplot(ereactome_dmps_TvsL, showCategory = 20, title = "Reactome — DMPs TUMOR vs LESÃO")
dev.new(); emapplot(ereactome_dmps_TvsL, showCategory = 30)
dev.new(); cnetplot(ereactome_dmps_TvsL, showCategory = 5)

write.xlsx(as.data.frame(ereactome_dmps_TvsL), "Reactome_DMPs_TUMORvsLESAO.xlsx")


# ============================================================
# ENRIQUECIMENTO — DMRs TUMOR vs LESÃO
# ============================================================

# --- GO BP — DMRs ----------------------------------------------
ego_dmrs_TvsL <- enrichGO(
  gene          = genes_dmrs_entrez$ENTREZID,
  universe      = background_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
ego_dmrs_TvsL <- pairwise_termsim(ego_dmrs_TvsL)
cat("GO BP DMRs — termos significativos:", nrow(as.data.frame(ego_dmrs_TvsL)), "\n")
#GO BP DMRs — termos significativos: 864 

dev.new(); barplot(ego_dmrs_TvsL, showCategory = 20, title = "GO BP - DMRs TUMOR vs LESÃO")
dev.new(); dotplot(ego_dmrs_TvsL, showCategory = 20, title = "GO BP - DMRs TUMOR vs LESÃO")
write.xlsx(as.data.frame(ego_dmrs_TvsL), "GO_BP_DMRs_TUMORvsLESAO.xlsx")

# --- KEGG — DMRs -----------------------------------------------
ekegg_dmrs_TvsL <- enrichKEGG(
  gene          = genes_dmrs_entrez$ENTREZID,
  universe      = background_entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
ekegg_dmrs_TvsL <- setReadable(ekegg_dmrs_TvsL, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
cat("KEGG DMRs — pathways significativos:", nrow(as.data.frame(ekegg_dmrs_TvsL)), "\n")
#KEGG DMRs — pathways significativos: 72 

dev.new(); barplot(ekegg_dmrs_TvsL, showCategory = 20, title = "KEGG - DMRs TUMOR vs LESÃO")
dev.new(); dotplot(ekegg_dmrs_TvsL, showCategory = 20, title = "KEGG - DMRs TUMOR vs LESÃO")
write.xlsx(as.data.frame(ekegg_dmrs_TvsL), "KEGG_DMRs_TUMORvsLESAO.xlsx")

# --- Reactome — DMRs -------------------------------------------
ereactome_dmrs_TvsL <- enrichPathway(
  gene          = genes_dmrs_entrez$ENTREZID,
  universe      = background_entrez$ENTREZID,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
ereactome_dmrs_TvsL <- pairwise_termsim(ereactome_dmrs_TvsL)
cat("Reactome DMRs — pathways significativos:", nrow(as.data.frame(ereactome_dmrs_TvsL)), "\n")
#Reactome DMRs — pathways significativos: 71

dev.new(); barplot(ereactome_dmrs_TvsL, showCategory = 20, title = "Reactome - DMRs TUMOR vs LESÃO")
dev.new(); dotplot(ereactome_dmrs_TvsL, showCategory = 20, title = "Reactome - DMRs TUMOR vs LESÃO")
dev.new(); emapplot(ereactome_dmrs_TvsL, showCategory = 10)
dev.new(); cnetplot(ereactome_dmrs_TvsL, showCategory = 5)

write.xlsx(as.data.frame(ereactome_dmrs_TvsL), "Reactome_DMRs_TUMORvsLESAO.xlsx")

# ============================================================
# BLOCO 12 — DMPs: RECIDIVA (só TUMOR)
# ============================================================

cat("\nDistribuição de Relapse_status nos tumores:\n")
print(table(targetsFlt$Sample_type, targetsFlt$Relapse_status))

idx_tumor    <- targetsFlt$Sample_type == "TUMOR" &
  !is.na(targetsFlt$Relapse_status) &
  targetsFlt$Relapse_status != ""

targetsTumor <- targetsFlt[idx_tumor, ]
mValsTumor   <- mVals[, idx_tumor]
bValsTumor   <- bVals[, idx_tumor]

cat("Tumores para análise de recidiva:", nrow(targetsTumor), "\n")
print(table(targetsTumor$Relapse_status))

relapse <- factor(targetsTumor$Relapse_status)
cat("Levels de relapse:", levels(relapse), "\n")

design_rel          <- model.matrix(~0 + relapse)
colnames(design_rel) <- levels(relapse)  # "NO" "YES"
colnames(design_rel)

fit_rel  <- lmFit(mValsTumor, design_rel)

cont_rel <- makeContrasts(
  YESvsNO = YES - NO,
  levels  = design_rel
)
fit2_rel <- contrasts.fit(fit_rel, cont_rel)
fit2_rel <- eBayes(fit2_rel)

cat("\nDMPs Recidiva (FDR < 0.05):\n")
print(summary(decideTests(fit2_rel)))
# DMPs Recidiva (FDR < 0.05):
# YESvsNO
# Down         0
# NotSig  780174
# Up           0


# coef = "YESvsNO" 
dmps_rel <- topTable(fit2_rel, number = Inf,
                     coef = "YESvsNO", adjust.method = "BH")

# Ver quantos CpGs ficam próximos de significância
cat("p < 0.05 (sem correção):", sum(dmps_rel$P.Value < 0.05), "\n") #p < 0.05 (sem correção): 55278 
cat("p < 0.01 (sem correção):", sum(dmps_rel$P.Value < 0.01), "\n")#p < 0.01 (sem correção): 11042
cat("FDR < 0.10:", sum(dmps_rel$adj.P.Val < 0.10), "\n")
cat("FDR < 0.20:", sum(dmps_rel$adj.P.Val < 0.20), "\n")

# Ver distribuição dos p-values brutos (deve ser uniforme se não há sinal)
hist(dmps_rel$P.Value, breaks = 100,
     main = "Distribuição p-values - Status de Recidiva",
     xlab = "p-value")





rel_idx   <- which(relapse == "YES")
norel_idx <- which(relapse == "NO")
beta_rel   <- rowMeans(bValsTumor[, rel_idx],   na.rm = TRUE)
beta_norel <- rowMeans(bValsTumor[, norel_idx], na.rm = TRUE)

dmps_rel <- cbind(
  dmps_rel,
  data.frame(
    Beta_YES   = beta_rel[rownames(dmps_rel)],
    Beta_NO    = beta_norel[rownames(dmps_rel)],
    Delta_Beta = (beta_rel - beta_norel)[rownames(dmps_rel)]
  )
)

manifest_rel                   <- manifest2[rownames(dmps_rel), manifest_cols]
manifest_rel$UCSC_RefGene_Name <- clean_genes(manifest_rel$UCSC_RefGene_Name)
dmps_rel_anotadas              <- cbind(dmps_rel, manifest_rel)

dmps_rel_sig <- dmps_rel_anotadas[
  dmps_rel_anotadas$adj.P.Val < 0.05 &
    abs(dmps_rel_anotadas$Delta_Beta) > 0.2, ]

cat("DMPs recidiva significativas:", nrow(dmps_rel_sig), "\n")
write.xlsx(dmps_rel_sig, "DMPs_Relapse.xlsx", asTable = TRUE)

# ============================================================
# BLOCO 13 — DMRs: RECIDIVA (só TUMOR)
# ============================================================

ALLMs_rel      <- rmSNPandCH(mValsTumor, rmcrosshyb = FALSE)
design_dmr_rel <- model.matrix(~relapse)
design_dmr_rel
# Com ~relapse e NO como referência, coef 2 = YES vs NO
colnames(design_dmr_rel)[2] <- "YES"
colnames(design_dmr_rel)

myannotation_rel <- cpg.annotate(
  "array",
  object        = ALLMs_rel,
  what          = "M",
  arraytype     = "EPICv2",
  epicv2Filter  = "mean",
  epicv2Remap   = TRUE,
  analysis.type = "differential",
  design        = design_dmr_rel,
  coef          = 2,
  fdr           = 0.05
)

dmrcoutput_rel <- dmrcate(myannotation_rel, lambda = 1000, C = 2)
dmrs_rel       <- extractRanges(dmrcoutput_rel, genome = "hg38")
dmrdf_rel      <- as.data.frame(dmrs_rel)

cat("DMRs recidiva encontradas:", nrow(dmrdf_rel), "\n")
write.xlsx(dmrdf_rel, "DMRs_Relapse.xlsx", asTable = TRUE)

# ============================================================
# BLOCO 14 — ENRIQUECIMENTO: DMPs E DMRs RECIDIVA
# ============================================================

# --- Preparar genes das DMPs de recidiva ---
sigCpGs_rel <- rownames(dmps_rel_sig)
allCpGs_rel <- rownames(dmps_rel)

# --- Preparar genes das DMRs de recidiva ---
genes_dmrs_rel <- unique(trimws(unlist(
  strsplit(na.omit(dmrdf_rel$overlapping.genes), ", ")
)))
genes_dmrs_rel <- genes_dmrs_rel[genes_dmrs_rel != ""]
cat("Genes DMRs recidiva:", length(genes_dmrs_rel), "\n")

genes_dmrs_rel_entrez <- bitr(
  genes_dmrs_rel,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
genes_dmrs_rel_entrez <- genes_dmrs_rel_entrez[
  !duplicated(genes_dmrs_rel_entrez$SYMBOL), ]
cat("Genes DMRs recidiva com Entrez:", nrow(genes_dmrs_rel_entrez), "\n")

# --- gometh: DMPs recidiva — GO e KEGG ---
# CORRIGIDO: array.type = "EPICv2"
gst_rel_GO <- gometh(
  sig.cpg    = sigCpGs_rel,
  all.cpg    = allCpGs_rel,
  collection = "GO",
  array.type = "EPICv2"
)
cat("\nTop 10 GO — Recidiva (DMPs):\n")
print(topGSA(gst_rel_GO, number = 10))
write.csv(gst_rel_GO[gst_rel_GO$FDR < 0.05, ],
          "Enrichment_DMPs_GO_Relapse.csv", row.names = FALSE)

gst_rel_KEGG <- gometh(
  sig.cpg    = sigCpGs_rel,
  all.cpg    = allCpGs_rel,
  collection = "KEGG",
  array.type = "EPICv2"   # CORRIGIDO
)
cat("\nTop 10 KEGG — Recidiva (DMPs):\n")
print(topGSA(gst_rel_KEGG, number = 10))
write.csv(gst_rel_KEGG[gst_rel_KEGG$FDR < 0.05, ],
          "Enrichment_DMPs_KEGG_Relapse.csv", row.names = FALSE)

# --- enrichGO / KEGG / Reactome: DMRs recidiva ---
ego_dmrs_rel <- enrichGO(
  gene          = genes_dmrs_rel_entrez$ENTREZID,
  universe      = background_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
ego_dmrs_rel <- pairwise_termsim(ego_dmrs_rel)
cat("GO BP DMRs recidiva — termos significativos:",
    nrow(as.data.frame(ego_dmrs_rel)), "\n")

dev.new(); dotplot(ego_dmrs_rel, showCategory = 20,
                   title = "GO BP - DMRs Recidiva")
dev.new(); emapplot(ego_dmrs_rel, showCategory = 20)
dev.new(); cnetplot(ego_dmrs_rel, showCategory = 5)
write.xlsx(as.data.frame(ego_dmrs_rel),
           "GO_BP_DMRs_Relapse.xlsx")

ekegg_dmrs_rel <- enrichKEGG(
  gene          = genes_dmrs_rel_entrez$ENTREZID,
  universe      = background_entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
ekegg_dmrs_rel <- setReadable(ekegg_dmrs_rel,
                              OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
cat("KEGG DMRs recidiva — pathways significativos:",
    nrow(as.data.frame(ekegg_dmrs_rel)), "\n")

dev.new(); dotplot(ekegg_dmrs_rel, showCategory = 20,
                   title = "KEGG - DMRs Recidiva")
dev.new(); cnetplot(ekegg_dmrs_rel, showCategory = 5)
write.xlsx(as.data.frame(ekegg_dmrs_rel),
           "KEGG_DMRs_Relapse.xlsx")

ereactome_dmrs_rel <- enrichPathway(
  gene          = genes_dmrs_rel_entrez$ENTREZID,
  universe      = background_entrez$ENTREZID,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
ereactome_dmrs_rel <- pairwise_termsim(ereactome_dmrs_rel)
cat("Reactome DMRs recidiva — pathways significativos:",
    nrow(as.data.frame(ereactome_dmrs_rel)), "\n")

dev.new(); dotplot(ereactome_dmrs_rel, showCategory = 20,
                   title = "Reactome - DMRs Recidiva")
dev.new(); emapplot(ereactome_dmrs_rel, showCategory = 20)
dev.new(); cnetplot(ereactome_dmrs_rel, showCategory = 5)
write.xlsx(as.data.frame(ereactome_dmrs_rel),
           "Reactome_DMRs_Relapse.xlsx")

# ============================================================
sessionInfo()
