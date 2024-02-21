
# HDCA Heart Visium data analysis 2

# Data import
## All sections 

## Set directory containing folders of spaceranger output with samples included in the analysis
visium.dir <- "/home/st-analysis_home/zaneta.andrusivova/projects/hdca_heart/visium_data/"

## Set path for exporting results
export_path = "/home/st-analysis_home/zaneta.andrusivova/projects/hdca_heart/r_analysis/"

## Load packages
library(Matrix)
library(magrittr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(STutility)
library(Rcpp)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(clustree)
library(harmony)
library(gridExtra)

## Create infoTable and Seurat object
samples <- list.files(visium.dir, recursive = TRUE, full.names = TRUE, pattern = 'filtered_feature_bc_matrix.h5')
spotfiles <- list.files(visium.dir, recursive = TRUE, full.names = TRUE, pattern = 'tissue_positions_list.csv')
imgs <- list.files(visium.dir, recursive = TRUE, full.names = TRUE, pattern = 'tissue_hires_image.png')
json <- list.files(visium.dir, recursive = TRUE, full.names = TRUE, pattern = 'scalefactors_json.json')

section.name <- samples
section.name <- gsub(paste0(visium.dir, "/"),"", gsub("/filtered_feature_bc_matrix.h5","",section.name))

infoTable <- data.frame(section.name, samples, spotfiles, imgs, json, stringsAsFactors = FALSE)

se_heart <- InputFromTable(infotable = infoTable, platform="Visium")


## Data quality control
### number of Genes/UMIs
a <- VlnPlot(se_heart, features = "nFeature_RNA", group.by = "section.name", pt.size = 0) + theme(legend.position = 'none')
b <- VlnPlot(se_heart, features = "nCount_RNA", group.by = "section.name", pt.size = 0) + theme(legend.position = 'none')
a/b

ST.FeaturePlot(se_heart, features = "nFeature_RNA", dark.theme = F, cols = c("lightgrey", "red", "darkred"), pt.size = 1.2, ncol = 4)
ST.FeaturePlot(se_heart, features = "nCount_RNA", ncol = 4, cols = c("lightgrey", "red", "darkred"), pt.size = 1.2)

### hemoglobin presence
c <- VlnPlot(se_heart, features = "HBA1", group.by = "section.name", pt.size = 0) + theme(legend.position = 'none')
d <- ST.FeaturePlot(se_heart, features = "HBA1", ncol = 4, cols = c("lightgrey", "red", "darkred"))
grid.arrange(c, d, ncol=2)

### % of MT and ribo genes
se_heart$percent.mt <- PercentageFeatureSet(se_heart, pattern = "^MT")
se_heart$percent.ribo <- PercentageFeatureSet(se_heart, pattern = "RP")
VlnPlot(se_heart, features = c("percent.mt", "percent.ribo"), ncol = 2, pt.size = 0, cols = ("darkred"))


## Data filtering
### removing low quality spots, MT genes, ribosomal genes and hemoglobin genes
se_heart <- SubsetSTData(se_heart, expression = nFeature_RNA > 200)

enids <- read.table(file = "/home/st-analysis_home/zaneta.andrusivova/projects/hdca_spinalcord/visium_data/GRCh38_genes.tsv", header = T, stringsAsFactors = T)
enids <- data.frame(apply(enids, 2, as.character), stringsAsFactors = F)
rownames(enids) <- enids$gene_id
enids <- subset(enids, gene_biotype %in% c("protein_coding", "lncRNA"))
keep.genes <- intersect(rownames(se_heart), enids$gene_name)
se_heart <- se_heart[keep.genes, ]

mt.genes.heart <- grep(pattern = c("^MT-|^HB|^RP|MALAT1"), x = rownames(se_heart), value = T)
keep.genes.heart1 <- setdiff(rownames(se_heart), mt.genes.heart)
se_heart_filter <- se_heart[keep.genes.heart1, ]

## Data normalization
### SCTrabsform per section

### export staffli object
staffli_heart <- se_heart_filter@tools$Staffli

### split the dataset into a list of objects -> one pre section
se_section_list <- SplitObject(se_heart_filter, split.by = "section.name")

### normalize and identify variable features for each dataset independently
se_section_list <- lapply(se_section_list, SCTransform)

### select features for integration
features <- SelectIntegrationFeatures(object.list = se_section_list, nfeatures = 3000)
se_section_list <- PrepSCTIntegration(object.list = se_section_list, anchor.features = features)

### merge into one object
se_heart_merged <- merge(se_section_list[[1]], y = c(se_section_list[2:length(se_section_list)]), merge.data = T)
VariableFeatures(se_heart_merged) <- features

## PCA
se_norm <- RunPCA(se_heart_merged)

## Section integration with Harmony
se_harmony <- RunHarmony(se_norm, group.by.vars = "section.name", reduction = "pca", assay.use = "SCT")

### add Staffli object to seurat object
se_harmony@tools$Staffli <- staffli_heart

## Gene plots
FeaturePlot(se_harmony, features = c("BMP10"))
ST.FeaturePlot(se_harmony, features ="BMP10", ncol = 5)

## Fig. 1E - sc. 11
feature_list <- c("MT3", "NREP", "CKM", "FHL2", "LGALS3", "MASP1", "RELN", "PITX2", "ANGPT1", "PAM", "COL2A1", "BMP10", "NPPA", "ADAMTS8", "DKK3")
for (feature in feature_list){
  p <-   ST.FeaturePlot(se_harmony, features = feature, indices = 16, pt.size = 2.5, cols = c("white", "turquoise", "blue", "black"))
  pdf(paste0(export_path, "feature_plots/", feature, ".pdf"), width = 10, height = 10)
  print(p)
  dev.off()
}

## Fig. 1F - sc. 2 and 20
feature_list <- c("ACTA1", "NPPB", "ENO1")
for (feature in feature_list){
  p <-   ST.FeaturePlot(se_harmony, features = feature, indices = 2, pt.size = 2.5, cols = c("white", "turquoise", "blue", "black"))
  pdf(paste0(export_path, "feature_plots/", feature, ".pdf"), width = 10, height = 10)
  print(p)
  dev.off()
}
for (feature in feature_list){
  p <-   ST.FeaturePlot(se_harmony, features = feature, indices = 20, pt.size = 2.5, cols = c("white", "turquoise", "blue", "black"))
  pdf(paste0(export_path, "feature_plots/", feature, ".pdf"), width = 10, height = 10)
  print(p)
  dev.off()
}

##Fig. 4E - sc. 16
feature_list <- c("CHGB", "PRPH", "TH", "SLC18A3", "SST", "ADRB1", "CHRM2")
for (feature in feature_list){
  p <-   ST.FeaturePlot(se_harmony, features = feature, indices = 16, pt.size = 2.5, cols = c("white", "turquoise", "blue", "black"))
  pdf(paste0(export_path, "feature_plots/", feature, ".pdf"), width = 10, height = 10)
  print(p)
  dev.off()
}

##Suppl. Fig. 1E - sc. 11
feature_list <- c("MYH7", "MYH6", "HEY2", "MB", "IRX2", "ELN", "COL3A1", "ITLN1", "LYVE1", "ACTA1", "JAG1", "A2M", "IGFBP3", "CDH11", "APCDD1")
for (feature in feature_list){
  p <-   ST.FeaturePlot(se_harmony, features = feature, indices = 11, pt.size = 2.5, cols = c("white", "turquoise", "blue", "black"))
  pdf(paste0(export_path, "feature_plots/", feature, ".pdf"), width = 10, height = 10)
  print(p)
  dev.off()
}

##Suppl. Fig. 3D - sc. 3
feature_list <- c("HS3ST3A1", "RCAN1", "IRX1", "IRX2")
for (feature in feature_list){
  p <-   ST.FeaturePlot(se_harmony, features = feature, indices = 3, pt.size = 2.5, cols = c("white", "turquoise", "blue", "black"))
  pdf(paste0(export_path, "feature_plots/", feature, ".pdf"), width = 10, height = 10)
  print(p)
  dev.off()
}

##Suppl. Fig. 5G - sc. 2, 20, 27
feature_list <- c("APCDD1", "MSX1", "DKK22", "BMP4")
for (feature in feature_list){
  p <-   ST.FeaturePlot(se_harmony, features = feature, indices = 2, pt.size = 2.5, cols = c("white", "turquoise", "blue", "black"))
  pdf(paste0(export_path, "feature_plots/", feature, ".pdf"), width = 10, height = 10)
  print(p)
  dev.off()
}
for (feature in feature_list){
  p <-   ST.FeaturePlot(se_harmony, features = feature, indices = 20, pt.size = 2.5, cols = c("white", "turquoise", "blue", "black"))
  pdf(paste0(export_path, "feature_plots/", feature, ".pdf"), width = 10, height = 10)
  print(p)
  dev.off()
}
for (feature in feature_list){
  p <-   ST.FeaturePlot(se_harmony, features = feature, indices = 27, pt.size = 2.5, cols = c("white", "turquoise", "blue", "black"))
  pdf(paste0(export_path, "feature_plots/", feature, ".pdf"), width = 10, height = 10)
  print(p)
  dev.off()
}
