
# HDCA Heart Visium data analysis

# Sample selection and data import
## Heart sections with 3 and 4 chambers imported for clustering analysis

## Set directory containing folders of spaceranger output with samples included in the analysis
visium.dir <- "/home/st-analysis_home/zaneta.andrusivova/projects/hdca_heart/visium_data_selection/smaller_selection/"

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

infoTable <- data.frame(section.name, samples, spotfiles, imgs, json, stringsAsFactors = FALSE, 
                        age = c("w8", "w8", "w9", "w9",
                                "w6", "w6", "w7", "w7",
                                "w10", "w10", "w12", "w8", "w8",
                                "w8", "w8", "w10.5", "w10.5"))

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
ElbowPlot(se_norm, ndims = 50)

## Section integration with Harmony
se_harmony <- RunHarmony(se_norm, group.by.vars = "section.name", reduction = "pca", assay.use = "SCT")

### add Staffli object to seurat object
se_harmony@tools$Staffli <- staffli_heart

## Clustering
se_harmony <- FindNeighbors(se_harmony, reduction = "harmony", dims = 1:30)
se_harmony <- FindClusters(se_harmony)
se_harmony <- RunUMAP(se_harmony, dims = 1:30, reduction = "harmony")

DimPlot(se_harmony, reduction = "umap", dims = 1:2, pt.size = 1, label = T)
ST.FeaturePlot(se_harmony, features = "seurat_clusters", pt.size = 1.2, ncol = 4)

k <- DimPlot(se_harmony, group.by = "section.name", pt.size = 1)
l <- ST.FeaturePlot(se_harmony, features = "seurat_clusters", pt.size = 1.2, ncol = 4)
grid.arrange(k, l, ncol=2)

## DE Analysis
se_harmony <- PrepSCTFindMarkers(se_harmony)
de.markers.heart <- FindAllMarkers(se_harmony, only.pos = TRUE, logfc.threshold = 0.5)

### Rename clusters
annotation <- read.csv2("~/metadata/HDCA_heart_ST_annotations.csv", header = 1, sep = ";")
cells_clusters <- setNames(annotation$cluster_name, nm = annotation$cluster)
de.markers.heart$cluster <- cells_clusters[as.character(de.markers.heart$cluster)]

### Save the detable as a csv file. (Supp Tab. 1)
write.csv2(de.markers.heart, "~/supplementary_tables/supplementary_table_1_DEG_ST.csv")

top10 <- de.markers.heart %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.01) %>%
  filter(row_number() %in% 1:8)

DotPlot(se_harmony, features = unique(top10$gene) %>% rev()) + coord_flip()

## Gene plots
FeaturePlot(se_harmony, features = c("BRINP3", "PLXNA4"))
ST.FeaturePlot(se_harmony, features = c("BRINP3", "PLXNA4"), dark.theme = F, split.labels = T, indices = 1, show.sb = FALSE, ncol = 4)


### Save clustering object
saveRDS(se_harmony, "~data/ST/HDCA_heart_clustering.rds")

