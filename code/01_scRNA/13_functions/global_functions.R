CreateNNAdjacencyMatrix <- function (
  data, 
  reduction = "pca",
  k = 10,
  dims = 1:50
) {
  stopifnot(inherits(data, what = "Seurat"),
            reduction %in% names(data@reductions),
            max(dims) <= ncol(data[[reduction]]@cell.embeddings))
  knn <- RcppHNSW::hnsw_knn(data[[reduction]]@cell.embeddings[ ,dims], 
                            k = k, 
                            distance = "cosine")
  i <- rep(1:nrow(knn$idx), ncol(knn$idx))
  j <- c(knn$idx)
  p <- c(knn$dist>0)*1
  nn <- as(sparseMatrix(i = i, j=j, x = p, dims = c(nrow(knn$idx),nrow(knn$idx)),
                        dimnames = list(colnames(data), colnames(data))), "dgCMatrix")
  return(nn)
}


LouvainClustering <- function (
  nn,
  seed = 1,
  resolution = 1
) {
  g <- igraph::graph_from_adjacency_matrix(nn, mode = "undirected")
  set.seed(seed)
  cl <- igraph::cluster_louvain(g, resolution = resolution)
  return(cl$membership)
}


QC_03_plot <- function (
  data,
  filename,
  to_filter
) {
  png(filename = filename, width = 1200*5, height = 1200*4, res = 300)
  par(mfrow = c(4, 4), mar = c(4, 4, 4, 4))
  hist(data@meta.data$percent_mito, breaks = 400)
  hist(data@meta.data$percent_ribo, breaks = 400)
  hist(data$percent_hb, breaks = 100)
  hist(data$percent_hsp, breaks = 100)
  hist(data@meta.data$nCount_RNA, breaks = 400)
  hist(data@meta.data$nFeature_RNA, breaks = 400)
  # hist(data@meta.data$DoubletScore, breaks = 400)
  plot(data@meta.data$percent_mito, data@meta.data$percent_ribo, cex = 0.2, pch = 16, col = ifelse(to_filter, "red", "black"))
  plot(data@meta.data$nCount_RNA, data@meta.data$nFeature_RNA, cex = 0.2, pch = 16, col = ifelse(to_filter, "red", "black"))
  plot(log(data@meta.data$nCount_RNA), log(data@meta.data$nFeature_RNA), cex = 0.2, pch = 16, col = ifelse(to_filter, "red", "black"))
  plot(data$percent_mito, data$percent_hb, log = "xy", cex = .1, col = ifelse(to_filter, "red", "black"))
  abline(v=30, h=10)
  plot(data$percent_ribo, data$percent_hb, log = "xy", cex = .1, col = ifelse(to_filter, "red", "black"))
  abline(v=3, h=10)
  plot(data$nFeature_RNA, data$percent_hb, log = "xy", cex = .1, col = ifelse(to_filter, "red", "black"))
  abline(v=250, h=10)
  plot(data$percent_hsp, data$percent_mito, cex = .1, col = ifelse(to_filter, "red", "black"))
  plot(data$percent_hsp, data$percent_ribo, cex = .1, col = ifelse(to_filter, "red", "black"))
  dev.off()
}


QC_04_plot <- function (
  filename,
  features_num,
  to_filter
) {
  png(filename, width = 1200*4, height = 1200*2, res = 300)
  par(mfrow=c(1,2), mar=c(7,10,1,1))
  violist(data, genes = features_num, clustering = "sampleID", transparency = 50, srt = 45, pt.col = ifelse(to_filter, "red", "black"), pt.cex = 0.1)
  barlist(data, genes = features_num, clustering = "sampleID", srt = 45, orderby = "percent_mito", col = ifelse(to_filter, "red", "black"))
  dev.off()
}



FindDoublets <- function (
  data,
  myfolder
) {
  stopifnot("sampleID" %in% colnames(data[[]]))
  # Doublet Finding for each dataset separately
  datasets <- unique(data$sampleID)
  dataset_size <- table(data$sampleID)
  
  stopifnot("sampleID" %in% colnames(data[[]]))
  cli::cli_alert_info("Splitting Seurat object by sampleID")
  data_list <- lapply(
    datasets, function(x) {
      cli::cli_alert("Creating a subset for sample '{x}'")
      subset <- data[, data@meta.data$sampleID == x]
      cli::cli_alert("Normalizing data")
      subset <- NormalizeData(subset, verbose = FALSE)
      cli::cli_alert("Calculating top variable features")
      subset <- FindVariableFeatures(subset, nfeatures = 4000, verbose = FALSE) # Could be reduced to 2000 to save time
      cli::cli_alert("Scaling data")
      subset <- ScaleData(subset, verbose = FALSE)
      cli::cli_alert("Running PCA")
      subset <- RunPCA(subset, npcs = 50, verbose = FALSE, seed.use = 42)
      # knn <- RcppHNSW::hnsw_knn(subset@reductions$pca@cell.embeddings[,1:50], 
      #                           k = 10, 
      #                           distance = "cosine")
      # i <- rep(1:nrow(knn$idx), ncol(knn$idx))
      # j <- c(knn$idx)
      # p <- c(knn$dist>0)*1
      # nn <- as(sparseMatrix(i = i, j=j, x=p, dims = c(nrow(knn$idx),nrow(knn$idx)),
      #                       dimnames = list(colnames(subset), colnames(subset))), "dgCMatrix")
      cli::cli_alert("Creating adjacency matrix from PCA")
      nn <- CreateNNAdjacencyMatrix(subset, reduction = "pca", dims = 1:50)
      # g <- igraph::graph_from_adjacency_matrix(nn,mode = "undirected")
      # set.seed(1)
      # cl <- igraph::cluster_louvain(g, resolution = 1)
      cli::cli_alert("Clustering of adjacency matrix")
      cl <- LouvainClustering(nn, resolution = 1)
      subset$clusters_louvain_dataset <- factor(cl)
      cli::cli_alert("Running UMAP")
      subset <- RunUMAP(subset, 
                        dims = 1:50, 
                        reduction = "pca",
                        n.neighbors = 10,
                        verbose = FALSE,
                        n.epochs = 200,
                        reduction.key = "unfilt",
                        reduction.name = "umap_unfilt",
                        seed.use = 42)
      png(filename = paste0(myfolder, "QC/Doublet_Finding/dataset_clusters_optimized_individual_datasets_prefiltered_adjusted/", x, "_clusters.png"), width = 2500, height = 2500, res = 300)
      print(DimPlot(subset, group.by = "clusters_louvain_dataset", reduction = "umap_unfilt"))
      dev.off()
      cli::cli_alert("Performing pN-pK parameter sweeps (DoubletFinder)")
      sweep.res.list <- paramSweep_v3(subset, PCs = 1:50, sct = FALSE)
      cli::cli_alert("Summarizing sweeps (DoubletFinder)")
      sweep.stats <- summarizeSweep(sweep.res.list , GT = FALSE)
      bcmvn <- find.pK(sweep.stats)
      pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK) # Could look more into this line and the one below
      pK <- as.numeric(as.character(pK[[1]]))
      cli::cli_alert("Selected pK: {pK}")
      ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
        geom_point() +
        geom_line()
      ggsave(paste0(myfolder, "QC/Doublet_Finding/pK_estimation_optimized_individual_datasets_prefiltered_adjusted/", x, "_pK.png"), width = 12, height = 4)
      annotations <- subset$clusters_louvain_dataset
      cli::cli_alert("Model proportions of homotypic doublets (DoubletFinder)")
      homotypic.prop <- modelHomotypic(annotations) # Percentage of doublets in this dataset that can be expected to be homotypic rather than heterotypic
      # nExp <- (nrow(subset@meta.data))/130000
      nExp <- dataset_size[x]/125000 # Percentage of expected doublets in the whole dataset before filtering
      cli::cli_alert("Percentage of expected doublets: {nExp} (DoubletFinder)")
      nExp_poi <- round(nExp*dataset_size[x]) # Number of expected doublets in the whole dataset before filtering
      nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # Number of expected heterotypic doublets (which is what we are predicting here) in the whole dataset before filtering.
      n_filtered <- dataset_size[x] - nrow(subset@meta.data) # Number of cells already filtered out for this dataset
      prop_filtered <- n_filtered/dataset_size[x] # Proportion of dataset already removed by filtering
      n_doublets_filtered <- round(prop_filtered*nExp_poi.adj) # Number of doublets already expected to be filtered out by chance
      nExp_poi.adj.adj <- nExp_poi.adj - n_doublets_filtered # Number of doublets still to be filtered out from this dataset
      cli::cli_alert("Number of doublets to be filtered out: {nExp_poi.adj.adj} (DoubletFinder)")
      cli::cli_alert("Model proportions of homotypic doublets (DoubletFinder)")
      cli::cli_alert("Running DoubletFinder with selected parameters (DoubletFinder)")
      subset <- doubletFinder_v3(subset,
                                 pN = 0.25,
                                 pK = pK,
                                 nExp = nExp_poi.adj.adj,
                                 PCs = 1:50)
      cli::cli_alert("Saving results to sample {x} Seurat object")
      subset@meta.data$pANN <- subset@meta.data[,length(colnames(subset@meta.data)) - 1]
      subset@meta.data[,length(colnames(subset@meta.data))-2] <- NULL
      subset@meta.data$DF.classifications <- subset@meta.data[,length(colnames(subset@meta.data)) - 1]
      subset@meta.data[,length(colnames(subset@meta.data))-2] <- NULL
      return(subset)
    }
  )
  #cli::cli_alert_info("Merge Seurat objects")
  #data <- merge(data_list[[1]], data_list[-1])
  #rm(data_list)
  #return(data)
  return(data_list)
}
