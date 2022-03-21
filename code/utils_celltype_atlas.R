## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-09-10
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Util functions for cell type specific preprocessing

library(cluster)


process_data= function(seu,
                       dis_filter= F ,
                       default.assay= "RNA"){

  # Initialize the Seurat object

  DefaultAssay(seu) = default.assay
  # Process the data --------------------------------------------------------------------

  seu <- NormalizeData(seu,
                       normalization.method = 'LogNormalize',
                                 scale.factor = 10000,
                                 verbose = FALSE)

  seu <- FindVariableFeatures(seu,
                              selection.method = 'vst',
                              nfeatures = 3000,
                              verbose = FALSE)

  seu <- ScaleData(seu,
                             verbose = FALSE,
                             features = rownames(seu))

  seu <- RunPCA(seu,
                          verbose = FALSE)
  # DimHeatmap(seu, dims = 1:20, cells = 500, balanced = TRUE)
  # JackStrawPlot(seu, dims = 1:10)
  # ElbowPlot(seu, ndims = 30)
  seu <- RunUMAP(seu, reduction = 'pca', dims = 1:30, verbose = FALSE)

  seu@reductions
  ## RE Cluster
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30)

  seq_res <- seq(0.5, 1.5, 0.1)

  seu <- FindClusters(seu,
                                resolution = seq_res,
                                verbose = F)

  DimPlot(seu)
  cell_dists <- dist(seu@reductions$pca@cell.embeddings,
                     method = "euclidean")

  cluster_info <- seu@meta.data[,grepl("RNA_snn_res",
                                                 colnames(seu@meta.data))] %>%
    dplyr::mutate_all(as.character) %>%
    dplyr::mutate_all(as.numeric)

  silhouette_res <- apply(cluster_info, 2, function(x){
    si <- silhouette(x, cell_dists)
    if(!is.na(si)) {
      mean(si[, 'sil_width'])
    } else {
      NA
    }
  })

  seu[["opt_clust_new"]] <- seu[[names(which.max(silhouette_res))]]

  spam_cols <- grepl(paste0(DefaultAssay(seu), "_snn_res"),
                     colnames(seu@meta.data)) |
    grepl("seurat_clusters",colnames(seu@meta.data))

  seu@meta.data <- seu@meta.data[,!spam_cols]

  DimPlot(seu, group.by = "opt_clust_new")

  final_embedding <- DimPlot(seu, group.by = "opt_clust", reduction = "umap_original") +
    ggtitle(paste0("n cells ", ncol(seu)))

  print("Generating outputs")

  # Save data
  Idents(seu) = "opt_clust"

  saveRDS(seu, file = paste0("data/QC_processing/samplewise/seu_processed_", sample_name, ".rds"))

  # Plot QC files

  pdf(file = paste0("output/figures/samplewiseQC/", sample_name, ".pdf"), width = 8, height = 8)

  plot(cowplot::plot_grid(nrow = 1, ncol = 2, filt_p1, filt_p2))
  plot(quality_plt)
  plot(quality_plt_bis)
  #plot(quality_plt_bis_2)
  plot(final_embedding)

  dev.off()

}

