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
## Subset to fibros ,reintegrate and find stable cluster

library(Seurat)
library(tidyverse)
library(harmony)
library(cluster)

source("analysis/utils.R")

seu = readRDS( file = "output/seu.objs/integrated_cellstate_noECs.rds")

# Function to process a single celltype into cell states --------------------------------------

process_cell_type= function(integrated_data,
                            cell_type= "T.cells",
                            class_label= "celltype",
                            label= "X",
                            start_res= 0.4){

  if(label=="X"){label= cell_type}
  # Cell types to include ----------------------------------------------------------------------------

  # Read object to subset -----------------------------------------------------

  true_ix <- integrated_data@meta.data[, class_label] %in% cell_type

  # subset based on filtering and quickly get characteristic profile
  integrated_data <- integrated_data[ , true_ix]

  # subset the object in 2 main categories: batch and patient

  integrated_data <- SplitObject(integrated_data, split.by = "orig.ident")

  # First get all the most variable genes per patient per batch

  hvg_list<- map(integrated_data, function(x) {

      DefaultAssay(x) <- "RNA"

      x <- FindVariableFeatures(x, selection.method = "vst",
                                nfeatures = 3000)

      x@assays[["RNA"]]@var.features

    }) %>% unlist()

  hvg_list <- table(hvg_list) %>%
      sort(decreasing = TRUE)

  gene_selection <- hvg_list[1:2000] %>% names()

  #   # Genes that are stable among batches
  # gene_selection <- hvg_list[hvg_list == 2] %>% names()
  #
  # Re-integrate data

  integrated_data <- reduce(integrated_data,
                            merge,
                            merge.data = TRUE)

  integrated_data <- integrated_data %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(features = gene_selection,
           #npcs = 400,
           verbose = FALSE)

  # Original PCs
  original_pca_plt <- DimPlot(object = integrated_data,
                              reduction = "pca",
                              pt.size = 1,
                              group.by = "orig.ident")


  # Integrate the data -----------------------
  integrated_data <- RunHarmony(integrated_data,
                                "orig.ident",
                                plot_convergence = TRUE,
                                max.iter.harmony = 20)

  # Corrected dimensions -----------------------
  corrected_pca_plt <- DimPlot(object = integrated_data,
                               reduction = "harmony",
                               pt.size =1,
                               group.by = "orig.ident")

  # Create the UMAP with new reduction -----------
  integrated_data <- integrated_data %>%
    RunUMAP(reduction = "harmony", dims = 1:30)

  # Clustering and optimization -------------------------
  print("Optimizing clustering")

  integrated_data <- FindNeighbors(integrated_data,
                                   reduction = "harmony",
                                   dims = 1:30)

  seq_res <- seq(start_res, 1, 0.1)

  # Delete previous clustering
  integrated_data@meta.data <- integrated_data@meta.data[, !grepl("RNA_snn_res",
                                                                  colnames(integrated_data@meta.data))]

  # Create new clustering
  integrated_data <- FindClusters(integrated_data,
                                  resolution = seq_res,
                                  verbose = F)

  # Optimize clustering ------------------------------------------------------
  cell_dists <- dist(integrated_data@reductions$harmony@cell.embeddings,
                     method = "euclidean")

  cluster_info <- integrated_data@meta.data[,grepl("RNA_snn_res",
                                                   colnames(integrated_data@meta.data))] %>%
    dplyr::mutate_all(as.character) %>%
    dplyr::mutate_all(as.numeric)

  silhouette_res <- apply(cluster_info, 2, function(x){
    si <- silhouette(x, cell_dists)
    mean(si[, 'sil_width'])
  })

  integrated_data[["opt_state"]] <- integrated_data[[names(which.max(silhouette_res))]]

  Idents(integrated_data) = "opt_state"

  # Delete others clustering
  integrated_data@meta.data <- integrated_data@meta.data[, !grepl("RNA_snn_res",
                                                                  colnames(integrated_data@meta.data))]

   # Save object ------------------------------------------------------

  saveRDS(integrated_data, file = paste0("output/cell_specific_obj/", label, ".rds")
          )

  # Print QC file ------------------------------------------------------

  umap_corrected_sample <- DimPlot(object = integrated_data,
                                   reduction = "umap",
                                   pt.size = 1,
                                   group.by = "orig.ident")

  umap_corrected_clustering <- DimPlot(object = integrated_data,
                                       reduction = "umap",
                                       pt.size = 1,
                                       group.by = "opt_state")

  # Print proportions ------------------------------------------------------
  x= integrated_data[[]]
  x= x %>% mutate(group = ifelse(grepl(pattern = "hf", orig.ident), "hf", "ct"))
  integrated_data= AddMetaData(integrated_data, x$group, "group")

  p1= calc_props(integrated_data@meta.data, cluster.col = "orig.ident", group.col = "opt_state")
  p2= calc_mean_proportions_per_group(integrated_data@meta.data,
                                      cluster.col = "opt_state",
                                      group.col= "group")

  p2= plot_mean_proportions(p2$groupwise, cluster.col = "opt_state", main = "opt_state mean percentage")

  cell_count_sample = rownames_to_column(x, "cellid") %>%
     group_by(orig.ident) %>%
     count(cellid) %>%
     ggplot(aes(x = n,
                y = orig.ident)) +
     geom_bar(stat="identity") +
     theme(legend.position = "bottom") +
     ylab("") + xlab("n_cells")

  cell_count_group = rownames_to_column(x, "cellid") %>%
    group_by(group) %>%
    count(cellid) %>%
    ggplot(aes(x = n,
               y = group)) +
    geom_bar(stat="identity") +
    theme(legend.position = "bottom") +
    ylab("") + xlab("n_cells")

  pdf(file = paste0("output/cell_specific_obj/",label, ".pdf"), height = 7, width = 9)

  print(original_pca_plt)
  print(corrected_pca_plt)
  print(umap_corrected_sample)
  print(umap_corrected_clustering)
  print(cell_count_sample)
  print(cell_count_group)
  print(p1)
  print(p2)

  dev.off()

}


# call differnet cell types: ------------------------------------------------------------------

#cell_class
seu.meta= seu[[]]
celltypes= unique(seu.meta$celltype)

map(celltypes, function(x){
  print(x)
  process_cell_type(seu, x)})

# run by cell type:

process_cell_type(seu, cell_type = c("T.cells", "NK.cells", "B.cells") , label= "lymphocytes")
process_cell_type(seu, cell_type = c("macrophages", "granulocytes" ) , label= "myeloid.cells")

process_cell_type(seu, cell_type = c("Endothelial"))
process_cell_type(seu, cell_type = "fibroblasts")
process_cell_type(seu, cell_type = "macrophages")
process_cell_type(seu, cell_type = "B.cells")
process_cell_type(seu, cell_type= "T.cells")
