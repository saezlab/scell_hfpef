## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-04-06
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## myofib vs matrifib analysis
# reintegrate and find stable cluster

library(Seurat)
library(tidyverse)
library(harmony)
library(cluster)

source("code/utils.R")


integrated_data = readRDS( "output/seu.objs/study_integrations/harmony_fib_filt.rds")

cell_type = c("2", "3")
class_label = "opt_clust_integrated"

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
  integrated_data <- RunHarmony(integrated_data,group.by.vars = c("orig.ident","study"),
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

    saveRDS(integrated_data, file = paste0("output/seu.objs/myofib_matrifib", ".rds")
  )

# intepret ------------------------------------------------------------------------------------

  integrated_data= readRDS( file = paste0("output/seu.objs/myofib_matrifib", ".rds")
  )

  FeaturePlot(integrated_data, feature= c("Acta2", "Postn", "Cilp", "Comp", "Thbs4"))

  DimPlot(integrated_data, group.by = c("study", "opt_state", "group"), label =  T)

  Idents(integrated_data)

  marks= FindAllMarkers(integrated_data)

  saveRDS(marks, "output/fib_integration/marker_list/myofib.matrifib.subclustermarker.rds")

  FeaturePlot(integrated_data, feature= c("Cd74", "Lyz2", "Acta2", "Mfap4", "Scara5","Cd248"))

  marks %>%as_tibble() %>% filter(cluster== 3) %>% arrange(desc(avg_log2FC)) %>% print(n=40)

  x= integrated_data[[]]
  x= x %>% mutate(group = ifelse(grepl(pattern = "hf", orig.ident), "hf", "ct"))

  integrated_data= AddMetaData(integrated_data, x$group, "group")

  p1= calc_props(integrated_data@meta.data, cluster.col = "study", group.col = "opt_state")
  p2= calc_mean_proportions_per_group(integrated_data@meta.data,
                                      cluster.col = "opt_state",
                                      group.col= "group")


  NABA_SETS_mouse= readRDS("data/prior_knowledge/NABA_mouse.rds")
  integrated_data =add_geneset_scores(integrated_data, NABA_SETS_mouse)

  meta= integrated_data[[]]
  p.deg.up=  ggplot(data= meta, aes(x= opt_state, y= Core_matrisome1,fill =group))+
    scale_fill_manual(values= hfpef_cols)+
    geom_boxplot()+
    theme_minimal()+
    #ggtitle(paste0("top ", cutoff, " genes from DEG" ))+
    ggtitle("upregulated DEG")+
    theme(axis.text.x = element_text(angle= 40, hjust=1),
          axis.text= element_text(color ="black"),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = "none",
          plot.title = element_text(size = 12)
    )
