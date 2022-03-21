## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-09-06
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Take single file ouput form run_sample_wise_preprocessing.R and integrate samples
## via harmony!

library(tidyverse)
library(Seurat)
library(harmony)
library(cluster)
library(clustree)
library(cowplot)
library(purrr)


path= "output/seu.objs/samplewise/"
slide_files= list.files("output/seu.objs/samplewise/")
def_assay= "RNA"

# Because of incompatibility with objects objects should be appended manually
slide_files_path <- set_names(paste0(path,slide_files), gsub(pattern = "[.]rds",
                                                             replacement = "",
                                                             slide_files))

integrated_data <- map(slide_files_path, readRDS)

print("You managed to load everything")

# Calculate HVG per sample - Here we assume that batch and patient effects aren't as strong
# since cell-types and niches should be greater than the number of batches
hvg_list <- map(integrated_data, function(x) {

  DefaultAssay(x) <- def_assay

  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 3000)

  x@assays[[def_assay]]@var.features

}) %>% unlist()

hvg_list <- table(hvg_list) %>%
  sort(decreasing = TRUE)

gene_selection_plt <- hvg_list %>% enframe() %>%
  group_by(value) %>%
  mutate(value = as.numeric(value)) %>%
  summarize(ngenes = length(name)) %>%
  ggplot(aes(x = value, y = ngenes)) +
  geom_bar(stat = "identity")

gene_selection <- hvg_list[1:3000] %>% names()

# Create merged object:
integrated_data <- purrr::reduce(integrated_data,
                          merge,
                          merge.data = TRUE)

print("You managed to merge everything")
print("Object size")
print(object.size(integrated_data))

DefaultAssay(integrated_data) <- def_assay

# Process it before integration -----------------------
integrated_data <- integrated_data %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(features = gene_selection,
         npcs = 30,
         verbose = FALSE)


original_pca_plt <- DimPlot(object = integrated_data,
                            reduction = "pca",
                            pt.size = .1,
                            group.by = "orig.ident")


# Add batch info
# must contain orig.ident
# batch_info <- read_csv(batch_file)
# batch_covars <- colnames(batch_info)
#
# tmp_meta <- integrated_data@meta.data %>%
#   left_join(batch_info, by = "orig.ident")
#
# integrated_data@meta.data <- bind_cols(integrated_data@meta.data, tmp_meta[, batch_covars[-1]])

# Integrate the data -----------------------
integrated_data <- RunHarmony(integrated_data,
                              group.by.vars = c("orig.ident"),
                              plot_convergence = TRUE,
                              assay.use = def_assay,
                              max.iter.harmony = 20)

corrected_pca_plt <- DimPlot(object = integrated_data,
                             reduction = "harmony",
                             pt.size = .1,
                             group.by = "orig.ident")

# Create the UMAP with new reduction -----------
integrated_data <- integrated_data %>%
  RunUMAP(reduction = "harmony", dims = 1:30,
          reduction.name = "umap_harmony") %>%
  RunUMAP(reduction = "pca", dims = 1:30,
          reduction.name = "umap_original")

integrated_data <- FindNeighbors(integrated_data,
                                 reduction = "harmony",
                                 dims = 1:30)

DefaultAssay(integrated_data) <- def_assay

optimize= T
if(optimize) {

  # Clustering and optimization -------------------------
  print("Optimizing clustering")

  seq_res <- seq(0.4, 1.5, 0.1)

  integrated_data <- FindClusters(object = integrated_data,
                                  resolution = seq_res,
                                  verbose = F)

  clustree_plt <- clustree(integrated_data,
                           prefix = paste0(DefaultAssay(integrated_data), "_snn_res."))

  # Optimize clustering ------------------------------------------------------
  cell_dists <- dist(integrated_data@reductions$harmony@cell.embeddings,
                     method = "euclidean")


  cluster_info <- integrated_data@meta.data[,grepl(paste0(DefaultAssay(integrated_data),"_snn_res"),
                                                   colnames(integrated_data@meta.data))] %>%
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

  integrated_data[["opt_clust_integrated"]] <- integrated_data[[names(which.max(silhouette_res))]]

  Idents(integrated_data) = "opt_clust_integrated"

  # Reduce meta-data -------------------------------------------------------------------------
  spam_cols <- grepl(paste0(DefaultAssay(integrated_data), "_snn_res"),
                     colnames(integrated_data@meta.data)) |
    grepl("seurat_clusters",colnames(integrated_data@meta.data))

  integrated_data@meta.data <- integrated_data@meta.data[,!spam_cols]

} else {

  print("Not Optimizing clustering")

  seq_res <- seq(0.4, 1.6, 0.2)

  integrated_data <- FindClusters(integrated_data,
                                  resolution = seq_res,
                                  verbose = F)

  clustree_plt <- clustree(integrated_data,
                           prefix = paste0(DefaultAssay(integrated_data),
                                           "_snn_res."))

  integrated_data <- FindClusters(integrated_data,
                                  resolution = default_resolution,
                                  verbose = F)

  integrated_data[["opt_clust_integrated"]] <- integrated_data[["seurat_clusters"]]

  spam_cols <- grepl(paste0(DefaultAssay(integrated_data), "_snn_res"),
                     colnames(integrated_data@meta.data)) |
    grepl("seurat_clusters",colnames(integrated_data@meta.data))

  integrated_data@meta.data <- integrated_data@meta.data[,!spam_cols]

}

# add group label to meta data ----------------------------------------------------------------

x= integrated_data[[]]
x= x %>% mutate(group= ifelse(grepl("ct", orig.ident), "ct", "hf"))
integrated_data= AddMetaData(integrated_data, x$group, col.name = "group")

# Save object ------------------------------------------------------

saveRDS(integrated_data, file = "output/seu.objs/seu_harmony.rds")
integrated_data= readRDS("output/seu.objs/seu_harmony.rds")

# Print QC file ------------------------------------------------------

umap_corrected_sample <- DimPlot(object = integrated_data,
                                 reduction = "umap_harmony",
                                 pt.size = .1,
                                 group.by = "orig.ident",
                                 label=T)

umap_corrected_clustering <- DimPlot(object = integrated_data,
                                     reduction = "umap_harmony",
                                     pt.size = .1,
                                     group.by = "opt_clust_integrated",
                                     label = T)

umap_sample <- DimPlot(object = integrated_data,
                       reduction = "umap_original",
                       pt.size = .1,
                       group.by = "orig.ident")

umap_clustering <- DimPlot(object = integrated_data,
                           reduction = "umap_original",
                           pt.size = .1,
                           group.by = "opt_clust_integrated",
                           label = T)

x= integrated_data[[]]

precent.mt.p= ggplot(x, aes(x= opt_clust_integrated, y= percent.mt))+
  geom_point()+
  geom_boxplot()
percent.rb.p= ggplot(x, aes(x= opt_clust_integrated, y= percent.rb))+
  geom_point()+
  geom_boxplot()
precent.dis.p= ggplot(x, aes(x= opt_clust_integrated, y= dissociation_s1))+
  geom_point()+
  geom_boxplot()
s.score.p= ggplot(x, aes(x= opt_clust_integrated, y= S.Score))+
  geom_point()+
  geom_boxplot()
g2m.score.p= ggplot(x, aes(x= opt_clust_integrated, y= G2M.Score))+
  geom_point()+
  geom_boxplot()
count.p= ggplot(x, aes(x= opt_clust_integrated, y= nCount_RNA))+
  geom_point()+
  geom_boxplot()
nfeature.p= ggplot(x, aes(x= opt_clust_integrated, y= nFeature_RNA))+
  geom_point()+
  geom_boxplot()
doublet.p= ggplot(x, aes(x= opt_clust_integrated, y= doublet_score))+
  geom_point()+
  geom_boxplot()

sanity_p1= plot_grid(precent.mt.p,percent.rb.p, precent.dis.p,
                    ncol = 1, align = "v")

sanity_p2= plot_grid(g2m.score.p ,
                     s.score.p,
                     nfeature.p,
                     count.p, ncol = 1, align = "v")

pdf(file = "output/figures/samplewiseQC/integrated_overview.pdf", height = 10, width = 12)

print(gene_selection_plt)
print(original_pca_plt)
print(corrected_pca_plt)
print(umap_sample)
print(umap_corrected_sample)
print(clustree_plt)
print(umap_clustering)
print(umap_corrected_clustering)
print(sanity_p1)
print(sanity_p2)
dev.off()

# # Give reductions to ease future analysis
#
# reductions_list <-  list(meta_data = integrated_data@meta.data,
#                          reduction = integrated_data@reductions[["umap_harmony"]])
#
# saveRDS(reductions_list,
#         file = gsub("[.]rds", "_umap.rds", out_file))
