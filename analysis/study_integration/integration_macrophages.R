## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-10-19
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Integration of macrophages between different studies.

library(tidyverse)
library(Seurat)
library(harmony)
library(cluster)
library(clustree)
library(cowplot)
library(purrr)

# # Load SC objects and unifydate meta data: --------------------------------------------------
# Ang II model from circulation:

circ= readRDS( file = "output/seu.objs/study_integrations/macs/circ_macs.rds")
forte= readRDS( file = "output/seu.objs/study_integrations/macs/mi_macs.rds")
hfpef= readRDS("output/cell_specific_obj/macrophages.rds")

hfpef@meta.data= hfpef@meta.data %>%
  mutate(group = ifelse(group== "hfpef", "hf", "ct"),
         study= "hfpef")

# Integration via Harmony ---------------------------------------------------------------------

## 1. Rename cells to facilitate back mapping later:
rename_cells= function(seu, string_add_to_cellid){
  seu = RenameCells(seu,  string_add_to_cellid)
}

hfpef= rename_cells(hfpef, "a")
forte= rename_cells(forte, "b")
circ= rename_cells(circ, "c")

DefaultAssay(hfpef)= "RNA"
DefaultAssay(forte) = "RNA"
DefaultAssay(circ)= "RNA"

DimPlot(circ)
DimPlot(hfpef)
DimPlot(forte)


# Calculate HVGs  ------------------------------------------------
## Calculate HVGs for each sample in each study, and select those that appear most often

seu.obj= list(circ,
              hfpef,
              forte)
rm(circ)
rm(hfpef)
rm(forte)

calculate_hvgs= function(seu, sample){

  integrated_data <- SplitObject(seu, split.by = sample)

  # First get all the most variable genes per patient per batch

  hvg_list<- map(integrated_data, function(x) {

    DefaultAssay(x) <- "RNA"

    x <- FindVariableFeatures(x, selection.method = "vst",
                              nfeatures = 3000)

    x@assays[["RNA"]]@var.features

  }) %>% unlist()

  hvg_list <- table(hvg_list) %>%
    sort(decreasing = TRUE)

  hist(hvg_list)

  gene_selection= hvg_list[which(hvg_list==max(hvg_list))]
  print(paste0("hvg genes with max number=", length(gene_selection)))
  i= 1

  while(length(gene_selection)<2000){
    gene_selection= c(gene_selection,hvg_list[which(hvg_list==(max(hvg_list)-i))])
    i=i+1
    print(paste0("hvg genes incresed to ", length(gene_selection)))
  }

  return(names(gene_selection))

}

x= map(seu.obj, function(x){calculate_hvgs(x, "orig.ident")})

# These are now the hvg genes per study:
hvg_list <- table(unlist(x)) %>%
  sort(decreasing = TRUE)

hist(hvg_list)
# we will select only those hvg genes that ware hvgs in each study: for 3 studies-> 3
gene_selection= hvg_list[which(hvg_list>1)]
gene_selection= names(gene_selection)

# Read object to subset -----------------------------------------------------

#   # Genes that are stable among batches
# gene_selection <- hvg_list[hvg_list == 2] %>% names()
#
# Re-integrate data

integrated_data <- reduce(seu.obj,
                          merge,
                          merge.data = TRUE)

rm(seu.obj)

#integrated_data

DefaultAssay(integrated_data) = "RNA"

integrated_data@meta.data <- integrated_data@meta.data %>% select(orig.ident, nCount_RNA, nFeature_RNA, percent.mt, group, study)
colnames(integrated_data@meta.data)

## Scale -> PCA -> Harmony
integrated_data <- integrated_data %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(features = gene_selection,
          npcs = 30,
         verbose = FALSE)

# Original PCs
original_pca_plt <- DimPlot(object = integrated_data,
                            reduction = "pca",
                            pt.size = 1,
                            group.by = "orig.ident")

original_study <- DimPlot(object = integrated_data,
                            reduction = "pca",
                            pt.size = 0.4,
                            group.by = "study")

#saveRDS(integrated_data, "output/seu.objs/study_integrations/harmony_fib_pre.rds")

integrated_data@meta.data %>% group_by(study) %>% count
integrated_data@meta.data %>% group_by(group, study) %>% count

# Integrate the data -----------------------

integrated_data <- RunHarmony(integrated_data,
                              group.by.vars= c("study", "orig.ident"),
                              plot_convergence = TRUE,
                              max.iter.harmony = 20,
                              project.dim = F,
                              assay.use= "RNA")

# saveRDS(integrated_data, "output/seu.objs/study_integrations/harmony_macs.rds")
# integrated_data = readRDS( "output/seu.objs/study_integrations/harmony_macs.rds")

corrected_pca_plt <- DimPlot(object = integrated_data,
                             reduction = "harmony",
                             pt.size =0.4,
                             group.by = "study")

DefaultAssay(integrated_data)

# RUN UMAP and clustering:

integrated_data <- integrated_data %>%
  RunUMAP(reduction = "harmony", dims = 1:30,
          reduction.name= "umap_harmony_study")%>%
  RunUMAP(reduction = "pca", dims = 1:30,
          reduction.name = "umap_original")


integrated_data <- FindNeighbors(integrated_data,
                                 reduction = "harmony",
                                 dims = 1:30)



# Create new clustering
print("Optimizing clustering")

seq_res <- seq(0.1, 1, 0.1)

integrated_data <- FindClusters(object = integrated_data,
                                resolution = seq_res,
                                verbose = F)

#DimPlot(integrated_data, group.by = "RNA_snn_res.0.2")

clustree_plt <- clustree(integrated_data,
                         prefix = paste0(DefaultAssay(integrated_data), "_snn_res."))

cell_dists <- dist(integrated_data@reductions$harmony@cell.embeddings,
                   method = "euclidean")


cluster_info <- integrated_data@meta.data[,grepl(paste0("RNA_snn_res"),
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

rm(cell_dists)

integrated_data[["opt_clust_integrated"]] <- integrated_data[[names(which(silhouette_res == max(silhouette_res)))]]
#integrated_data[["opt_clust_integrated"]] = integrated_data[["RNA_snn_res.0.5"]]
Idents(integrated_data) = "opt_clust_integrated"

# Reduce meta-data -------------------------------------------------------------------------
spam_cols <- grepl(paste0(DefaultAssay(integrated_data), "_snn_res"),
                   colnames(integrated_data@meta.data))
integrated_data@meta.data <- integrated_data@meta.data[,!spam_cols]


saveRDS(integrated_data, "output/seu.objs/study_integrations/macs/harmony_macs.rds")
integrated_data = readRDS( "output/seu.objs/study_integrations/macs/harmony_macs.rds")

p.harmony.pca= DimPlot(integrated_data, reduction = "harmony")
p.harmony.umap= DimPlot(integrated_data, reduction = "umap_harmony_study", group.by = "opt_clust_integrated", label = T)
p.harmony.study= DimPlot(integrated_data, reduction = "umap_harmony_study", group.by = "study")
p.harmony.sample= DimPlot(integrated_data, reduction = "umap_harmony_study", group.by = "orig.ident")

source("code/utils.R")
seu.meta= integrated_data[[]]

fibs.counts= seu.meta%>% group_by(study) %>% count()
x= seu.meta%>% group_by(opt_clust_integrated, study) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
x=x %>%
  left_join(fibs.counts %>%
              rename(n.all= n),
            by= "study") %>%
  mutate(freq2= n/n.all)

p.int.fib.prop= ggplot(x, aes(x= study, y= freq2, fill =opt_clust_integrated))+
  geom_col()+
  scale_fill_manual(values = col_vector)+
  theme_minimal()

x2= x %>% mutate(freak= sum(freq2)) %>% ungroup() %>% mutate(sacaled.prop= freq2/freak)
x2 %>% print(n=100)
p.proportion.sclaed= ggplot(x2, aes(x= opt_clust_integrated, y= sacaled.prop, fill =study))+
  geom_col()+
  geom_hline(yintercept = c(1/3, 2/3))+
  theme_minimal()
FeaturePlot(integrated_data, features = c("Lyz2", "Cd68", "Cx3cr1", "Postn"))

p.marker = VlnPlot(integrated_data, c("Postn", "Meox1", "Cilp", "Igfbp3", "Dcn", "Col1a1","Fap"))
p.study=  calc_props(seu.meta, cluster.col = "opt_clust_integrated", "study")
p.sample=  calc_props(seu.meta, cluster.col = "opt_clust_integrated", "orig.ident")

pdf("output/figures/mac_integration_report.pdf")
original_pca_plt
original_study
p.harmony.pca
p.harmony.sample
p.harmony.study
#clustree_plt
p.harmony.umap
p.study
p.sample
p.int.fib.prop
p.proportion.sclaed
p.marker
dev.off()

