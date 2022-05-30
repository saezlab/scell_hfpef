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
## after identifying problematic cluster downstream in cluster_label.R,
## we are going to delete those, reintegrate the data , recluster.
## Also calculation of proportion is done

library(tidyverse)
library(Seurat)
library(harmony)
library(cluster)
library(clustree)
library(cowplot)
library(purrr)

integrated_data= readRDS("output/seu.objs/integrated_cellstate.rds")
def_assay= "RNA"

unique(integrated_data$cellstate)

# remove three cell clusters that appear to be low Q:
integrated_data= subset(integrated_data,
                        cellstate %in% c("mesenchymal.activated.Ec"),
                        invert= T)

# from downstream analysis we will add more problematic cells for the removal
source("analysis/utils.R")
seu.meta= integrated_data[[]]
# qc_plots= get_QC_plots(seu.meta, "celltype2")
#
# pdf("output/figures/QCs/celltype2.QC.plots.pdf")
# print(qc_plots)
# dev.off()

#-> Cluster in Question was Endothelial2 but appears fine

# NEW INT -------------------------------------------------------------------------------------

#split by sample=
integrated_data <- SplitObject(integrated_data, split.by = "orig.ident")

# First get all the most variable genes per patient

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
         pcs= 30,
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
                              max.iter.harmony = 20,
                              reduction.save = "harmony_filt")


# Corrected dimensions -----------------------
corrected_pca_plt <- DimPlot(object = integrated_data,
                             reduction = "harmony_filt",
                             pt.size =1,
                             group.by = "orig.ident")

# Create the UMAP with new reduction -----------
integrated_data <- integrated_data %>%
  RunUMAP(reduction = "harmony_filt", dims = 1:30,
          reduction.name = "umap_harmony_filt") %>%
  RunUMAP(reduction = "pca", dims = 1:30,
          reduction.name = "umap_original_filt")

#compare both UMAPs
DimPlot(integrated_data, reduction= "umap_original_filt", label= T, group.by = "opt_clust_integrated")
DimPlot(integrated_data, reduction= "umap_harmony_filt", label= T, group.by = "opt_clust_integrated")

# we need to remove one celltype2 labeling
integrated_data@meta.data =integrated_data@meta.data %>% mutate(celltype2= ifelse(celltype2== "mesenchymal.endothelial", "endothelial", celltype2))
unique(integrated_data$celltype2)
#integrated_data= FindVariableFeatures(integrated_data)
saveRDS(integrated_data, "output/seu.objs/integrated_cellstate_noECs.rds")
#integrated_data= readRDS( "output/seu.objs/integrated_cellstate_noECs.rds")

##plot new UMAPS:
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p1= DimPlot(integrated_data, reduction= "umap_harmony_filt", label= T, group.by = "opt_clust_integrated", cols = col_vector)
p2= DimPlot(integrated_data, reduction= "umap_harmony_filt", label= T, group.by = "celltype", cols = col_vector)
p3= DimPlot(integrated_data, reduction= "umap_harmony_filt", label= T, group.by = "celltype2", cols = col_vector)
p4= DimPlot(integrated_data, reduction= "umap_harmony_filt", label= T, group.by = "cellstate", cols = col_vector)


pdf("output/figures/UMAPs_overview_noEC.pdf")
p1
p3
p2
p4
dev.off()
