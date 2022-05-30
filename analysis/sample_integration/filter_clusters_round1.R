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

integrated_data= readRDS("output/seu.objs/integrated_annotated.rds")
def_assay= "RNA"

unique(integrated_data$celltype)

# remove three cell clusters that appear to be low Q:
integrated_data= subset(integrated_data,
                        celltype %in% c("fibro+mitochondrial?", "lowQ.cells?", "endothelial/fibro?"),
                        invert= T)

# from downstream analysis we will add more problematic cells for the removal
source("analysis/utils.R")
seu.meta= integrated_data[[]]
qc_plots= get_QC_plots(seu.meta, "celltype2")

pdf("output/figures/QCs/celltype2.QC.plots.pdf")
print(qc_plots)
dev.off()

#-> Cluster in Question was Endothelial2 but appears fine


#Cluster0 in macrophage cell state is in question:
macs= readRDS("output/cell_specific_obj/macrophages.rds")
DimPlot(macs)
macs.met= macs[[]]
qc_plots_macs= get_QC_plots(macs.met, "opt_state")

#-> cluster 0 in macrophages subcluster has high mt, low count an low features, for that we will remove those cells.

pdf("output/figures/QCs/mac_lowQ_detection.pdf")
print(qc_plots_macs)
dev.off()

cells_to_remove= rownames(macs.met[macs.met$opt_state== 0, ])

integrated_data= subset(integrated_data, cells = cells_to_remove, invert = T)

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
                             reduction = "harmony",
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

#integrated_data= FindVariableFeatures(integrated_data)
saveRDS(integrated_data, "output/seu.objs/integrated_annotated_filtered.rds")
integrated_data= readRDS( "output/seu.objs/integrated_annotated_filtered.rds")

##plot new UMAPS:
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p1= DimPlot(integrated_data, reduction= "umap_harmony_filt", label= T, group.by = "opt_clust_integrated", cols = col_vector)
p2= DimPlot(integrated_data, reduction= "umap_harmony_filt", label= T, group.by = "celltype", cols = col_vector)
p3= DimPlot(integrated_data, reduction= "umap_harmony_filt", label= T, group.by = "celltype2", cols = col_vector)

pdf("output/figures/UMAPs_overview.pdf")
p1
p3
p2
dev.off()
