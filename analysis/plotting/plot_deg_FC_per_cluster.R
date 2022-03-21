## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-26
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
##  plot deg expression of in fibs




library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

source("analysis/utils.R")


# # load data ---------------------------------------------------------------------------------


seu.objs=  readRDS("output/seu.objs/cell.state.obj.list.rds")

seu= seu.obj$fibroblasts
seu@meta.data= seu@meta.data%>% mutate(group_state= paste0(group,"_", cellstate))


pb.studies= readRDS( "output/fib_integration/pseudobulked.per.study.sample.rds")
#signature = readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
signature= readRDS("output/DEG_downsampled_intersect_hfpef.rds")


# plot FC -------------------------------------------------------------------------------------
cells_states= unique(seu@meta.data$cellstate)
Idents(seu)

fc.per.state= map(cells_states, function(x){
  seu2= subset(seu, cellstate==x)
  Idents(seu2)= "group"
  fc= Seurat::FoldChange(seu2, ident.1= "hf", ident.2= "ct")%>%
    rownames_to_column("gene")%>%
    mutate(cellstate= x)

}) %>% do.call(rbind, .)



Hmaps= map(signature$fibroblasts, function(x){
  df= fc.per.state %>% filter(gene %in% x)%>%
    select(gene, cellstate, avg_log2FC)%>%
    pivot_wider(names_from = cellstate, values_from = avg_log2FC)
  df= column_to_rownames(df, "gene")
  col_fun = colorRamp2(c(-1.5, 0, 1.5), c("darkblue", "white", "darkred"))
  col_fun(seq(-3, 3))
  Heatmap(df, col = col_fun, name= "avg_logFC")

})

pdf("output/figures/cell_type_analysis/fibros/fibroblast_deg_heatmap_up.pdf",
    width= 3,
    height= 12)
print(Hmaps[[1]])
dev.off()

pdf("output/figures/cell_type_analysis/fibros/fibroblast_deg_heatmap_dn.pdf",
    width= 3,
    height= 7)
print(Hmaps[[2]])
dev.off()


Idents(seu)= "group"
fc= Seurat::FindMarkers(seu, ident.1= "hf", ident.2= "ct")%>%
    rownames_to_column("gene")#%>%
    mutate(cellstate= x)

df= fc %>% filter(gene %in% unlist(signature$fibroblasts))%>%
  mutate(dummy= 1)%>%
select(gene,dummy , avg_log2FC)%>%
  pivot_wider(names_from = dummy, values_from = avg_log2FC)
df= column_to_rownames(df, "gene")
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("darkblue", "white", "darkred"))
col_fun(seq(-3, 3))
Heatmap(df, col = col_fun, name= "avg_logFC")

# plot_single_genes ---------------------------------------------------------------------------
genes_oi= c("Angptl4", "Slit3", "Sparc", "Adamtsl2", "Loxl1", "Mt1", "Mt2", "Rbm3")

df= map(genes_oi, function(x){
  p= FeaturePlot(seu, features = x, split.by = "group")
})

