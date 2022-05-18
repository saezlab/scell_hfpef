## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-09
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## plot marker genes for celltype assignment

library(Seurat)
library(tidyverse)
source("analysis/utils.R")

# process:
marker= read.csv("data/Celltypes_marker_plot.csv")
marker = str_replace_all(unlist(str_split(marker$Markers, ",")), " ", "")
marker= marker[!marker %in% c("", "Vcam1")]


seu= readRDS("output/seu.objs/integrated_cellstate_noECs.rds")

seu@meta.data= seu@meta.data %>%
  mutate(celltype = str_replace_all(celltype, "fibroblasts", "Fibroblasts"),
         celltype = str_replace_all(celltype, "Smooth.muscle_pericytes", "SMC/Pericytes"),
         celltype = str_replace_all(celltype, "granulocytes", "Granulocytes"),
         celltype = str_replace_all(celltype, "macrophages", "Macrophages"),
         celltype2 = str_replace_all(celltype2, "macrophages", "Macrophages"),
         celltype2 = str_replace_all(celltype2, "fibroblasts", "Fibroblasts"),
         celltype2 = str_replace_all(celltype2, "Smooth.muscle_pericytes", "SMC/Pericytes"),
         celltype2 = str_replace_all(celltype2, "granulocytes", "Granulocytes"),
         celltype2 = str_replace_all(celltype2, "endo", "Endo"),
         celltype2 = str_replace_all(celltype2, "Cd", "CD"),
         celltype2 = str_replace_all(celltype2, "T.cell", "T.cells")
         )

saveRDS(seu, "output/seu.objs/integrated_cellstate_nameupdated.rds")

seu= readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")

# PLOT cell marker features: ------------------------------------------------------------------
##update factor levels
seu@meta.data= seu@meta.data %>%
  mutate(celltype  = factor(celltype,levels = c("Fibroblasts",
                                                "Endothelial",
                                                "NK.cells",
                                                "Macrophages",
                                                "T.cells",
                                                "B.cells",
                                                "Granulocytes",
                                                "SMC/Pericytes")))%>%
  mutate(celltype2  = factor(celltype2,levels = c("Endothelial",
                                                  "Fibroblasts.1",
                                                "Fibroblasts.2",
                                                "NK.cells",
                                                "Macrophages",
                                                "CD8.T.cells",
                                                "CD4.T.cells",
                                                "B.cells",
                                                "Granulocytes",
                                                "SMC/Pericytes")))

Idents(seu)= "celltype2"
DefaultAssay(seu)= "RNA"

p.marker= DotPlot(seu, features = marker,
                  cols= c("#0353A4", "#C41E3D"),
                  cluster.idents = F)+
  coord_flip()+
  theme(axis.text.x= element_text(angle= 40, hjust = 1))+
  labs(x= "", y= "")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)
        )

p.marker


cols= c("#EF946C", "#BEAED4", "#77B6EA","#785474","#A4B494", "#848FA5", "#58B09C", "#66B3BA", "#AF1B3F","#F6EFA6")

umap1= DimPlot(seu,
               group.by= "celltype2",
               #cols = sample(col_vector, 50 ,replace = F),
               cols= cols,
               label = T,
               label.size = 4,
               label.col= "black",
)+
  NoLegend()+
  ggtitle("")+
  labs(x= "UMAP1",
       y= "UMAP2")

umap1= umap1+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size= 10),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
)

umap1



# Compose for figure 1 ------------------------------------------------------------------------
p.distance_ratio= readRDS("output/figures/main/p.distance_ratio_obj.rds")
p.mean.= readRDS( "output/figures/main/proportions_obj.rds")
p.mean.= p.mean.+theme( panel.border = element_rect(colour = "black", fill=NA, size=1))
p1= cowplot::plot_grid(umap1+
                     theme(line = element_blank()),
                   #p.marker,
                   p.mean.,
                   #p.distance_ratio,
                   ncol =1,
                   rel_widths = c(1,1.5),
                   align = "hv",
                   axis= "lr",
                   labels= c("C", "D", "E","F"))
?plot_grid

p2= cowplot::plot_grid(p.marker,
                   p.distance_ratio,
                   ncol =1,
                   rel_widths = c(1,1.5),
                   align = "hv",
                   axis= "lr",
                   rel_heights = c(1,0.6),
                   labels= c("C", "D", "E","F"))

cowplot::plot_grid(umap1,
                   p.marker,

                   ncol =2,
                   rel_widths = c(1,1.5),
                   align = "h",
                   axis= "t",
                   rel_heights = c(1,0.6),
                   labels= c("C", "D", "E","F"))
library(cowplot)

p.mean. +theme(axis.text = element_text(size= 12, color = "black"))

p.mean. = unify_axis(p.mean.)
p.distance_ratio= unify_axis(p.distance_ratio)
p.marker= unify_axis(p.marker)

left= align_plots(umap1, p.mean., axis= "l", align ="v")
right= align_plots(p.marker,p.distance_ratio, axis= "r", align ="v")
top.row= plot_grid(left[[1]], right[[1]], align= "h", axis = "t", ncol = 2, rel_widths = c(1,1.4))
bottom.row= plot_grid(left[[2]], right[[2]], align= "h", axis = "t", ncol = 2)

p_comb= cowplot::plot_grid(top.row,   bottom.row,

                   rel_heights = c(1.6, 1),
                   ncol =1)
p_comb
#
# plot_grid(plotlist= left, ncol= 1)
#
# plot_grid(umap1, p.mean., axis= "v", align ="l", ncol= 1)
# plot_grid(p.marker, p.distance_ratio, axis= "l", align ="l", ncol= 1)

# save ----------------------------------------------------------------------------------------


pdf("output/figures/UMAPs_overview2.pdf",
    width= 4,
    height=4)
umap1
dev.off()

pdf("output/figures/Celltype_marker.pdf",
    height= 5)
p.marker
dev.off()

pdf("output/figures/main/figure_1_panel.pdf",
    height= 8,
    width= 10)
p_comb
dev.off()
