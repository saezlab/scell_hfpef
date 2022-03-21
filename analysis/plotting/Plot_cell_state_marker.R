## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-16
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## cell state marker plot



library(Seurat)
library(tidyverse)
library(WriteXLS)

source("analysis/utils.R")

# load sc data
seu.objs=  readRDS("output/seu.objs/cell.state.obj.list.rds")
seu= seu.objs$fibroblasts
#rm(seu.objs)

statemarker= readRDS("output/cell_state_marker.rds")




# FIBROBLAST_markerplot -----------------------------------------------------------------------


statemarker$fibroblasts%>% filter(cluster==4)%>% filter(p_val_adj<0.05)%>% arrange(desc(avg_log2FC))%>%
  mutate(base= ifelse(gene %in% NABA_SETS_mouse$Basement_membranes,1,0),
         col= ifelse(gene %in% NABA_SETS_mouse$Collagens,1,0),
         core= ifelse(gene %in% NABA_SETS_mouse$Core_matrisome,1,0),
         matri= ifelse(gene %in% NABA_SETS_mouse$Matrisome,1,0),
         prote= ifelse(gene %in% NABA_SETS_mouse$Proteoglycans,1,0)
  )


## get marker genes to plot :

fib.state.maker.laura= c("Col15a1", "Hsd11b1", "Cxcl14", "Igfbp3","Comp", "Timp3", "Mgp","Cd248", "Gfpt2", "Pi16", "Cxcl1",
"Isg15", "Ifit3", "Cilp", "Thbs4", "Ltbp2", "Wif1", "Dkk3")

fib.state.maker.me =c("Col15a1", "Pi16", "Comp", "Cxcl12", "Coch", "Ccl19")

state_0= statemarker$fibroblasts %>% filter(cluster==0)%>% filter(p_val_adj<0.05)%>% arrange(desc(avg_log2FC))%>% pull(gene)
#saveRDS(state_0[1:99], "output/fib_state_vasc?.rds")

marker_max_FC= statemarker$fibroblasts %>%
  group_by(cluster)%>%
  filter(p_val_adj<0.05)%>% top_n(5,(avg_log2FC)) %>%
  pull(gene)

seu@meta.data= seu@meta.data%>%
mutate(cellstate= factor(cellstate, levels = c("Col15a1+", "Igfbp3+",
                                              "Pi16+","Cxcl1+","Cilp+","Wif1+")))

p.marker= DotPlot(seu, group.by = "cellstate", features = c("Col4a1", "Col15a1", unique(marker_max_FC)), cols= c("#0353A4", "#C41E3D"),
                  cluster.idents = F)+
  coord_flip()+
  theme(axis.text.x= element_text(angle= 40, hjust = 1))+
  labs(x= "", y= "")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p.marker

pdf("output/figures/cell_type_assignment/fibroblast_cellstate_marker.pdf",
    width = 5,
    height= 6)
p.marker
dev.off()


p.marker= DotPlot(seu, group.by = "cellstate", features = unique(unlist(markx)), cols= c("#0353A4", "#C41E3D"),
                  cluster.idents = F)+
  coord_flip()+
  theme(axis.text.x= element_text(angle= 40, hjust = 1))+
  labs(x= "", y= "")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


# FIBROBLAST_UMAP----------------------------------------------------------------------------------------

cols= c("#F6EFA6","#EF946C", "#BEAED4", "#77B6EA","#ab87a6", "#66B3BA","#848FA5", "#58B09C", "#AF1B3F")

umap1= DimPlot(seu,
               group.by= "cellstate",
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
  theme(#axis.line = element_blank(),
    axis.text = element_blank(),
    axis.title= element_text(size= 9),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )
umap1


pdf("output/figures/UMAP_fibs.pdf",
    height= 3,
    width= 3)
umap1
dev.off()



# MACS_MARKER ---------------------------------------------------------------------------------

macs= seu.objs$macrophages

marker_max_FC= statemarker$macrophages %>%
  group_by(cluster)%>%
  filter(p_val_adj<0.05)%>% top_n(5,(avg_log2FC)) %>%
  pull(gene)
Laura_marker= c("Pglyrp1", "Itgal", "Ear2","Mrc1", "Cd163", "Cbr2","Ly6c2", "Chil3","Spp1", "Cxcl2","H2-Eb1", "H2-Ab1")

macs@meta.data= macs@meta.data%>%
  mutate(cellstate= factor(cellstate, levels = c("Col15a1+", "Igfbp3+",
                                                 "Pi16+","Cxcl1+","Cilp+","Wif1+")))

p.marker= DotPlot(macs, group.by = "cellstate",
                  #features = c("Ccr2", "MhcII", "Cxcl2", unique(marker_max_FC)),
                  features= Laura_marker,
                  cols= c("#0353A4", "#C41E3D"),
                  cluster.idents = F)+
  coord_flip()+
  theme(axis.text.x= element_text(angle= 40, hjust = 1))+
  labs(x= "", y= "")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p.marker

pdf("output/figures/cell_type_assignment/macrophages_cellstate_marker.pdf",
    width= 4
    height= 10  )
p.marker
dev.off()
# MACS_UMAP -----------------------------------------------------------------------------------

cols= c("#EF946C", "#785474", "#AF1B3F","#77B6EA","#66B3BA","#ab87a6", "#66B3BA","#848FA5", "#58B09C", "#AF1B3F")

umap1= DimPlot(macs,
               group.by= "cellstate",
               #cols = sample(col_vector, 50 ,replace = F),
               cols= cols,
               label = F,
               label.size = 4,
               label.col= "black",
)+
  #NoLegend()+
  ggtitle("")+
  labs(x= "UMAP1",
       y= "UMAP2")

umap1= umap1+
  theme(#axis.line = element_blank(),
    axis.text = element_blank(),
    axis.title= element_text(size= 9),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )
umap1


pdf("output/figures/UMAP_macs.pdf",
    height= 3,
    width= 6)
umap1
dev.off()

