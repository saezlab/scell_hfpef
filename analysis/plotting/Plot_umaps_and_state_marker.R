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


library(ggplot2)
library(Seurat)
library(tidyverse)
library(WriteXLS)

source("code/utils.R")

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
  #coord_flip()+
  theme(axis.text.x= element_text(angle= 90, hjust = 1))+
  labs(x= "", y= "")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p.marker

pdf("output/figures/main/fibroblast_cellstate_marker.pdf",
    width = 8.5,
    height= 2.5)
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
  theme(axis.line = element_blank(),
    axis.text = element_blank(),
    axis.title= element_text(size= 9),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
umap1
saveRDS(umap1, "output/figures/main/umap_fibs.rds")

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
plot_marker= c("Lyve1", "Mrc1", "C1qc", "Adgre1",
               "Itgal", "Ear2", "Eno3", "Ace",
               "Cxcl2", "Spp1", "Gpnmb", "Ctsd","Ccl2",
               "H2-Ab1", "H2-Eb1", "H2-Aa","Cd74",
               "Ccr2",
                "Ly6c2", "Fn1", "Plac8"
               )
p.marker= DotPlot(macs, group.by = "cellstate",
                  #features = c("Ccr2", "MhcII", "Cxcl2", unique(marker_max_FC)),
                  #features= Laura_marker,
                  features= plot_marker,
                  cols= c("#0353A4", "#C41E3D"),
                  cluster.idents = F)+
  coord_flip()+
  theme(axis.text.x= element_text(angle= 40, hjust = 1))+
  labs(x= "", y= "")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p.marker

pdf("output/figures/main/macrophages_cellstate_marker.pdf",
    width= 4.5,
    height= 6 )
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


# UMAP of integrated fibroblasts --------------------------------------------------------------

integrated_data= readRDS( "output/seu.objs/study_integrations/harmony_fib_filt.rds")
integrated_data@meta.data=  integrated_data@meta.data%>%
  mutate(study= ifelse(study=="circ", "AngII",
                       ifelse(study=="forte", "MI", "HFpEF")),
         group = ifelse(group== "ct", "Control", "HF"))
source("code/utils.R")
#cols= c( "#7FC97F" ,"#BEAED4", "#FDC086" )
#cols= c("#EF946C", "#785474", "#AF1B3F","#77B6EA","#66B3BA","#ab87a6", "#66B3BA","#848FA5", "#58B09C", "#AF1B3F")
umap1= DimPlot(integrated_data,
               group.by= "opt_clust_integrated",
               #cols = sample(col_vector, 50 ,replace = F),
               cols= col_vector,
               label = F,
               label.size = 4,
               label.col= "black",
)+
  #NoLegend()+
  ggtitle("")+
  labs(x= "UMAP1",
       y= "UMAP2")+
  theme(axis.line = element_blank(),
    legend.text=element_text(size=11),
    axis.title= element_text(size= 9),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))


umap1
#cols= c("#EF946C", "#785474", "#AF1B3F","#77B6EA","#66B3BA","#ab87a6", "#66B3BA","#848FA5", "#58B09C", "#AF1B3F")
umap.disease= DimPlot(integrated_data,
               group.by= "group",
               pt.size = 0.01,
               #cols = sample(col_vector, 50 ,replace = F),
               cols= hfpef_cols,
               label = F,
               label.size = 4,
               label.col= "black",
)+
  #NoLegend()+
  ggtitle("")+
  labs(x= "UMAP1",
       y= "UMAP2")+
  theme(axis.line = element_blank(),
    axis.text = element_blank(),
        legend.text=element_text(size=11),
    axis.title= element_text(size= 9),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.margin = unit(c(-1, -1, -1, -1), "cm")
  )
umap.disease

umap.study= DimPlot(integrated_data,
               group.by= "study",
               #cols = sample(col_vector, 50 ,replace = F),
               cols= c(  "#EE7674",
                         "#494304",
                         col_vector[5],
                         col_vector[8]),
               label = F,
               label.size = 4,
               pt.size = 0.1,
               label.col= "black",
               order = rev(c("MI", "AngII", "HFpEF"))
)+
  #NoLegend()+
  ggtitle("")+
  labs(x= "UMAP1",
       y= "UMAP2")+
  theme(axis.line = element_blank(),
    axis.text = element_blank(),

    legend.text=element_text(size=11),
    axis.title= element_text(size= 9),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.margin = unit(c(-1, -1, -1, -1), "cm")
  )
umap.study
#


# plot umap for study sep
#load time data separately since it is not stored in integrated object
x= readRDS("output/seu.objs/study_integrations/meta_with_time.rds")
x =  x %>%
  mutate(time= ifelse(is.na(time), 7, time))%>%
  mutate(study= ifelse((time == 14 & study=="forte") | (time== 7 & group == "ct"),
                       "MI_late", study),
         study= ifelse(study == "forte" & time != 14,
                       "MI_early",
                       study))%>%
  mutate(study= ifelse(study=="circ", "AngII",
                       ifelse(study=="hfpef","HFpEF", study)),
         group = ifelse(group== "ct", "Control", "HF"))%>%
  rename(study2= study)

table(x$time, x$study2)
meta= integrated_data@meta.data
meta= meta %>% as.data.frame() %>% rownames_to_column("cellid") %>% left_join(x%>%
                                                                          select(cellid, time, study2), by= "cellid")
integrated_data= AddMetaData(integrated_data, meta$study2, "study2")


umap.study2= DimPlot(integrated_data,
                     group.by= "study2",
                     #cols = sample(col_vector, 50 ,replace = F),
                     cols= c("#F9F19A", "#FABF87", "#82C180",  "#BEAED4"),
                     label = F,
                     label.size = 4,
                     pt.size = 0.1,
                     label.col= "black",
                     order = rev(c("MI_early", "MI_late", "AngII", "HFpEF"))
)+
  #NoLegend()+
  ggtitle("")+
  labs(x= "UMAP1",
       y= "UMAP2")+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),

        legend.text=element_text(size=11),
        axis.title= element_text(size= 9),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.margin = unit(c(-1, -1, -1, -1), "cm")
  )


umap.study2
# saveRDS(list(umap.disease, umap.study, umap1), "output/figures/main/Fig3/")
p.c= cowplot::plot_grid(
                        umap.disease,#+coord_equal(),
                        umap.study,#+coord_equal(),

                        umap1,#+coord_equal(),
                        rel_heights = c(1,  1),
                        ncol =3,
                        align="hv",
                        axis= "tplr")
p.c
pdf("output/figures/main/Fig3/fib_integrated_umaps.pdf",
    height= 3, width= 12)
p.c
dev.off()


pdf("output/figures/main/fib_integrated_umaps.pdf",
    height= 6, width= 6)
umap1+coord_equal()
umap.disease+coord_equal()
umap.study+coord_equal()
dev.off()

FeaturePlot(integrated_data,
                   features = c("Lyz2", "Spp1","Col1a1"
                   ),
                   coord.fixed = T,
                   combine = T)



umap3= FeaturePlot(integrated_data,
                   features = c("Cilp", "Postn","Acta2", "Actb", "Col4a1",
                                "Col15a1","Cxcl14", "Smoc2", "Dpep1",
                                "Pi16", "Igfbp3", "Wif1", "Ifit1", "Ccl2"
                                ),
                   coord.fixed = T,
                   combine = F

)

umap3= map(umap3, function(x){
  x+
  #ggtitle("")+
  labs(x= "UMAP1",
       y= "UMAP2")+
  theme(#axis.line = element_blank(),
    axis.text = element_blank(),
    axis.title= element_text(size= 11),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )
})

umap3.c= cowplot::plot_grid(plotlist = umap3, ncol = 3)

pdf("output/figures/main/UMAP_feature_int_fibs.pdf",
    height= 12,
    width= 10)
umap3.c
dev.off()
