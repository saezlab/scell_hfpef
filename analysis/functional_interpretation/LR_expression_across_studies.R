## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2023-01-26
##
## Copyright (c) Jan D. Lanzer, 2023
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## compare t-values from hfpef found mac-fibro interactions between studies

library(tidyverse)
library(magrittr)
library(Seurat)

#load LR pairs in HFpEF
lr= readRDS(file= "output/liana_candidates_MF.rds")
lr$MF


get_logFC= function(seu,
                    feats= unlist(lr$MF)){
  names(feats)= NULL

  seu@meta.data= seu@meta.data %>% mutate(c.g= paste0(celltype,"_", group))
  Idents(seu)= "c.g"
  fibs= FindMarkers(seu,

                    assay = "RNA",
                    slot= "data",
                    ident.1 = "fibroblasts_hf",
                    ident.2= "fibroblasts_ct",
                    logfc.threshold = 0,
                    min.pct = 0.1
  )

  macs= FindMarkers(seu,
                    assay = "RNA",
                    slot= "data",
                    ident.1 = "macrophages_hf",
                    ident.2= "macrophages_ct",
                    logfc.threshold = 0,
                    min.pct = 0.1


  )


  Ma= macs[,] %>% as.data.frame()%>%
    rownames_to_column("gene")%>%
    mutate(celltype = "macrophage")

  Fi= fibs[,] %>% as.data.frame()%>%
    rownames_to_column("gene")%>%
    mutate(celltype = "fibroblasts")

  rbind(Ma, Fi)


}

#Ang II  ---------------------------------------------------------------

directory= "/home/jan/R-projects/sc-exploration/"
seu= readRDS( paste0(directory, "data/circ_MI/integrated_seurat_circMI2.rds"))
seu_meta= seu[[]]

mac.maker=  c("Lyz2", "Cd68", "Cx3cr1")
Idents(seu)= "integrated_snn_res.0.1"
p1= DimPlot(seu, group.by = "integrated_snn_res.0.1",label = T)
p2=DotPlot(seu, features = mac.maker)
p3= FeaturePlot(seu,
                features =  mac.maker
)

seu@meta.data= seu@meta.data %>%
  mutate(celltype= ifelse(integrated_snn_res.0.1 %in% c(0,8), "fibroblasts", ""),
         celltype= ifelse(integrated_snn_res.0.1 %in% c(1), "macrophages", celltype))%>%
  mutate(group = ifelse(treatment== "angiotensin II", "hf", "ct"),
         study= "angii")

Idents(seu)= "celltype"
DefaultAssay(seu)= "RNA"
seu2= subset(seu, celltype %in% c("macrophages", "fibroblasts"))
rm(seu)
angii=  get_logFC(seu2)
rm(seu2)
saveRDS(angii, "output/lgfc_angii.rds")
angii=readRDS("output/lgfc_angii.rds")

# MI ------------------------------------------------------------------------------------------

int_cell_MI = readRDS(file= paste0(directory, "output/cell_mi_files/integrateint_cell_MI.rds"))

## add meta data about
mac.maker=  c("Lyz2", "Cd68", "Cx3cr1")

marker= read.csv("data/Celltypes_marker_plot.csv")
marker = str_replace_all(unlist(str_split(marker$Markers, ",")), " ", "")

p1= DimPlot(int_cell_MI, group.by = "integrated_snn_res.0.1",label = T)
p2=DotPlot(int_cell_MI, features = marker)+
  coord_flip()
p3= FeaturePlot(int_cell_MI,
                features =  mac.maker
)

# add celltype group
int_cell_MI@meta.data= int_cell_MI@meta.data %>%
  mutate(celltype= ifelse(integrated_snn_res.0.1 %in% c(0,5,6), "fibroblasts", ""),
         celltype= ifelse(integrated_snn_res.0.1 %in% c(1, 8), "macrophages", celltype))%>%
  mutate(group = ifelse(OP== "myocardial infarction", "hf", "ct"),
         study= "forte")

Idents(int_cell_MI)= "celltype"
DefaultAssay(int_cell_MI)= "RNA"
int_cell_MI= subset(int_cell_MI,  celltype %in% c("macrophages", "fibroblasts"))

genes_mi= rownames(int_cell_MI)

#now create two subset, (early and late MI)
forte_early= subset(x = int_cell_MI, subset = time != 14 )
forte_late= subset(int_cell_MI, time %in% c(0,14) | group == "ct" )
rm(int_cell_MI)

#
# get_logFC= function(seu,
#                     feats= unlist(lr$MF)){
#   names(feats)= NULL
#
#   seu@meta.data= seu@meta.data %>% mutate(c.g= paste0(celltype,"_", group))
#   Idents(seu)= "c.g"
#   fibs= FindMarkers(seu,
#
#               assay = "RNA",
#               slot= "data",
#               ident.1 = "fibroblasts_hf",
#               ident.2= "fibroblasts_ct",
#               logfc.threshold = 0,
#               min.pct = 0
#               )
#
#   macs= FindMarkers(seu,
#                     assay = "RNA",
#                     slot= "data",
#                     ident.1 = "macrophages_hf",
#                     ident.2= "macrophages_ct",
#                     logfc.threshold = 0,
#                     min.pct = 0
#
#
#                       )
#
#
#   Ma= macs[feats,] %>% as.data.frame()%>%
#     rownames_to_column("gene")%>%
#     mutate(celltype = "macrophage")
#
#   Fi= fibs[feats,] %>% as.data.frame()%>%
#     rownames_to_column("gene")%>%
#     mutate(celltype = "fibroblasts")
#
#   rbind(Ma, Fi)
#
#
# }

MI_late=  get_logFC(forte_late)
saveRDS(MI_late, "output/lgfc_MI_late.rds")

rm(forte_late)


MI_early=  get_logFC(forte_early)
saveRDS(MI_early, "output/lgfc_MI_early.rds")
rm(forte_early)

MI_early= readRDS("output/lgfc_MI_late.rds")
MI_late= readRDS("output/lgfc_MI_late.rds")

# hfpef ---------------------------------------------------------------------------------------
seu= readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")
seu= subset(seu,  celltype %in% c("Macrophages", "Fibroblasts"))


seu@meta.data= seu@meta.data%>%
  mutate(celltype=tolower(celltype) ,
         group= ifelse(group== "hfpef", "hf", group))
HFpef_fc= get_logFC(seu)

saveRDS(HFpef_fc, "output/lgdc_hfpef.rds")

HFpef_fc= readRDS( "output/lgdc_hfpef.rds")
angii=readRDS("output/lgfc_angii.rds")
MI_late=readRDS("output/lgfc_MI_late.rds")
MI_early=readRDS("output/lgfc_MI_early.rds")

df= list("AngII"= angii,
     "MI_early"= MI_early,
     "MI_late"= MI_late,
     "HFpEF"= HFpef_fc)

df= map(names(df), function(x) (df[[x]]%>% mutate(study=x )))%>% do.call(rbind, .)

df %>%
  filter(gene %in% unlist(lr$MF))%>%
  ggplot(., aes(x= study, y=gene, fill = avg_log2FC))+
  facet_grid(rows= vars(celltype))+
  geom_tile()


library(ComplexHeatmap)

map.m= df %>%
  filter(gene %in% unlist(lr$MF$l))%>%
  filter(celltype== "macrophage")%>%
  select(gene,study, avg_log2FC)%>%
  pivot_wider(names_from=study , values_from = avg_log2FC)%>%
  select(gene, HFpEF, AngII, MI_early, MI_late)%>%
  as.data.frame()%>%
  column_to_rownames("gene")
map.m_p= df %>%
  filter(gene %in% unlist(lr$MF$l))%>%
  filter(celltype== "macrophage")%>%
  select(gene,study, p_val_adj)%>%
  pivot_wider(names_from=study , values_from = p_val_adj)%>%
  select(gene, HFpEF, AngII, MI_early, MI_late)%>%
  as.data.frame()%>%
  column_to_rownames("gene")

map.m%>%
  Heatmap(cluster_columns = F,
          name= "logFC", border = "black",rect_gp = gpar(col = "darkgrey"),  )

map.f= df %>%
  filter(gene %in% unlist(lr$MF$r))%>%
  filter(celltype== "fibroblasts")%>%
  select(gene,study, avg_log2FC)%>%
  pivot_wider(names_from=study , values_from = avg_log2FC)%>%
  select(gene, HFpEF, AngII, MI_early, MI_late)%>%
  as.data.frame()%>%
  column_to_rownames("gene")%>%
  Heatmap(cluster_columns = F,
          name= "logFC", border = "black",rect_gp = gpar(col = "darkgrey") )

pdf("output/figures/supp/MAC_ligand_comp.pdf",
    width= 2.5,
    height=3.5)
map.m
dev.off()
pdf("output/figures/supp/FIB_ligand_comp.pdf",
    width= 2.5,
    height=5.5)
map.f
dev.off()

