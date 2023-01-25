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
## Prepare subsetting of study atlases to macrophages prior integration between different studies. HFpEF and AngII model

library(tidyverse)
library(Seurat)
library(cowplot)
library(purrr)

# # Load SC objects and unifydate meta data: --------------------------------------------------
# Ang II model from circulation:
directory= "/home/jan/R-projects/sc-exploration/"

circ= readRDS( paste0(directory, "data/circ_MI/integrated_seurat_circMI.rds"))
meta = read_csv(paste0(directory, "data/circ_MI/meta_circMI.csv"))
files = readRDS(paste0(directory, "data/circ_MI/file_meta.rds"))

meta= meta %>% left_join(files %>% select(-file.path) %>% rename(Extract.Name = filename.short))

seu_meta= circ[[]]
seu_meta= seu_meta %>%
  left_join(meta %>% rename(orig.ident= filename.long ), by= "orig.ident")%>%
mutate(group= ifelse(Factor.Value.compound. == "angiotensin II", "hf", "ct"))

# add info to seu, mainly treatment info
circ = AddMetaData(circ, seu_meta$group, col.name = "group")
circ@meta.data= circ@meta.data %>%
  mutate( study= "AngII")

Idents(circ)= "seurat_clusters"
DefaultAssay(circ) = "RNA"

mac.maker=  c("Lyz2", "Cd68", "Cx3cr1")

marker= read.csv("data/Celltypes_marker_plot.csv")
marker = str_replace_all(unlist(str_split(marker$Markers, ",")), " ", "")

p1= DimPlot(circ, group.by = "seurat_clusters",label = T)
p2=DotPlot(circ, features = marker)+
  coord_flip()
p3= FeaturePlot(circ,
            features =  mac.maker
            )

p.c= plot_grid(plot_grid(p1, p3), p2, ncol = 1)

spam_cols <- grepl("snn_res",
                   colnames(circ@meta.data))
circ@meta.data <- circ@meta.data[,!spam_cols]

macs= subset(circ,subset = seurat_clusters %in%  c(2,7,8,15) )

saveRDS(macs, file = "output/seu.objs/study_integrations/macs/circ_macs.rds")

macs

rm(circ)


# MI  -----------------------------------------------------------------------------------------

#meta data is missing
forte= readRDS( file = "../sc-exploration/output/cell_mi_files/fibro_subset_integrated.rds")
forte@meta.data= forte@meta.data %>%
  mutate(group = ifelse(OP== "myocardial infarction", "hf", "ct"),
         study= "forte")

fibs_meta = forte[[]]
rm(forte)

df1= fibs_meta %>% distinct(orig.ident, group, study)

# load full altas and add meta data
forte= readRDS( file = paste0(directory, "output/cell_mi_files/integrated_seurat.rds"))

df2= forte@meta.data %>% left_join(df1, by= "orig.ident")
forte= AddMetaData(forte, df2$study, "study")
forte= AddMetaData(forte, df2$group, "group")

DefaultAssay(forte)= "RNA"

p1= DimPlot(forte, group.by = "seurat_clusters",label = T)
p2=DotPlot(forte, features = marker)+
  coord_flip()
p3= FeaturePlot(forte,
                features =  mac.maker
)


p.c= plot_grid(plot_grid(p1, p3), p2, ncol = 1)

p.c

macs.c= c(0,8,10,11,14,15)

spam_cols <- grepl("snn_res",
                   colnames(forte@meta.data))
forte@meta.data <- forte@meta.data[,!spam_cols]

macs= subset(forte,subset = seurat_clusters %in%  macs.c )

saveRDS(macs, file = "output/seu.objs/study_integrations/macs/mi_macs.rds")

