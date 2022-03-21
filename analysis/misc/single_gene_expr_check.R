/## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-10-05
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## quickly check gene expression


library(Seurat)
library(tidyverse)

cell.types= list.files("output/cell_specific_obj/")
cell.types = str_replace_all(cell.types, ".pdf", "")
cell.types = str_replace_all(cell.types, ".rds", "") %>% unique()
cell.types = cell.types[!grepl(pattern = "statemarker", x=  cell.types)]
cell.types = cell.types[!grepl(pattern = "DE", x=  cell.types)]
seu.obj.paths= list.files("output/cell_specific_obj/", full.names = T)
seu.obj.paths= seu.obj.paths[grepl(".rds", seu.obj.paths)]

# load seurats
seu.objs= map(seu.obj.paths, readRDS)
names(seu.objs) = cell.types


seu_int = readRDS( file = "output/seu.objs/integrated_cellstate_nameupdated.rds")

seu.objs= readRDS(file ="output/seu.objs/cell.state.obj.list.rds")
# marker ident --------------------------------------------------------------------------------

fibs = seu.objs$fibroblasts
DimPlot(fibs, pt.size = 3)

FeaturePlot(fibs, features = "Fap", split.by = "group")
Idents(fibs)= "cellstate"
VlnPlot(fibs, features = "Fap", split.by = "group")
##FAP

fibs= readRDS("output/cell_specific_obj/fibroblasts.rds")

pdf("plot.fap.pdf")
DimPlot(seu_int, pt.size = 1, label= T)
gset_myo(seu_int, c("Fap", "Postn"))
DimPlot(fibs, label =T)
FeaturePlot(fibs, "Fap")
VlnPlot(fibs, c("Fap", "Postn", "Meox1", "Ltbp2"))
dev.off()

# Check genes for Florian --------------------------------------------------------------------------


p1=  FeaturePlot(seu_int, c("Lrp1", "Tmem219"))
p2= DimPlot(seu_int, group.by = "celltype",label= T)

seu.meta= seu_int[[]]
x= GetAssayData(seu_int, assay = "RNA", slot= "data")
y= as.data.frame(x[genesoi, ])

feat= seu.meta[colnames(y), "celltype"]

df= cbind(t(y), feat)
df= as_tibble(df) %>% mutate(Lrp1 = as.numeric(Lrp1),
                             Tmem219 = as.numeric(Tmem219))

p.lrp = df%>%
  ggplot(., aes(x= feat, y= Lrp1))+
  geom_boxplot()+
  labs(x= "celltype",
       y= "Lrp1 (log cpm)")

p.tmem =df%>%
  ggplot(., aes(x= feat, y= Tmem219))+
  geom_boxplot()+
  labs(x= "celltype",
       y= "Tmem219 (log cpm)")

p_comb1= cowplot::plot_grid(p2, p1, nrow= 2)
p_comb2= cowplot::plot_grid(p.lrp, p.tmem, align = "v", ncol = 1)

pdf(file= "output/figures/DEA/lrp_tmem_expression_overview.pdf",
    height= 15)
cowplot::plot_grid(p_comb1,
                   p_comb2, ncol = 1, rel_heights = c(1,0.8))
dev.off()


# check regulons ------------------------------------------------------------------------------
library(dorothea)
data("dorothea_mm")

dorothea_mm%>% filter(tf== "Ppara") %>% print(n=100)

dorothea_mm%>% filter(tf== "Pparg") %>% print(n=100)
dorothea_mm%>% filter(target=="Angptl4")
