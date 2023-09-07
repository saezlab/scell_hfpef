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

###s1pr3 laura

rm(seu_int)

seu= readRDS("output/seu.objs/study_integrations/harmony_fib_filt_2.rds")

p1=  FeaturePlot(seu, c("S1pr3"))
p2= VlnPlot(seu,features = "S1pr3", group.by = c("group", "study"))
unique(seu@meta.data$study)
seu@meta.data= seu@meta.data%>% mutate(exp.group= paste0(study, group, sep= "_"))
p2= VlnPlot(seu,features = "S1pr3", group.by = c("exp.group"))
p3= VlnPlot(seu,features = "S1pr3", group.by = c("opt_clust_integrated"))
p2

library(cowplot)

cowplot::plot_grid(p1,p2,p3, ncol = 2)


# check regulons ------------------------------------------------------------------------------
library(dorothea)
data("dorothea_mm")

dorothea_mm%>% filter(tf== "Ppara") %>% print(n=100)

dorothea_mm%>% filter(tf== "Pparg") %>% print(n=100)
dorothea_mm%>% filter(target=="Angptl4")


# check regen fibs ----------------------------------------------------------------------------

int.fibs= readRDS("output/seu.objs/study_integrations/harmony_fib_filt.rds")

features = c("Col11a1", "Col12a1", "Nppc")
toupper(features)
Idents(int.fibs)= "group"
int.fibs@meta.data= int.fibs@meta.data %>% mutate(study= ifelse(study== "circ", "AngII",
                                                                ifelse(study=="forte", "MI", "HFpEF")))%>%
  mutate(study.group= paste0(study, "_", group))
p1= VlnPlot(int.fibs, features = c("Col11a1", "Col12a1", "Nppc", "Fn1"), group.by = "study.group", ncol = 4)
p2= FeaturePlot(int.fibs, features = c("Col11a1", "Col12a1", "Nppc", "Fn1"))
int.fibs$opt_clust_integrated
DimPlot(int.fibs, group.by = "opt_clust_integrated")
pdf("output/figures/regenerative_markrer.pdf",
    width = 9, height= 9)
cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 1.5))
dev.off()

p1= VlnPlot(subset(int.fibs, study== "hfpef" ), features, group.by = "group")
p2= VlnPlot(subset(int.fibs, study== "forte" ), features, group.by = "group")
p3= VlnPlot(subset(int.fibs, study== "circ" ), features, group.by = "group")

df= readRDS("output/fib_integration/marker_list/DEG_per_study_LOGFC.rds")

map(names(df), function(x){
  df[[x]][features,]
})



# get cellstate raw numbers -------------------------------------------------------------------

dF= map(unique(meta$study), function(x){
  meta %>% filter(study == x)%>% group_by(group, opt_clust_integrated)%>% count%>%
    mutate(study= x)
}) %>% do.call(rbind, .)

pls= map(c("ct", "hf"), function(x){

    dF %>% filter(group ==x) %>%
    ggplot(., aes(x= study, y= opt_clust_integrated, fill = n)+)
    geom_tile()+
    scale_fill_gradient(low= "grey", high = "white")+
    geom_text(aes(label =n))+
    labs(y= "IFS")

})

cowplot::plot_grid(plotlist= pls, ncol = 2)


# angptl4 across cell types -------------------------------------------------------------------

seu = readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")
Idents(seu) = "celltype2"
seu@meta.data= seu@meta.data%>%
  mutate(celltype2= ifelse(celltype2== "Fibroblasts.2" | celltype2 == "Fibroblasts.1",
                           "Fibroblasts",
                           celltype2))
DimPlot(seu, group.by = "celltype2")
p1= VlnPlot(seu, features = "Angptl4", split.by = "group", pt.size = 0.0000, flip = T)
p2= VlnPlot(seu, features = "Angptl4", split.by = "group", pt.size = 0.1, flip = T)
p3= DotPlot(seu, features = "Angptl4", split.by = "group")

pdf("output/figures/main/Fig5/Angptl4_expr.pdf",
    width= 4.5,
    height= 4)
p1
p2
p3

dev.off()
p3
