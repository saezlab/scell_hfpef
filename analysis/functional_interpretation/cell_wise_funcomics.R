## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-10-05
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## add funcomics per cell


library(Seurat)
library(tidyverse)

source(file = "analysis/misc/add_funcomics.R")

seu_int = readRDS( file = "output/seu.objs/study_integrations/harmony_fib_filt.rds")


seu_int= add_tf_activities(seu_int, species= "mouse")

DefaultAssay(seu_int)= "RNA"

seu_int= add_path_activities(seu_int, species = "mouse")
prog.M= GetAssayData(seu_int, assay= "progeny")

saveRDS(prog.M, "output/seu.objs/study_integrations/prog.M.rds")

prog.M= readRDS( "output/seu.objs/study_integrations/prog.M.rds")
seu_int= readRDS("output/seu.objs/study_integrations/harmony_fib_filt_funcomics.rds")


# plot features -------------------------------------------------------------------------------

DefaultAssay(seu_int)= "progeny"

seu_int
pathways= rownames(GetAssayData(seu_int,assay = "progeny"))

progeny_plot= FeaturePlot(seu_int, features = pathways)

pdf("output/figures/funcomics_full.pdf",
    width= 20, height = 15)
progeny_plot
dev.off()


DefaultAssay(seu_int)= "RNA"

FeaturePlot(seu_int, c("Nox2", "Nox", "Nox1"))

# TF check:
DefaultAssay(seu_int)= "dorothea"

genesoi= c("Lrp1", "Tmem219")



# explore TF space within cell type -----------------------------------------------------------

cell_oi <- seu_int@meta.data[, "celltype"] %in% c("Endothelial")

# subset based on filtering and quickly get characteristic profile
seu_int <- seu_int[ , cell_oi]

DefaultAssay(seu_int)= "dorothea"


cells_ct= rownames(seu_int[[]])[seu_int[[]]$group == "ct"]
cells_hf=rownames(seu_int[[]])[seu_int[[]]$group == "hf"]

tf_act= GetAssayData(seu_int)
rownames(tf_act)

colnames(tf_act)

FeaturePlot(seu_int, features = "Hif1a")



# correlate pathways with angptl4 -------------------------------------------------------------

studies= unique(seu_int@meta.data$study)

prog.M= GetAssayData(seu_int, assay= "progeny")

pls= map(studies, function(x){
  cell_oi <- seu_int@meta.data[, "study"] %in% x
  seu_int2 <- seu_int[ , cell_oi]
  prog.M= GetAssayData(seu_int2, assay= "progeny")
  #TF.M = GetAssayData(seu_int, assay= "dorothea")
  gex.M= GetAssayData(seu_int2,assay = "RNA", slot= "data")
  print(paste0("cell ids agree:", all(colnames(gex.M)== colnames(prog.M))))


  cor.tests= apply(prog.M, 1, FUN = function(x){
    res= cor.test(x,gex.M["Angptl4", ] )
    c("p"= res$p.value, res$estimate)
  })
  cor.df.p= t(cor.tests)%>%as.data.frame()%>% rownames_to_column("pathway")%>%
    as_tibble()

  p.pathways= cor.df.p%>% ggplot(., aes(x= reorder(pathway, cor), y = cor, col= -log10(p)))+
    geom_point()+
    scale_color_gradient(low= "grey", high = "darkred")+
    theme(axis.text.x = element_text(angle= 90, hjust= 1))

  p.pathways

})
studies
cell_oi <- seu_int@meta.data[, "celltype"] %in% c("fibroblasts")

# subset based on filtering and quickly get characteristic profile
seu_int2 <- seu_int[ , cell_oi]

DefaultAssay(seu_int)
prog.M= GetAssayData(seu_int, assay= "progeny")
TF.M = GetAssayData(seu_int, assay= "dorothea")
gex.M= GetAssayData(seu_int,assay = "RNA", slot= "data")

table(colnames(gex.M)== colnames(prog.M))
dim(prog.M)
dim(gex.M)


VlnPlot(object = seu_int, features = "Angptl4", group.by= "study",split.by= "group")

FeaturePlot(seu_int, "Angptl4")
DimPlot(seu_int)

map(studies, function(x){
  cell_oi <- seu_int@meta.data[, "study"] %in% x
  seu_int2 <- seu_int[ , cell_oi]
  #prog.M= GetAssayData(seu_int2, assay= "progeny")
  #TF.M = GetAssayData(seu_int, assay= "dorothea")
  gex.M= GetAssayData(seu_int2,assay = "RNA", slot= "data")

  hist(gex.M["Angptl4", ])

})
cor.tests= apply(prog.M, 1, FUN = function(x){
  res= cor.test(x,gex.M["Angptl4", ] )
  c("p"= res$p.value, res$estimate)
})
cor.df.p= t(cor.tests)%>%as.data.frame()%>% rownames_to_column("pathway")%>%
  as_tibble()

p.pathways= cor.df.p%>% ggplot(., aes(x= reorder(pathway, cor), y = cor, col= -log10(p)))+
  geom_point()+
  scale_color_gradient(low= "grey", high = "darkred")+
  theme(axis.text.x = element_text(angle= 90, hjust= 1))

p.pathways


# compare angptl4 positive cells --------------------------------------------------------------


# define a cut off for angptl4 expression and compare progeny scores.
seu_int = readRDS( file = "output/seu.objs/integrated_cellstate.rds")

seu_ = readRDS( file = "output/seu.objs/integrated_funcomics.rds")


seu.objs=  readRDS("output/seu.objs/cell.state.obj.list.rds")
seu= seu.objs$fibroblasts
seu= add_path_activities(seu, species = "mouse")

progeny_scores = progeny::progeny(expr = as.matrix(seu[["RNA"]]@data), scale=TRUE,
                                  organism="Mouse", top=500, perm=1,z_scores= F)
rownames(progeny_scores)= str_replace_all(rownames(progeny_scores), pattern = "\\.", replacement = "-")
seu[['progeny']] = CreateAssayObject(counts = t(progeny_scores))

prog.M2= GetAssayData(seu,assay = "progeny")

VlnPlot(seu, "Angptl4",group.by = "cellstate")

gex.M= GetAssayData(seu,assay = "RNA", slot= "data")

#check distribution:
hist(gex.M["Angptl4",], breaks=100)
cutoff=1


cell.pos= colnames(gex.M)[gex.M["Angptl4",]>cutoff]
cell.neg= colnames(gex.M)[gex.M["Angptl4",]<=cutoff]

Idents(seu)= "group"
cell.neg= CellsByIdentities(seu, idents= "ct")
cell.pos= CellsByIdentities(seu, idents= "hf")

plot.tab= prog.M2 %>%as.data.frame() %>% rownames_to_column("pathway")%>% as_tibble() %>%
  pivot_longer(names_to = "cellid", values_to= "p.score", -pathway)%>%
  mutate("Angptl4+"= ifelse(cellid %in% cell.pos,"y",
                            ifelse(cellid %in% cell.neg, "n", NA)
                            )
         )%>%
  filter(!is.na('Angptl4+'))

ps= sapply(unique(plot.tab$pathway), function(x){
  wilcox.test(p.score ~`Angptl4+`, data= plot.tab%>% filter(pathway==x))$p.value

})

names(ps[ps<0.05])

enframe(ps)%>% mutate(logp= -log10(value))%>%
  ggplot(., aes(x= reorder(name,logp),  y= logp))+
  geom_col()+
  geom_hline(yintercept = -log10(0.05))

ggplot(plot.tab%>% filter(pathway %in%names(ps[ps<0.05])), aes( x = `Angptl4+`, y= p.score))+
  geom_violin()+
  geom_boxplot(width= 0.3)+
  facet_grid(. ~ pathway, scales = "free")
  stat_compare_means()



prog.M= GetAssayData(seu_int2, assay= "progeny")
  #TF.M = GetAssayData(seu_int, assay= "dorothea")
  gex.M= GetAssayData(seu_int2,assay = "RNA", slot= "data")
  print(paste0("cell ids agree:", all(colnames(gex.M)== colnames(prog.M))))

  prog.M= prog.M2
  cor.tests= apply(prog.M, 1, FUN = function(x){
    res= cor.test(x,gex.M["Angptl4", ] )
    c("p"= res$p.value, res$estimate)
  })
  cor.df.p= t(cor.tests)%>%as.data.frame()%>% rownames_to_column("pathway")%>%
    as_tibble()

  p.pathways= cor.df.p%>% ggplot(., aes(x= reorder(pathway, cor), y = cor, col= -log10(p)))+
    geom_point()+
    scale_color_gradient(low= "grey", high = "darkred")+
    theme(axis.text.x = element_text(angle= 90, hjust= 1))

  p.pathways

  prog.M2
  table(colnames(gex.M)== colnames(prog.M2))
  comb= rbind(prog.M2 ,  "Angptl4"= gex.M["Angptl4",])
  comb= t(comb)
  plot(comb[,"Hypoxia"],comb[,"Angptl4"])

##TF:

  cor.tests= apply(TF.M, 1, FUN = function(x){
    res= cor.test(x,gex.M["Angptl4", ] )
    c("p"= res$p.value, res$estimate)
  })
cor.df= t(cor.tests)%>%as.data.frame()%>% rownames_to_column("tf")%>%
  as_tibble()

x= dorothea::dorothea_mm
angptl4tfs= x %>% filter(target =="Angptl4") %>% pull(tf)

p.TFs= cor.df%>% filter(p<0.001) %>%
  mutate(reg= ifelse(tf %in% angptl4tfs, "reg", "not.reg"))%>%
  ggplot(., aes(x= reorder(tf, cor), y = cor, col= -log10(p)))+
  geom_point(aes(size=reg))+
  scale_color_gradient(low= "grey", high = "darkred")+
  theme(axis.text.x = element_text(angle= 90, hjust= 1))

cowplot::plot_grid(p.pathways, p.TFs, ncol = 1)

