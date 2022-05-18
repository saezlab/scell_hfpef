## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-05-14
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## test for set of DEGs whether rather cell state marker or general response.
library(tidyverse)
library(Seurat)
library(lsr)

fib= readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled2.rds")
seu= readRDS("output/seu.objs/cell.state.obj.list.rds")
seu= seu$fibroblasts
degs= fib$total$hfpef




# Anovas for cell state label ------------------------------------------------------------

#get normalized expr data
log.data= GetAssayData(seu, slot = "data")
# get cell state labels:
Idents(seu)= "cellstate"
labels.state= CellsByIdentities(seu, )
labels.state= enframe(labels.state, name="state", value= "cellid") %>% unnest(cellid)

Idents(seu) = "group"
labels.group= CellsByIdentities(seu)
labels.group= enframe(labels.group, name="group", value= "cellid" )%>% unnest(cellid)

dim(as.matrix(log.data[x,]))

#prepare dataframe to run anovas on
model.df= (log.data[degs,])%>%as.data.frame()%>%
  t() %>%
  as.data.frame()%>%
  rownames_to_column("cellid")%>%
  as_tibble()%>%
  left_join(labels.state, by= "cellid")%>%
  left_join(labels.group, by= "cellid")
#
# x= model.df[, c("Angptl4", "state", "cellid")]
# x= x %>% drop_na
# aov.1=aov(Angptl4~state, data = model.df)
# aov.2=aov(Angptl4~group, data = model.df)
# eta1= etaSquared(aov.1)[1,1]
# eta2= etaSquared(aov.2)[1,1]
# eta1/eta2
# x= degs[1]

#run anovas per gene and extragct etasquared
df.aov= sapply(degs, function(x){
  print(x)

  aov.1=aov(as.formula(paste0(x ,"~ state")), data = model.df)
  aov.2=aov(as.formula(paste0(x ,"~ group")), data = model.df)
  eta1= etaSquared(aov.1)[1,1]
  eta2= etaSquared(aov.2)[1,1]

  c(eta1, eta2, eta1*eta2)


})
rownames(df.aov)= c("eta.state", "eta.group", "sg.ratio")
df.aov= t(df.aov)

library(ComplexHeatmap)
h1= Heatmap(df.aov[,1:2])= Heatmap(df.aov[,3])
aov.tidy= rownames_to_column(as.data.frame(df.aov), "gene")%>% as_tibble()

topstate= aov.tidy%>% arrange(desc(eta.state))%>%pull(gene)
topgroup= aov.tidy%>% arrange((eta.state))%>%pull(gene)

VlnPlot(seu, features = c(topstate[1:5]),split.by = "group")
Idents(seu)= "cellstate"
VlnPlot(seu, features = c(topstate[1:5]),split.by = "cellstate")
VlnPlot(seu, features = c(topgroup[1:5]),split.by = "group")
VlnPlot(seu, features = c(topgroup[1:5]),split.by = "cellstate")
FeaturePlot(seu, features = c(topgroup[1:5]), split.by = "group")


Idents(seu)= "cellstate"

VlnPlot(seu, features = c(topstate[1:6]))

Idents(seu)= "cellstate"
VlnPlot(seu, features = c(topgroup[1:6]))



h2= Heatmap(df.aov[,1], name = "eta.sq_cellstate", column_title = "",column_labels = "")
h2
?Heatmap

fib$overlap$all
df.aov[,1]

library(fgsea)

fgseaSimple(pathways= list("fibrosis_core"= fib$overlap$all,
                           "BM"= naba$Basement_membranes) ,stats= df.aov[,1], nperm= 1000 )
p1= fgsea::plotEnrichment(pathway= fib$overlap$all, stats= df.aov[,1])
p2= fgsea::plotEnrichment(pathway= naba$Basement_membranes, stats= df.aov[,1])

cowplot::plot_grid(p1+ggtitle("shared fibrosis"),
                   p2+ggtitle("BM"))


topgroup
naba= readRDS("data/prior_knowledge/NABA_mouse.rds")
fgseaSimple(pathways= naba ,stats= df.aov[,1], nperm= 1000 )

msigDB= readRDS("data/prior_knowledge/msigDB.mouse_translated.rds")
msigDB$MSIGDB_KEGG

fgsea.res= fgseaSimple(pathways= list("glucose"= msigDB$MSIGDB_REACTOME$REACTOME_GLUCOSE_METABOLISM,
                                      "FA"= msigDB$MSIGDB_KEGG$KEGG_FATTY_ACID_METABOLISM),
                       stats= df.aov[,1], nperm= 2000 )

fgsea.res= fgseaSimple(pathways= msigDB$MSIGDB_REACTOME,
                       stats= df.aov[,1], nperm= 2000 )
fgsea.res%>% as_tibble()%>% arrange(NES, padj)



# metabolism annotation -----------------------------------------------------------------------


