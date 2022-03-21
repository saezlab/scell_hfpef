## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-09-06
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
##  explore integrated and assign cell types

library(tidyverse)
library(Seurat)
library(cowplot)

## read sc data:
#
# cell_study= readRDS("~/R-projects/sc-exploration/output/cell_mi_files/integrateint_cell_MI_annotated.rds")
#
#   all.markers_forteMI= FindAllMarkers(cell_study,
#                             assay = "RNA")
#
# saveRDS(all.markers_forteMI, file = "output/marker_forte.rds")
#rm(cell_study)

all.markers_forteMI= readRDS("output/marker_forte.rds")

## read sc hfpef:
seu= readRDS("output/seu.objs/seu_harmony.rds")
Idents(seu) = seu$opt_clust_integrated

all.markers= FindAllMarkers(seu,  assay = "RNA")

saveRDS(all.markers, file = "output/marker_hfpef_v1.rds")
all.markers= readRDS(file = "output/marker_hfpef_v1.rds")

#save marker list as xls:
library(WriteXLS)
markers.list= split(x=all.markers, f = all.markers$cluster)
WriteXLS(markers.list, ExcelFileName = "output/figures/cell_type_assignment/marker_all_1.xlsx")

#select only upregulated marker:
all.markers= all.markers %>% filter(avg_log2FC>0)
markers.list= split(x=all.markers, f = all.markers$cluster)
WriteXLS(markers.list, ExcelFileName = "output/figures/cell_type_assignment/marker_up_1.xlsx")




# perform ORA with both datasets --------------------------------------------------------------


## get top genesets per cluster from elvira:
reduced.marker= all.markers_forteMI %>% group_by(cluster) %>% filter(avg_log2FC>0)%>%
  top_n(x = .,n = 300, wt = p_val_adj)

genesets= split(x=reduced.marker$gene, f = reduced.marker$cluster)

##get genes from Monicas list:
genesets2 = read.csv("data/prior_knowledge/human_and_mouse_cell_markers - Markers.csv") %>%
  as_tibble() %>% filter(Tissue.of.Origin== "Heart")

genesets2= split(x=genesets2$Mouse.Gene, f = genesets2$Cell.Type)
#genesets2= split(x=genesets2$Human.Gene, f = genesets2$Cell.Type)

saveRDS(genesets2, "data/prior_knowledge/human_and_mouse_cell_markers_processed.rds")

## get top DE genes from our cluster:
stat= all.markers %>% group_by(cluster) %>% filter(avg_log2FC>0,
                                            p_val_adj<0.05)%>%
  top_n(x = .,n =  400, wt = p_val_adj)
top.marker= split(x=stat$gene, f = stat$cluster)

#function for ORA:
GSE_analysis = function(geneList,Annotation_DB){
  library(dplyr)
  library(tidyr)
  library(tibble)

  geneList = geneList[geneList %in% unique(unlist(Annotation_DB))]

  ResultsDF = matrix(0,nrow = length(Annotation_DB),ncol = 5)
  rownames(ResultsDF) = names(Annotation_DB)
  colnames(ResultsDF) = c("GenesInPathway","GenesInList","GeneNames","p_value","corr_p_value")

  DB_genecontent = length(unique(unlist(Annotation_DB)))

  GenesDB = DB_genecontent
  SelectedGenes = length(geneList)

  for(gset in rownames(ResultsDF)){
    GP = length(Annotation_DB[[gset]])
    GL = length(intersect(Annotation_DB[[gset]],geneList))

    ResultsDF[gset,"GenesInList"] = GL
    ResultsDF[gset,"GenesInPathway"] = GP
    ResultsDF[gset,"GeneNames"] = paste(intersect(Annotation_DB[[gset]],geneList),collapse = ",")
    ResultsDF[gset,"p_value"] = phyper(q=GL - 1, m=GP, n=GenesDB-GP, k=SelectedGenes, lower.tail = FALSE, log.p = FALSE)
  }

  ResultsDF[,"corr_p_value"] = p.adjust(ResultsDF[,"p_value"],method = "BH")
  ResultsDF = data.frame(ResultsDF,stringsAsFactors = F)
  ResultsDF = ResultsDF[order(ResultsDF[,"p_value"]),]

  ResultsDF = ResultsDF %>%
    rownames_to_column("gset") %>%
    mutate_at(c("GenesInPathway","GenesInList",
                "p_value","corr_p_value"),
              as.numeric) %>%
    dplyr::arrange(corr_p_value,GenesInList)

  return(ResultsDF)

}

# function to plot ORA results:
plot_ORA= function(Ora_res, top.marker ){

  for (i in names(Ora_res)){
    Ora_res[[i]] = Ora_res[[i]] %>% mutate(cluster= i)
  }

  Ora_res= as_tibble(do.call(rbind, Ora_res))

  p.p.ora= ggplot(Ora_res, aes(x= factor(cluster, levels = seq(0,18,1)),
                               y= gset, fill= -log10(corr_p_value)))+
    geom_tile()+
    scale_fill_gradient2(low= "white" , high= "black")

  return(p.p.ora)
}

Ora_res= lapply(top.marker, function(x){
  GSE_analysis(x, genesets)

})

p.p.ora= plot_ORA(Ora_res, top.marker)

Ora_res2= lapply(top.marker, function(x){
  GSE_analysis(x, genesets2)

})

p.p.ora2= plot_ORA(Ora_res2, top.marker)

pdf("output/figures/cell_type_assignment/ora_with_forte_1.pdf",
    width = 15,
    height= 10)

plot_grid(p.p.ora+ggtitle("forte_genemarker"),
          p.p.ora2+ggtitle("internal_marker_list"),
          DimPlot(seu, reduction = "umap_harmony", label= T, group.by = "opt_clust_integrated")
)
dev.off()

#based on this plot we can assing preliminary cell types:
## Hard coded!

x= seu[[]]
DimPlot(seu, label= T)
x= as_tibble(x) %>%
  mutate(celltype= ifelse(opt_clust_integrated %in% c(0,9), "Fibroblasts",NA),
         celltype= ifelse(opt_clust_integrated %in% c(1,11), "Endothelial",celltype),
         celltype= ifelse(opt_clust_integrated %in% c(2), "B.cells",celltype),
         celltype= ifelse(opt_clust_integrated %in% c(3,6), "T.cells",celltype),
         celltype= ifelse(opt_clust_integrated %in% c(4), "Macrophages",celltype),
         celltype= ifelse(opt_clust_integrated %in% c(12), "Granulocytes",celltype),
         celltype= ifelse(opt_clust_integrated %in% c(13), "SMC/Pericytes",celltype),
         celltype= ifelse(opt_clust_integrated %in% c(5), "NK.cells",celltype),
         celltype= ifelse(opt_clust_integrated %in% c(10), "lowQ?",celltype),
         celltype= ifelse(opt_clust_integrated %in% c(8), "endothelial/fibro?",celltype),
         celltype= ifelse(opt_clust_integrated %in% c(7), "fibro/mito?",celltype)
         ) %>%
  mutate(celltype2= ifelse(opt_clust_integrated %in% c(0), "Fibroblasts",NA),
         celltype2= ifelse(opt_clust_integrated %in% c(9), "Fibroblasts.activated",celltype2),
         celltype2= ifelse(opt_clust_integrated %in% c(1), "Endothelial.1",celltype2),
         celltype2= ifelse(opt_clust_integrated %in% c(11), "Endothelial.2",celltype2),
         celltype2= ifelse(opt_clust_integrated %in% c(2), "B.cells",celltype2),
         celltype2= ifelse(opt_clust_integrated %in% c(3), "T.effector",celltype2),
         celltype2= ifelse(opt_clust_integrated %in% c(6), "T.helper",celltype2),
         celltype2= ifelse(opt_clust_integrated %in% c(4), "Macrophages",celltype2),
         celltype2= ifelse(opt_clust_integrated %in% c(12), "Granulocytes",celltype2),
         celltype2= ifelse(opt_clust_integrated %in% c(13), "SMC/Pericytes",celltype2),
         celltype2= ifelse(opt_clust_integrated %in% c(5), "NK.cells",celltype2),
         celltype2= ifelse(opt_clust_integrated %in% c(10), "lowQ.cells?",celltype2),
         celltype2= ifelse(opt_clust_integrated %in% c(8), "endothelial/fibro?",celltype2),
         celltype2= ifelse(opt_clust_integrated %in% c(7), "fibro+mitochondrial?",celltype2)
  )

ggplot(x, aes(x= opt_clust_integrated, y= percent.mt))+
  geom_point()+
  geom_boxplot()

ggplot(x, aes(x= celltype, y= percent.mt))+
  geom_point()+
  geom_boxplot()

ggplot(x, aes(x= celltype, y= doublet_score))+
  geom_point()+
  geom_boxplot()

# add celltypes to the meta data

seu = AddMetaData(seu, x$celltype, col.name = "celltype")
seu = AddMetaData(seu, x$celltype2, col.name = "celltype2")
DimPlot(seu, group.by = "celltype", label=T
        )
DimPlot(seu, group.by = "celltype2", label=T
)

saveRDS(seu, file = "output/seu.objs/integrated_annotated.rds")


# get cell IDs of problematic cells  ------------------------------------------

meta.seu %>% distinct(opt_clust_integrated, celltype, celltype2) %>% arrange(opt_clust_integrated)%>%
  as_tibble()%>%
  write_csv("output/figures/cell_type_assignment/suggested_annotation.csv")

meta.seu = seu[[]]
p1= calc_props(meta.seu, "opt_clust_integrated", "group")
p2= calc_props(meta.seu, "celltype2", "group")
p3= calc_props(meta.seu, "celltype", "group")

rownames_to_column(meta.seu, var = "cell.id")
