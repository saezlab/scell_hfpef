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
## try to interpret hfpef signatures with impc gene sets:


library(tidyverse)
library(Seurat)

hfpef.fib= readRDS(file = "output/cell_specific_obj/fibroblasts.rds")
hfpef.fib= seu

hfpef= readRDS("output/seu.objs/integrated_cellstate_noECs.rds")
impc= read.csv("data/prior_knowledge/IMPC_Cardiovascular_System.tsv", sep ='\t') %>% as_tibble
impc

genesets= split(impc$Gene, f= impc$Phenotype)

#remove eye related sets:
genesets= genesets[!grepl("eye", names(genesets))]
genesets= genesets[!grepl("retin", names(genesets))]
genesets= genesets[!grepl("placen", names(genesets))]
genesets= genesets[!grepl("corne", names(genesets))]
genesets= genesets[!grepl("vitel", names(genesets))]

length(genesets)
names(genesets)

saveRDS(genesets, "data/prior_knowledge/IMPC_processed.rds")
genesets= readRDS("data/prior_knowledge/IMPC_processed.rds")
de.list= readRDS("output/cell_specific_obj/cell_type_based_DE/cells_DEA.rds")


gr_net=
impc %>% filter(Phenotype %in% names(genesets)) %>%
  mutate(mor = 1,
        likelihood= 1) %>% distinct(Phenotype, Gene, mor, likelihood)
x= "Endothelial"

library(decoupleR)

xx= map(names(de.list), function(x){
  vec= de.list[[x]] %>% pull(avg_log2FC)
  vec= as.matrix(vec)
  rownames(vec)= de.list[[x]]$gene
  df = run_wmean(mat = vec,
          gr_net,.source = "Phenotype", .target = "Gene"
          ) %>% mutate(celltype = x)

})
x= do.call(rbind, xx)%>%filter(p_value<0.05,
                               statistic== "norm_wmean") %>%
  ggplot(., aes(x= source, y= score, fill= celltype))+
  geom_col(position = "dodge")+
  coord_flip()

impc%>% filter(Phenotype== "dilated heart left ventricle") %>% View()
