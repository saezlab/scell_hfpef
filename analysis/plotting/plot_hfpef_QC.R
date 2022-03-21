## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-03-15
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Plot QC metrics for supplement


library(Seurat)
library(tidyverse)

# full data set
seu= readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")

source("analysis/utils.R")

meta_seu= seu[[]]

p1= meta_seu %>%
  ggplot(aes(x= nCount_RNA, y= nFeature_RNA))+
  geom_point() +
  theme_classic()

p2= meta_seu %>%
  ggplot(aes(x= celltype, y= doublet_score))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))

p3= meta_seu %>%
  ggplot(aes(x= celltype, nCount_RNA))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))


p4= meta_seu %>%
  ggplot(aes(x= celltype, nFeature_RNA))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))


p5= meta_seu %>%
  ggplot(aes(x= celltype, y= dissociation_s1))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))

p6= meta_seu %>%
  ggplot(aes(x= celltype, y=percent.mt))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))

p.c1= cowplot::plot_grid(p1, p2, p3, p4, p5,p6, labels = c("F", "G", "H", "I", "J", "K"))


## by sample

p6= meta_seu %>%
  ggplot(aes(x= orig.ident, y= doublet_score))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))+
  theme_classic()

p7= meta_seu %>%
  ggplot(aes(x= orig.ident, nCount_RNA))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))+
  theme_classic()

p8= meta_seu %>%
  ggplot(aes(x= orig.ident, nFeature_RNA))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))+
  theme_classic()

p9= meta_seu %>%
  ggplot(aes(x= orig.ident, y= dissociation_s1))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))+
  theme_classic()

p10=  meta_seu %>%
  ggplot(aes(x= orig.ident))+
  geom_bar(stat="count")+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))+
  theme_classic()+
  labs(y= "cell count")

p.c2= cowplot::plot_grid(p6, p7, p8,p9,  p10, labels = "AUTO", nrow = 1)
p.c2

p.c=
cowplot::plot_grid(p.c2, p.c1, nrow = 2, rel_heights = c(1,1.7))

pdf(file= "output/figures/QCs/data_overview.pdf",
    height= 10,
    width= 10)
p.c
dev.off()
