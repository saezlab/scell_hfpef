## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2023-01-25
##
## Copyright (c) Jan D. Lanzer, 2023
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
## check for background genes (use fib signatures to see if they separte other cell types)
##

library(Seurat)
library(tidyverse)
library(pROC)
library(ComplexHeatmap)

source("code/utils.R")

## fib sigs
gene_signatures= readRDS( "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
gene_signatures$total= gene_signatures$total[names(gene_signatures$total)!= "MI"]

## hfpef study obj

seu = readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")
DimPlot(seu)

# add geneset scores for full data set --------------------------------------------------------

seu= add_geneset_scores(seu, gene_signatures$total)

meta= seu@meta.data

# loop over celltypes and calculate auroc based on different signature scores

celltypes= unique(meta$celltype)
genescores= paste0( names(gene_signatures$total),"1")

res= sapply(celltypes, function(x){
  sapply(genescores, function(y){
    df.test= meta %>% dplyr::filter(celltype==x)
    f.auc= pROC::auc( df.test$group, df.test[[y]])
    as.numeric(f.auc)
  })
})



Heatmap(res, name = "AUROC")

# most celltypes are not seperable, so assuming background should be uniform is unlikely to be background
