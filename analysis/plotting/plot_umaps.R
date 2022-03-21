## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-16
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## cell state marker plot



library(Seurat)
library(tidyverse)
library(WriteXLS)

source("analysis/utils.R")

# load sc data
seu.objs=  readRDS("output/seu.objs/cell.state.obj.list.rds")
seu= seu.objs$fibroblasts
rm(seu.objs)

p.fibs=


pdf("output/figures/cell_type_assignment/fibroblast_cellstate_marker.pdf",
    width = 5,
    height= 5.5)
p.marker
dev.off()
