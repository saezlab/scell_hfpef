## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-02-02
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## plotting functions for DEG


library(tidyverse)
library(Seurat)
source("analysis/utils.R")
intersects= readRDS("output/DEG_downsampled_intersect_hfpef.rds")
intersects$fibroblasts

seu.objs=  readRDS("output/seu.objs/cell.state.obj.list.rds")
seu= seu.objs$fibroblasts
seu2= readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")

fibrosis_genes= c("Col1a1","Col1a2", "Col4a1","Timp1", "Sparc", "Pcolce", "Loxl2", "Loxl1" )
hs= c("Hspa1a", "Hspa1b", "Hspb1")
met= c("Angptl4", "Ace", "Dpep1")
infl=c("Dpep1", "Ifi205", "Ccdc80")

p= VlnPlot(seu, features =fibrosis_genes,
        group.by = "group",
        cols = hfpef_cols, ncol = length(fibrosis_genes),
        combine = F,
        pt.size	= -1)

p= map(p, function(x){
  x +theme(axis.title = element_blank(),
           legend.position = "none",
           plot.title = element_text(size=10))


})

p.fib= cowplot::plot_grid(plotlist = p, ncol = length(p))
pdf("output/figures/DEA/hfpef_signature_vlns.pdf",
    width= 7.5,
    height= 2.5)
p.fib
dev.off()


####
p= VlnPlot(seu, features =infl,
           group.by = "group",
           cols = hfpef_cols, ncol = length(infl),
           combine = F,
           pt.size	= -1)

p= map(p, function(x){
  x +theme(axis.title = element_blank(),
           legend.position = "none",
           plot.title = element_text(size=10))


})

p.hs= cowplot::plot_grid(plotlist = p, ncol = length(p))
pdf("output/figures/DEA/hfpef_signature_vlns2.pdf",
    width= 3,
    height= 2.5)
p.hs
dev.off()


####
p= VlnPlot(seu, features =met,
           group.by = "group",
           cols = hfpef_cols, ncol = length(met),
           combine = F,
           pt.size	= -1)

p= map(p, function(x){
  x +theme(axis.title = element_blank(),
           legend.position = "none",
           plot.title = element_text(size=10))


})

p.infl= cowplot::plot_grid(plotlist = p, ncol = length(p))
pdf("output/figures/DEA/hfpef_signature_vlns3.pdf",
    width= 3,
    height= 2.5)
p.infl
dev.off()



VlnPlot(seu, features =hs,
        group.by = "group",
        cols = hfpef_cols, ncol = length(hs))
VlnPlot(seu, features =infl,
        group.by = "group",
        cols = hfpef_cols, ncol = length(hs))
VlnPlot(seu, features = c("Angptl4", "Col4a1", "Ace", "Slit3", "Pcoce", "Timp1"), group.by = "group",
        cols = hfpef_cols)
