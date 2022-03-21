## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-10-15
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Plot cell proportions

library(Seurat)
library(tidyverse)
library(ggsignif)


# CELL TYPES ----------------------------------------------------------------------------------


# full data set
seu= readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")

source("analysis/utils.R")

x= seu[[]]
x = x %>% mutate(celltype= ifelse(grepl("Fib", celltype2), "Fibroblasts", celltype2))
p.by.sample= calc_props(x, cluster.col = "orig.ident", "celltype")
p.by.sample2= calc_props(x, cluster.col = "orig.ident", "celltype2")
p.by.group= calc_props(x, cluster.col = "group", "celltype")
p.by.group2= calc_props(x, cluster.col = "group", "celltype2")

df.mean= calc_mean_proportions_per_group(x, "celltype", "group")

df. = df.mean$groupwise%>%
  mutate(celltype  = factor(celltype,levels = c("Fibroblasts",
                                                    "Endothelial",
                                                    "NK.cells",
                                                    "Macrophages",
                                                    "Cd8.T.cell",
                                                "Cd4.T.cell",
                                                    "B.cells",
                                                    "Granulocytes",
                                                    "SMC/Pericytes")))
p.mean= plot_mean_proportions( df., "celltype", "Celltype Composition")

p.mean. = p.mean+
  geom_signif(
    y_position = c(0.6, 0.25, 0.1,0.18),
    xmin = c(0.8, 1.8, 3.8, 6.8),
    xmax = c(1.2, 2.2, 4.2, 7.2),
    annotation = rep("*", 4), tip_length = 0
  ) +
  ylim(c(0,0.7))+
  labs(y= "mean %",
       x= "",
       fill= "")+
  ggtitle("")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(color= "black"))

p.mean.

pdf(file = "output/figures/proportion_per_group_final.pdf",
    width= 6,
    height= 4)
p.mean.
dev.off()


pdf(file = "output/figures/proportion_per_group.pdf",
    width= 6,
    height= 5)
p.by.sample
p.by.sample2
p.by.group
p.by.group2
p.mean

dev.off()



# CELL STATES ---------------------------------------------------------------------------------

seu.objs=  readRDS("output/seu.objs/cell.state.obj.list.rds")

##fibs=
seu= seu.objs$fibroblasts

y= seu[[]]
p.mean= calc_mean_proportions_per_group(y, "cellstate", "group")
p.mean= plot_mean_proportions(p.mean$groupwise, "cellstate", "Cellstate Composition")

p.fibs=p.mean +
  labs(y= "mean %",
       x= "",
       fill= "")+
  geom_signif(
    y_position = c(0.4, 0.2, 0.15), xmin = c(0.8, 3.8, 5.8), xmax = c(1.2, 4.2, 6.2),
    annotation = rep("*", 3), tip_length = 0
  ) +
  ylim(c(0,0.5))+
  ggtitle("")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(color= "black"))

##macs=
seu= seu.objs$macrophages
y= seu[[]]
p.mean= calc_mean_proportions_per_group(y, "cellstate", "group")
p.mean= plot_mean_proportions(p.mean$groupwise, "cellstate", "Cellstate Composition")

pmacs= p.mean +
  labs(y= "mean %",
       x= "",
       fill= "")+
  geom_signif(
    y_position = c(0.28, 0.28, 0.42, 0.4, 0.22),
    xmin = c(0.8, 1.8,2.8, 3.8,4.8),
    xmax = c(1.2,2.2, 3.2, 4.2, 5.2),
    annotation = rep("*",5), tip_length = 0
  ) +
  ylim(c(0,0.5))+
  ggtitle("")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(color= "black"))
pmacs


pdf("output/figures/cell_type_analysis/cellstate_macs_fibs_proportions.pdf",
    width = 5,
    height= 4)
p.fibs
pmacs
dev.off()

