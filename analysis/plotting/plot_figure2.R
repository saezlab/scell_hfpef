heade
library(cowplot)
library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
source("code/utils.R")
umap1 = readRDS( "output/figures/main/umap_fibs.rds")


ecm= readRDS("output/figures/main/Fig2/ecm.ora.rds")
ecm =unify_axis(ecm)+
  coord_equal()+
  theme(axis.title = element_blank())

progeny= readRDS("output/figures/main/Fig2/progeny.rds")
props= readRDS("output/figures/main/Fig2/props.rds")

mod.sc= readRDS("output/figures/main/Fig2/modscores.rds")

x= grid.grabExpr(draw(progeny[[2]][[2]]))

x= grid.grabExpr(draw(p.state[[2]], padding = unit(c(0, 0, 12, 0), "mm")))
unify_axis(props)
props

props= props+   theme(axis.text.x= element_text(size= 10, color = "black"),
                      axis.title= element_text(size= 9, color= "black"),
                      legend.title = element_text(size= 10),
                      legend.text = element_text(size= 9),
                      plot.margin = unit(c(0, 0, 0, 0), "cm")
                      )

p_row1= plot_grid(unify_axis(umap1)+
            theme(#axis.title = element_blank(),
                  axis.text= element_blank()),
          ecm,

          x,
          NULL,NULL,
          props,
          ncol = 3, align= "h", axis= "t",
          rel_widths = c(0.8,1,0.8),
          rel_heights = c(1, 0.8))

p_row1

x2= grid.grabExpr(draw(mod.sc[[1]]), padding = unit(c(0, 0, 12, 0), "mm"))
p.degs= mod.sc[[2]]
p.degs+theme(axis.text.x.bottom = element_text(size= 2))

p.bottom.row= plot_grid( p.degs+theme(axis.text.x.bottom = element_text(size= 2)), plot_grid(NULL,  x2),
          ncol = 1, rel_heights = c(1, 0.5))

pdf("output/figures/main/Fig2/bottom.pdf",
    width= 6,
    height= 4)
p.bottom.row
dev.off()

p_scores= plot_grid(NULL, NULL, p.degs+theme(axis.text.x.bottom = element_text(size= 2)),NULL, NULL, plot_grid(NULL,  x2),
                    ncol = 3, rel_heights = c(1, 0.5))
p_scores


p_all= plot_grid(p_row1,NULL, p_scores, ncol= 1)

pdf("output/figures/main/Fig2/top_row.pdf",
    width= 13,
    height= 7)
p_row1
dev.off()
pdf("output/figures/main/Fig2/full1.pdf",
    width= 13,
    height= 15)
p_all
dev.off()
