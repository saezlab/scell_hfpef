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
##  After filtering and reintegration we are labeling celltypes and states after supervised marker inspection

library(tidyverse)
library(Seurat)
library(cowplot)

## read sc hfpef:
seu= readRDS("output/seu.objs/integrated_annotated_filtered.rds")
Idents(seu) = seu$celltype2

## Hard coded!

x= seu[[]]
x= as_tibble(x) %>%
  mutate(celltype2= str_replace_all( celltype2,"T.cells.1", "Cd8.T.cell"),
         celltype2= str_replace_all(celltype2,"T.cells.2", "Cd4.T.cell"),
         celltype2= str_replace_all(celltype2,"Endothelial.2", "mesenchymal.endothelial"),
         celltype2= str_replace_all(celltype2,"Endothelial.1", "endothelial"))
x$celltype2
seu = AddMetaData(seu, x$celltype2, col.name = "celltype2")
DimPlot(seu, group.by = "celltype", label=T
        )
DimPlot(seu, group.by = "celltype2", label=T
)

## also add group label!
seu@meta.data=seu@meta.data%>%
  mutate(group = ifelse(grepl("ct", orig.ident),
                        "ct","hfpef"))

saveRDS(seu, file = "output/seu.objs/integrated_annotated_filtered.rds")


# get cell IDs of problematic cells  ------------------------------------------
p1= calc_props(x, "opt_clust_integrated", "group")
p2= calc_props(x, "celltype2", "group")
p3= calc_props(x, "celltype", "group")

rownames_to_column(meta.seu, var = "cell.id")
