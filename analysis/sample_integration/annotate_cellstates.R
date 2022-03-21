## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-09-10
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## calculate cell.state marker and manually annotate cell states


library(Seurat)
library(tidyverse)
library(WriteXLS)

cell.types= list.files("output/cell_specific_obj/")
cell.types = str_replace_all(cell.types, ".pdf", "")
cell.types = str_replace_all(cell.types, ".rds", "") %>% unique()
cell.types = cell.types[!grepl(pattern = "statemarker", x=  cell.types)]
cell.types = cell.types[!grepl(pattern = "DE", x=  cell.types)]
seu.obj.paths= list.files("output/cell_specific_obj/", full.names = T)
seu.obj.paths= seu.obj.paths[grepl(".rds", seu.obj.paths)]

# load seurats
seu.objs= map(seu.obj.paths, readRDS)
names(seu.objs) = cell.types

# marker ident --------------------------------------------------------------------------------
# perform DEA

de.list = map(seu.objs, function(x){
  Idents(x)= "opt_state"
  de = FindAllMarkers(x,
                 assay= "RNA",
                  )
})

names(de.list) = str_replace_all(names(de.list), "\\.", "_")

map(names(de.list), function(x){
  print(x)
  #y= str_replace(x, "\\.", "_")
  markers.list= split(x=de.list[[x]], f = de.list[[x]]$cluster)
  WriteXLS(markers.list,
           ExcelFileName = paste0("output/cell_specific_obj/", x, "_statemarker.xlsx")
           )

})

saveRDS(de.list, "output/cell_state_marker.rds")

# label cell states ---------------------------------------------------------------------------

## B.cells:
seu.objs$B.cells@meta.data= seu.objs$B.cells@meta.data %>%
  mutate(cellstate= ifelse(opt_state==0, "B-2.Bcell", opt_state),
         cellstate= ifelse(opt_state==1, "B-1.Bcell", cellstate))

## T.cells=
seu.objs$T.cells@meta.data= seu.objs$T.cells@meta.data %>%
  mutate(cellstate= ifelse(opt_state==0, "CD4.Th", opt_state),
         cellstate= ifelse(opt_state==1, "CD8.T.naive", cellstate),
         cellstate= ifelse(opt_state==2, "CD8.T.cytotoxic", cellstate),
         cellstate= ifelse(opt_state==3, "CD8.T.effector", cellstate),
         cellstate= ifelse(opt_state==4, "T.reg", cellstate))

## Macrophages=
seu.objs$macrophages@meta.data= seu.objs$macrophages@meta.data %>%
  mutate(cellstate= ifelse(opt_state==0, "Ccr2-/ MhcII-/ Ly6c2-", opt_state),
         cellstate= ifelse(opt_state==1, "Ccr2-/ MhcII-/ Lyve1+", cellstate),
         cellstate= ifelse(opt_state==2, "Ly6C++/ Ccr2++ monocytes", cellstate),
         cellstate= ifelse(opt_state==3, "Ccr2-/ MhcII-/ Cxcl2++", cellstate),
         cellstate= ifelse(opt_state==4, "Cccr2+/ MhcII+", cellstate))

## Fibroblasts=
seu.objs$fibroblasts@meta.data= seu.objs$fibroblasts@meta.data %>%
  mutate(cellstate= ifelse(opt_state==0, "Col15a1+", opt_state),
         cellstate= ifelse(opt_state==1, "Igfbp3+", cellstate),
         cellstate= ifelse(opt_state==2, "Pi16+", cellstate),
         cellstate= ifelse(opt_state==3, "Cxcl1+", cellstate),
         cellstate= ifelse(opt_state==4, "Cilp+", cellstate),
         cellstate= ifelse(opt_state==5, "Wif1+", cellstate))%>%
  mutate(cellstate= factor(cellstate, levels =  c("Col15a1+", "Igfbp3+",
                                                "Pi16+","Cxcl1+","Cilp+","Wif1+")))

## Endothelial.cells
seu.objs$Endothelial@meta.data= seu.objs$Endothelial@meta.data %>%
  mutate(cellstate= ifelse(opt_state==0, "Coronary.vasculature.ECs", opt_state),
         cellstate= ifelse(opt_state==1, "art.Ec", cellstate),
         cellstate= ifelse(opt_state==2, "senescence?", cellstate),
         cellstate= ifelse(opt_state==3, "capillary.Ec", cellstate),
         cellstate= ifelse(opt_state==4, "Endocardial ECs", cellstate)
         )

saveRDS(seu.objs, file = "output/seu.objs/cell.state.obj.list.rds")
seu.objs= readRDS(file = "output/seu.objs/cell.state.obj.list.rds")


# use cell state annotation to move to larger seu.obj -----------------------------------------
cell.types=names(seu.objs)
cell.types = cell.types[!cell.types %in% c("lymphocytes","myeloid.cells")]

x= lapply(cell.types, function(x){
  print(x)
  y= seu.objs[[x]]@meta.data
  as.data.frame(y) %>% rownames_to_column("cellID") %>% select(cellID, cellstate)
})

cellstate.df= do.call(rbind, x)

seu= readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")

seu.meta= seu[[]]
seu.meta= as.data.frame(seu.meta) %>%
  rownames_to_column("cellID") %>%
  select(-cellstate) %>%
  left_join(cellstate.df, by= "cellID")
seu = AddMetaData(seu, seu.meta$cellstate, "cellstate")
seu@meta.data
# now we add a celltate = celltype in the meta file for celltypes we didnt subcluster(
# nk, granulo and smooth muscle
seu@meta.data= seu@meta.data %>% mutate(cellstate = ifelse(is.na(cellstate),celltype, cellstate))
#unique(seu.objs$fibroblasts@meta.data$cellstate)

saveRDS(seu,"output/seu.objs/integrated_cellstate_nameupdated.rds")

