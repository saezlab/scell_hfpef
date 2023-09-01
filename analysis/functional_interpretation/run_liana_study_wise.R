## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2023-01-26
##
## Copyright (c) Jan D. Lanzer, 2023
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## run liana on the case-control contrast study wise and compare macrophage-fibroblast interactions


library(liana)
library(tidyverse)
library(magrittr)
library(Seurat)
library(ggrepel)
library(ComplexHeatmap)
library(biomaRt)

source("code/utils.R")
source("code/utils_funcomics.R")
source('~/R-projects/sc-exploration/R-scripts/utils.R')

directory= "/home/jan/R-projects/sc-exploration/"
gene_intersect= readRDS("output/gene_intersect.rds")

# prepare liana wrap --------------------------------------------------------------------------

op_resource <- select_resource("OmniPath")[[1]] %>%
  convert_to_murine()

## or load if biomart throws errors:
op_= readRDS("~/R-projects/sc_hfpef/output/figures/funcomics/liana_res_celltype.rds")
op_resource= op_[[1]]


aggregate_l= function(liana_res){
  liana_res= liana_res[names(liana_res)!="cellchat"]
  liana_res %<>% liana::liana_aggregate()
}

combine_liana_results= function(source.c,
                                target.c,
                                liana_1,
                                liana_2){

  df1= liana_1 %>% filter(source==source.c,
                          target== target.c) %>%
    arrange((median_rank))%>%
    dplyr::mutate(new_rank_ct= rank(aggregate_rank, ties.method = "average"))%>%
    dplyr::select(source, ligand, target, receptor, aggregate_rank, new_rank_ct)%>%
    dplyr::rename(aggregate_rank_ct = aggregate_rank)

  df2= liana_2 %>% filter(source== source.c,
                          target== target.c) %>%
    arrange((aggregate_rank))%>%
    mutate(new_rank_hf= rank(aggregate_rank, ties.method = "average"))%>%
    dplyr::select(source, ligand, target, receptor, aggregate_rank, new_rank_hf)%>%
    dplyr::rename(aggregate_rank_hf = aggregate_rank)

  comb.df= df2 %>%
    inner_join(., df1, by= c("ligand", "receptor", "source", "target"))%>%
    dplyr::mutate(hf_score= (1-aggregate_rank_hf)*aggregate_rank_ct)%>%
    arrange(desc(hf_score))

}


run_liana_case_control= function(seu, group_control= "ct", group_case= "hf",
                                 feats=NULL){
  if(is.null(feats)){
    feats= readRDS("output/gene_intersect.rds")
  }

  seu_ct= subset(seu, group== group_control, features = feats)
  seu_hf= subset(seu, group== group_case, features = feats)
  rm(seu)

  ## run LIANA between celltypes

  liana_res_ct <- liana_wrap(seu_ct,
                             cellchat.params=list(organism="mouse"),
                             resource = "custom",
                             external_resource = op_resource)# ,
                             #method =c("natmi", "connectome", "logfc")

  liana_res_hf <- liana_wrap(seu_hf,
                             cellchat.params=list(organism="mouse"),
                             resource = "custom",
                             external_resource = op_resource#,
                             #method =c("natmi", "connectome", "logfc") )
  )

  message("aggregating")

  liana_1= aggregate_l(liana_res_ct)
  liana_2= aggregate_l(liana_res_hf)

  message("combining")

  mac.fibs= combine_liana_results(source.c = "macrophages",
                                  target.c= "fibroblasts",
                                  liana_1 = liana_1,
                                  liana_2= liana_2)

  return(list("liana_ct"= liana_res_ct,
              "liana_hf"= liana_res_hf,
              "mac.fibs"= mac.fibs))

}

##load each study and run liana

#Ang II  ---------------------------------------------------------------

directory= "/home/jan/R-projects/sc-exploration/"
seu= readRDS( paste0(directory, "data/circ_MI/integrated_seurat_circMI2.rds"))
seu_meta= seu[[]]

mac.maker=  c("Lyz2", "Cd68", "Cx3cr1")
Idents(seu)= "integrated_snn_res.0.1"
p1= DimPlot(seu, group.by = "integrated_snn_res.0.1",label = T)
p2=DotPlot(seu, features = mac.maker)
p3= FeaturePlot(seu,
                features =  mac.maker
)

seu@meta.data= seu@meta.data %>%
  mutate(celltype= ifelse(integrated_snn_res.0.1 %in% c(0,8), "fibroblasts", ""),
         celltype= ifelse(integrated_snn_res.0.1 %in% c(1), "macrophages", celltype))%>%
  mutate(group = ifelse(treatment== "angiotensin II", "hf", "ct"),
         study= "angii")

Idents(seu)= "celltype"
DefaultAssay(seu)= "RNA"
seu2= subset(seu, celltype %in% c("macrophages", "fibroblasts"))
gene_angii= rownames(seu)

angii_l= run_liana_case_control(seu = seu2)

saveRDS(angii_l, "output/liana_angii.rds")

# MI ------------------------------------------------------------------------------------------


int_cell_MI = readRDS(file= paste0(directory, "output/cell_mi_files/integrateint_cell_MI.rds"))

## add meta data about
mac.maker=  c("Lyz2", "Cd68", "Cx3cr1")

marker= read.csv("data/Celltypes_marker_plot.csv")
marker = str_replace_all(unlist(str_split(marker$Markers, ",")), " ", "")

p1= DimPlot(int_cell_MI, group.by = "integrated_snn_res.0.1",label = T)
p2=DotPlot(int_cell_MI, features = marker)+
  coord_flip()
p3= FeaturePlot(int_cell_MI,
                features =  mac.maker
)

# add celltype group
int_cell_MI@meta.data= int_cell_MI@meta.data %>%
  mutate(celltype= ifelse(integrated_snn_res.0.1 %in% c(0,5,6), "fibroblasts", ""),
         celltype= ifelse(integrated_snn_res.0.1 %in% c(1, 8), "macrophages", celltype))%>%
  mutate(group = ifelse(OP== "myocardial infarction", "hf", "ct"),
         study= "forte")

Idents(int_cell_MI)= "celltype"
DefaultAssay(int_cell_MI)= "RNA"
int_cell_MI= subset(int_cell_MI,  celltype %in% c("macrophages", "fibroblasts"))

genes_mi= rownames(int_cell_MI)

#now create two subset, (early and late MI)
forte_early= subset(x = int_cell_MI, subset = time != 14 )
forte_late= subset(int_cell_MI, time %in% c(0,14) | group == "ct" )
rm(int_cell_MI)


mi_e= run_liana_case_control(seu = forte_early)
saveRDS(mi_e, "output/liana_mi_e.rds")

mi_l= run_liana_case_control(seu = forte_late)
saveRDS(mi_l, "output/liana_mi_l.rds")

mi_e$mac.fibs


# hfpef ---------------------------------------------------------------------------------------
seu= readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")
liana_res= readRDS("output/liana_combined_results.rds")

genes_hfpef= rownames(seu)

ggvenn::ggvenn(list("mi"=genes_mi, "ang"= gene_angii, "hfpef"= genes_hfpef))
gene_intersect= intersect(intersect(genes_hfpef, genes_mi), gene_angii)
length(gene_intersect)
saveRDS(gene_intersect, "output/gene_intersect.rds")
