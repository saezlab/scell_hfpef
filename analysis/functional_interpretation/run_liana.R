## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-15
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## run nichnet
library(liana)
library(tidyverse)
library(biomaRt)
library(magrittr)
library(Seurat)
library(ggrepel)
library(ComplexHeatmap)

source("code/utils.R")

source('~/R-projects/sc-exploration/R-scripts/utils.R')
# set up prior K ------------------------------------------------------------------------------

# load data -----------------------------------------------------------------------------------

#full seu:
seu= readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")
Idents(seu)= "celltype"

#celltype seu:
seu.obj= readRDS("output/seu.objs/cell.state.obj.list.rds")
macs= seu.obj$macrophages
fibs= seu.obj$fibroblasts
rm(seu.obj)

state_marker= readRDS("output/cell_state_marker.rds")
fib_sigs=readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")



#get O-path prior knowledge and translate to mouse:

op_resource <- select_resource("OmniPath")[[1]] %>%
  convert_to_murine()

liana_res <- liana_wrap(seu,
                        cellchat.params=list(organism="mouse"),
                        resource = "custom",
                        external_resource = op_resource)

saveRDS(list(op_resource, liana_res), "output/figures/funcomics/liana_res_celltype.rds")

op_= readRDS("~/R-projects/sc_hfpef/output/figures/funcomics/liana_res_celltype.rds")
op_resource= op_[[1]]

# liana ---------------------------------------------------------------------------------------
# we run liana in Control mice and HFpEF mice separately
unique(Idents(seu))

seu_ct= subset(seu, group== "ct")
seu_hf= subset(seu, group== "hfpef")
rm(seu)


## run LIANA between celltypes
liana_res_ct <- liana_wrap(seu_ct,
                        cellchat.params=list(organism="mouse"),
                        resource = "custom",
                        external_resource = op_resource# ,
                        #method =c("natmi", "connectome", "logfc")

                        )
saveRDS(liana_res_ct, "output/liana_ct.rds")

liana_res_hf <- liana_wrap(seu_hf,
                           cellchat.params=list(organism="mouse"),
                           resource = "custom",
                           external_resource = op_resource#,
                           #method =c("natmi", "connectome", "logfc") )
)



saveRDS(liana_res_hf, "output/liana_hf.rds")


liana_res_ct= readRDS("output/liana_ct.rds")
liana_res_hf= readRDS("output/liana_hf.rds")

liana_res_hf$natmi
## Aggregate Liana Results:
aggregate_l= function(liana_res){
  liana_res= liana_res[names(liana_res)!="cellchat"]
  liana_res %<>% liana::liana_aggregate()
  }

liana_1= aggregate_l(liana_res_ct)
liana_2= aggregate_l(liana_res_hf)

# filter for macs to fibs and rank based on aggregated rank
combine_liana_results= function(source.c,
                                target.c,
                                liana_1,
                                liana_2){

  df1= liana_1 %>% filter(source==source.c,
                          target== target.c) %>%
    arrange((median_rank))%>%
    mutate(new_rank_ct= rank(aggregate_rank, ties.method = "average"))%>%
    select(source, ligand, target, receptor, aggregate_rank, new_rank_ct)%>%
    rename(aggregate_rank_ct = aggregate_rank)

  df2= liana_2 %>% filter(source== source.c,
                          target== target.c) %>%
    arrange((aggregate_rank))%>%
    mutate(new_rank_hf= rank(aggregate_rank, ties.method = "average"))%>%
    select(source, ligand, target, receptor, aggregate_rank, new_rank_hf)%>%
    rename(aggregate_rank_hf = aggregate_rank)

  comb.df= df2 %>%
    inner_join(., df1, by= c("ligand", "receptor", "source", "target"))%>%
    mutate(hf_score= (1-aggregate_rank_hf)*aggregate_rank_ct)%>%
    arrange(desc(hf_score))

  }

mac.fibs= combine_liana_results(source.c = "macrophages",
                                target.c= "fibroblasts",
                                liana_1 = liana_1,
                                liana_2= liana_2)

fibs.mac= combine_liana_results(target.c = "macrophages",
                                source.c= "fibroblasts",
                                liana_1 = liana_1,
                                liana_2= liana_2)

saveRDS(list("fib.mac"= fibs.mac,"mac.fibs"= mac.fibs), "output/liana_combined_results.rds")
liana_res= readRDS("output/liana_combined_results.rds")


# Prioritize L-R pairs ----------------------------------------------------------------------
## 1. Selecting by HF-score
## 2. Effect size of Ligand upregulation and expression percentage


## enrich Ligands in
# source("analysis/utils.R")
# GSE_analysis(L, Annotation_DB = fib_sigs$total)

gex.fib= FoldChange(subset(seu, idents = c("Fibroblasts")),
                    ident.1= "hfpef",
                    ident.2= "ct",
                    group.by= "group")

gex.mac= FoldChange(subset(seu, idents = c("Macrophages")),
                    ident.1= "hfpef",
                    ident.2= "ct",
                    group.by= "group")

# we focus on macrophage -> fibroblast activation in HFpEF
top_LR = liana_res$mac.fibs%>%
  filter(hf_score >0.9)%>%
  arrange(desc(hf_score))#%>%
  #top_n(n= 40, hf_score)
hist(top_LR$hf_score)
L= top_LR %>% pull(ligand)%>% unique()
R= top_LR %>% pull(receptor)%>% unique()


fibs.candidates=  gex.fib%>%
  rownames_to_column("gene")%>%
  filter(gene %in% R,
         pct.1>0.1,
         avg_log2FC>0)%>%
  pull(gene)

macs.candidates= gex.mac%>%
  rownames_to_column("gene")%>%
  filter(gene %in% L,
         pct.1>0.15,
         avg_log2FC>0.1)%>%
  pull(gene)


##

#  fibroblast -> mac activation in HFpEF
top_LR = liana_res$fib.mac%>%
  filter(hf_score >=0.9)%>%
  arrange(desc(hf_score))#%>%
  top_n(n= 40, hf_score)

L= top_LR %>% pull(ligand)%>% unique()
R= top_LR %>% pull(receptor)%>% unique()


rec.candidataes=  gex.mac%>%
  rownames_to_column("gene")%>%
  filter(gene %in% R,
         pct.1>0.1)%>%
  pull(gene)

lig.candidates= gex.fib%>%
  rownames_to_column("gene")%>%
  filter(gene %in% L,
         pct.1>0.1,
         avg_log2FC>0.1)%>%
  pull(gene)

cand= list("FM"= list("l"= lig.candidates,
                      "r"= rec.candidataes),
           "MF"= list("l"= macs.candidates,
                      "r" = fibs.candidates))
saveRDS(cand, file= "output/liana_candidates_MF.rds")
