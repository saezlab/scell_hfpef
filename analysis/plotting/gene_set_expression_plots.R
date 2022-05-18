## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-03-30
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## plot ecm related genes sets
library(Seurat)
library(tidyverse)

source("code/utils.R")
fibs= readRDS("output/seu.objs/study_integrations/harmony_fib_filt.rds")
fibs= subset(fibs, group =="hf")


NABA_SETS= processNABA()

#use mapping dl earlier from biomart
gene_translate= readRDS("~/R-projects/sc-exploration/output/gene_translate.rds")

# translate gene sets to mouse
NABA_SETS_mouse= lapply(NABA_SETS, function(set){
  #print(set)
  tibble(set)  %>%
    dplyr::rename(Gene.name = set) %>%
    left_join(gene_translate, by= "Gene.name") %>%
    drop_na()%>%
    dplyr::pull(MGI.symbol)%>% unique()

})

names(NABA_SETS_mouse)= str_replace_all(names(NABA_SETS_mouse),"NABA_", "")
names(NABA_SETS_mouse)= str_to_title(names(NABA_SETS_mouse))
NABA_SETS_mouse

