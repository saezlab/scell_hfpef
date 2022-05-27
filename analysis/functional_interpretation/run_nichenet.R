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


library(parallelMap)
library(mlrMBO)
library(OmnipathR)
library(nichenetr)
library(tidyverse)
library(Seurat)
library(cowplot)
library(SingleCellExperiment)
library(biomaRt)

source("code/utils_nichenet.R")

source("code/utils.R")

source('~/R-projects/sc-exploration/R-scripts/utils.R')
# set up prior K ------------------------------------------------------------------------------

## 1. load networks and ligand_target matrix
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))

#if downloaded:
ligand_target_matrix = readRDS("data/prior_knowledge/nichenet/ligand_target_matrix.rds")
lr_network = readRDS("data/prior_knowledge/nichenet/lr_network.rds")
weighted_networks = readRDS("data/prior_knowledge/nichenet/weighted_networks.rds")

## 2. convert to mouse symobls

weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

lr_network = lr_network %>%
  mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
convert_human_to_mouse_symbols(colnames(ligand_target_matrix))
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

## 3. define receiver and sender cells and their ligand - target gene sets
# nichenet suggests to use expression cut-offs. I will use the DE-genes from above to in combination with expression cut-off


# load data -----------------------------------------------------------------------------------

seu= readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")

Idents(seu)= "celltype"
unique(Idents(seu))
#fib_sigs= readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
# get the fibroblast signature( up and dn)
hfpef_f= readRDS( "output/fib_integration/marker_list/DEG_hfpef.rds")

#get the candidates from LIANA
candidates = readRDS("output/liana_candidates_MF.rds")

# prepare nichnet input
input= prep_nichenet_input(sender = "Macrophages",
                    receiver = "Fibroblasts",
                    geneset_sender = candidates$MF$l,
                    geneset_receiver= hfpef_f$gene)
# MAIN
# predict ligand activity
ligand_activities = predict_ligand_activities(geneset = input$geneset_oi,
                                              background_expressed_genes = input$background_expressed_genes,
                                              ligand_target_matrix = ligand_target_matrix,
                                              potential_ligands = input$potential_ligands)

ligand_activities = ligand_activities %>%
  arrange(-pearson) %>%
  mutate(rank = rank((-pearson)))
p.ligand.activities= ligand_activities %>%
  mutate(cell= "macrophages")%>%
  ggplot(aes(x= cell, y= reorder(test_ligand,pearson), fill = pearson))+
  geom_tile()+
  scale_fill_gradient(low= "#2596be", high = "darkblue")+
  theme_minimal()+
  theme(axis.text = element_text(color= "black",size= 11),
        axis.text.x= element_blank(),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9),
        legend.position = "left")+
  labs(x= "",y= "Macrophage ligands")
pdf("output/figures/main/Fig5/nichenet_ligand.activities.pdf",
    width = 2,
    height= 3)
p.ligand.activities
dev.off()

ligand_activities
## plot - histogram of pearson corr for all ligands
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) +
  geom_histogram(color="black", fill="darkorange")  +
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(10, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) +
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity

# this plot can help to select the best top ligands. Here we only have a few ligands we tested,

best_upstream_ligands = ligand_activities %>%
  #top_n(10, pearson) %>%
  arrange(-pearson) %>%
  pull(test_ligand) %>%
  unique()

# Data transformation for heatmap plotting
active_ligand_target_links_df = best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links, geneset = input$geneset_oi, ligand_target_matrix = ligand_target_matrix) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                 ligand_target_matrix = ligand_target_matrix,
                                                                 cutoff = 0.0001)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

#  heatmap plotting

p_ligand_target_network = vis_ligand_target %>%
  make_heatmap_ggplot("","Fibroblast target genes",
                      color = "purple",
                      legend_position = "right",
                      x_axis_position = "bottom",
                                                                    legend_title = "Regulatory potential")  +
  theme(axis.text.x = element_text(face = "italic", color ="black"),
        axis.text= element_text(color= "black"))
  #scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
print(p_ligand_target_network)


pdf("output/figures/main/Fig5//nichenet_res.pdf",
    width= 5.5,
    height= 3.5)
p_ligand_target_network
dev.off()

p.comb= cowplot::plot_grid(p.ligand.activities, p_ligand_target_network, align= "h", axis= "tb", rel_widths = c(1,2.5))

pdf("output/figures/main/Fig5//niche_net_comb.pdf",
    width= 7.5,
    height= 3)
p.comb
dev.off()

  # omnipath with nichenet ----------------------------------------------------------------------

library(OmnipathR)
library(mlrMBO)
library(mlr)
library(nichenetr)
library(tidyverse)
library(Seurat)
library(magrittr)
library(parallelMap)
library(ParamHelpers)
library(smoof)
receiver= "fibroblasts"
sender= "macrophages"

expressed_genes_receiver = get_expressed_genes(receiver, seu, pct = 0.10, assay_oi = "RNA")
background_expressed_genes = expressed_genes_receiver %>%
  .[. %in% rownames(ligand_target_matrix)]

# sender:
list_expressed_genes_sender = sender  %>% unique() %>% lapply(get_expressed_genes, seu, 0.1, assay_oi= "RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

genes_oi= geneset_receiver = c(fib_sigs$unique$hfpef, unlist(fib_sigs$overlap))

t= nichenet_main(expressed_genes_receiver = expressed_genes_receiver,
              expressed_genes_transmitter = expressed_genes_sender,
              genes_of_interest = genes_oi
               )
nichenet_workarounds()
nichenet_test()
