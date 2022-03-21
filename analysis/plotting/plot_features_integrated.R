## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-02-04
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## plot integrated Feature plots


library(Seurat)
library(tidyverse)

source("analysis/utils.R")
source("analysis/utils_funcomics.R")
int.fibs= readRDS("output/seu.objs/study_integrations/harmony_fib_filt.rds")

gene_signatures= readRDS( "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
state_marker= readRDS("output/fib_integration/marker_list/integrated_marker.rds")

#prepare_cell state marker:
marker = state_marker %>% group_by(cluster) %>% filter(p_val_adj<0.05,
                                                       avg_log2FC>0)%>%
  top_n(n = 150, avg_log2FC)
marker=split(marker$gene, marker$cluster)

# use featureplots ------------------------------------------------------------------------
nx= map(marker , function(x){FeaturePlot(int.fibs, features = x[1:6])})

pdf("output/figures/integration_studies/marker_feature_plot.pdf",
    width= 7,
    height= 7)
nx
dev.off()

FeaturePlot(int.fibs, features = c("Acta2", "Postn", "Angptl4", "Thbs4", "Cthrc1"))
FeaturePlot(int.fibs, features = c("Igfbp3" , "Cxcl14", "Pi16", "Cilp"))

markx= map(marker , function(x){ x[1:6]})

p.marker= DotPlot(int.fibs, group.by = "opt_clust_integrated", features = unique(unlist(markx)), cols= up_dn_cols,
                  cluster.idents = F)+
  coord_flip()+
  theme(axis.text.x= element_text(angle= 40, hjust = 1))+
  labs(x= "", y= "")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

p.marker
colnames(int.fibs[[]])
