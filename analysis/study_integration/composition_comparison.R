## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-09
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Proportion analysis of integrated cell state cluster with single study cell state cluster


library(Seurat)
library(tidyverse)

source("analysis/utils.R")

int.fibs= readRDS("output/seu.objs/study_integrations/harmony_fib_filt.rds")
meta.data= int.fibs@meta.data
rm(int.fibs)

# load sc data
seu.objs=  readRDS("output/seu.objs/cell.state.obj.list.rds")
fibs= seu.objs$fibroblasts
rm(seu.objs)

#function that takes a meta file and creates a new cellID as universal mapping :
add_universal_cellid= function(seu.meta, vec= 1){

  #get clean barcode:
  new_cellid= str_split(rownames(seu.meta), pattern = "_")
  new_cellid= map(new_cellid, function(x){x[[vec]]})

  #add cleaned barcode with orig ident
  new_cellid= paste0(unlist(new_cellid), seu.meta$orig.ident)
  seu.meta$cellid_new= new_cellid
  return(seu.meta)
}


meta.data= add_universal_cellid(meta.data,vec= 2)

fibs.meta= add_universal_cellid(fibs@meta.data, 1)


proportion_of_cell_states= function(meta.study, study= "hfpef" ){

  #get meta.data from target.seu
  #meta.study = seu[[]]
  #meta.study= add_universal_cellid(meta.study)


  meta.data.2= meta.data %>%
    filter(study== "hfpef")%>%
    select(cellid_new, opt_clust_integrated)

  x= inner_join(meta.study %>% select(-opt_clust_integrated), meta.data.2, by= "cellid_new")


}
joined.meta= proportion_of_cell_states(fibs.meta)
y= calc_props(joined.meta, cluster.col ="opt_clust_integrated",  group.col = "cellstate")

y2= calc_props(joined.meta, group.col ="opt_clust_integrated",  cluster.col  = "cellstate")

pdf(file= "output/figures/cell_type_analysis/fibros/fib_integrated_cluster_mapping.pdf",
    height= 4, width=4)
y+labs(fill= "hfpef_cellstate",
       y= "proportion (%)",
        x= "integrated_cellstate")
y2+labs(fill= "integrated_cellstate",
        y= "proportion (%)",
        x= "hfpef_cellstate")
dev.off()



# add sankey ----------------------------------------------------------------------------------
library(networkD3)
library(dplyr)

links= joined.meta %>% select(cellid_new, opt_clust_integrated, cellstate)%>%
  group_by(cellstate, opt_clust_integrated)%>%
  count()%>% rename(source= cellstate,
                    target= opt_clust_integrated,
                    value= n)

nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>%
    unique()
)

links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1


# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name",
                   sinksRight=FALSE)

p
