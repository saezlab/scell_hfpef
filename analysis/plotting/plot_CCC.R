## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-03-16
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Plot sript for LIANA and Nichenet results 🐶

library(cowplot)
library(tidyverse)
library(Seurat)
library(igraph)

liana_res= readRDS("output/liana_combined_results.rds")
candidates= readRDS( file= "output/liana_candidates_MF.rds")
seu= readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")


#celltype seu:
seu.obj= readRDS("output/seu.objs/cell.state.obj.list.rds")
macs= seu.obj$macrophages
fibs= seu.obj$fibroblasts

state_marker= readRDS("output/cell_state_marker.rds")
fib_sigs=readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")

#

# aggregate rank overview ----------------------------------------------------------------------------------------


## plot overview
ps= map(liana_res, function(x){

  p.aggregate.rank= x %>%
    ggplot(., aes(x=aggregate_rank_hf, y= aggregate_rank_ct,  col = hf_score))+
    geom_point()+
    #geom_label_repel(aes(label= label), max.overlaps = 100)+
    labs(x= "aggr_ranks_hfpef",
         y= "aggr_ranks_ct",
         col= "HFpEF specificity")+
    theme_bw()

  p.hist= x %>%
    ggplot(., aes(x= hf_score))+
    geom_histogram()+
    theme_minimal_grid()

  list(p.aggregate.rank, p.hist)

})

p.comb= plot_grid(ps$fib.mac[[1]],
          ps$fib.mac[[2]],
          ps$mac.fibs[[1]],
          ps$mac.fibs[[2]],
          rel_widths =  c(1, 0.7),
          ncol= 2, labels= "AUTO")

pdf("output/figures/funcomics/liana_res_overview.pdf",
    height= 4)
p.comb
dev.off()

# add comparison of 3 groups ------------------------------------------------------------------

prio.df= liana_res$fib.mac %>%
  mutate(ct_prio= (1-aggregate_rank_ct)*aggregate_rank_hf)%>%
  mutate(both_prio= (1-aggregate_rank_ct)*(1-aggregate_rank_hf))%>% arrange(desc(both_prio))
prio.df%>%
  ggplot(., aes(x=aggregate_rank_hf, y= aggregate_rank_ct,  col = ct_prio))+
  geom_point()+
  #geom_label_repel(aes(label= label), max.overlaps = 100)+
  labs(x= "ranks_hfpef",
       y= "ranks_ct",
       col= "both specificity")+
  ggtitle("Liana LR pairs | macrophages->fibroblasts")+
  theme_minimal()
p.aggregate.rank

# add GEX data for candidate prioritization ---------------------------------------------------

## enrich Ligands in
# source("analysis/utils.R")
# GSE_analysis(L, Annotation_DB = fib_sigs$total)

# we focus on macrophage -> fibroblast activation in HFpEF

# plot L R genes to check for de.
plot_dots= function(L,
                    R,
                    seu,
                    sender,
                    receiver){
  Idents(seu) = "celltype"

  p.dot.ligands= DotPlot(subset(seu, idents = sender),
                         features = unique(L),
                         cols = "RdYlBu",
                         split.by = "group",
                         scale.min = 1,
                         scale.max =70,
                         ) +
    coord_flip()+
    theme(axis.text.x = element_text(angle= 60, hjust=1))+
    labs(y= "", x= "Ligands")+
    NoLegend()


  p.dot.receptor=
    DotPlot(subset(seu, idents = receiver),
            features = unique(R),
            cols = "RdYlBu",
            split.by = "group",
            scale.min = 1,
            scale.max = 70,
            scale.by = "radius") +
    coord_flip()+
    theme(axis.text.x = element_text(angle= 60, hjust=1))+
    labs(y= "", x= "Receptor")+
    NoLegend()

    cowplot::plot_grid(p.dot.ligands, p.dot.receptor)
}


#leg <- get_legend(fm)
fm= plot_dots(candidates$FM$l ,  candidates$FM$r, seu, sender= "Fibroblasts", receiver = "Macrophages")
mf= plot_dots(candidates$MF$l , candidates$MF$r, seu,receiver= "Fibroblasts", sender = "Macrophages")
dot.expression= plot_grid(fm , mf, ncol =1)

pdf("output/figures/funcomics/liana_dotplot_exprs_top25.pdf",
    heigh= 11,
    width= 6)
dot.expression
dev.off()

pdf("output/figures/supp/liana_dotplot_exprs_top25.pdf",
    heigh= 11,
    width= 6)
mf
dev.off()

Idents(seu) = "celltype"

gex.fib= FoldChange(subset(seu, idents = c("Fibroblasts")),
                    ident.1= "hfpef",
                    ident.2= "ct",
                    group.by= "group")

gex.mac= FoldChange(subset(seu, idents = c("Macrophages")),
                    ident.1= "hfpef",
                    ident.2= "ct",
                    group.by= "group")


# plot candidates ----------------------------------------------------------------------------
library(networkD3)
library(dplyr)

get_sankey= function(links){
  nodes= links%>% pivot_longer(c(ligand, receptor),  names_to = "type", values_to="name")%>%
    distinct(type, name)

  links$IDsource <- match(links$ligand, nodes$name)-1
  links$IDtarget <- match(links$receptor, nodes$name)-1

  # Make the Network

  my_color <- 'd3.scaleOrdinal() .domain(["ligand", "receptor"]) .range(["#69b3a2", "steelblue"])'
  p <- sankeyNetwork(Links = as.data.frame(links),
                     Nodes = as.data.frame(nodes),
                     Source = "IDsource",
                     Target = "IDtarget",
                     Value = "value",
                     NodeID = "name",
                     NodeGroup = "type",
                     sinksRight=FALSE)

  p

}

links= liana_res$mac.fibs %>%
  dplyr::filter(ligand %in% candidates$MF$l & receptor %in% candidates$MF$r)%>%
  left_join(gex.fib%>%
              rownames_to_column("receptor")%>%
              dplyr::rename(fc.fib= avg_log2FC), by= "receptor")%>%
  left_join(gex.mac%>%
              rownames_to_column("ligand")%>%
              dplyr::rename(fc.mac= avg_log2FC)%>%
              filter(ligand%in% candidates$MF$l), by= "ligand")%>%
  dplyr::select(ligand, receptor, hf_score, fc.mac, fc.fib)%>%
  mutate(label =  paste0(ligand, "->", receptor))

unique(links$ligand)
get_sankey(links %>% rename(value= hf_score))

#two macrophage ligands are also fibroblast receptors in our data:
double= links$ligand[links$ligand %in% links$receptor]
for (i in 1:length(double)){
  print(i)
  links$ligand=str_replace_all(links$ligand, double[i],  paste0(double[i], "_lig"))
  links$receptor= str_replace_all(links$receptor, double[i],  paste0(double[i], "_rec"))

}
links= links%>%
  mutate(label =  paste0(ligand, "->", receptor))

library(igraph)

## get Node
vertex= links %>% pivot_longer(c(ligand, receptor),  names_to = "type", values_to="node")%>%
  distinct(type, node)# %>% column_to_rownames("node")

vertex.df= vertex %>%#ownames_to_column(vertex, "name")  %>%
  left_join(links %>%
              distinct(ligand, fc.mac)%>%
              dplyr::rename(node= ligand, size= fc.mac), by= "node") %>%
  left_join(links %>%
              dplyr::distinct(receptor, fc.fib)%>%
                  dplyr::rename(node= receptor,
                         size2= fc.fib),
            by= "node") #%>%

vertex.df= vertex.df  %>%dplyr::mutate(size= ifelse(is.na(size), size2, size))%>%  dplyr::select(-size2)%>%
  dplyr::rename(name= node)


# save for plotting in cytoscape:
write.csv(vertex.df[,c(2,1,3)], "output/cytoscape_vertexdf.csv", row.names = F)
write.csv(links, "output/cytoscape_links.csv", row.names = F)

write_delim(vertex.df[,c(2,1,3)], "output/cytoscape_vertexdf.tsv")
write_delim(links, "output/cytoscape_links.tsv")


#plot with igrap
vertex.df
# create graph
#links= links%>%  filter(receptor!= "Adam15")
g= graph_from_data_frame(d = links, vertices = vertex.df[,c(2,3,1)])
V(g)$type= vertex[V(g)$name,1]
V(g)$type <- bipartite_mapping(g)$type
V(g)$size= V(g)$size * 20
V(g)$color <- ifelse(V(g)$type, "lightblue", "salmon")
E(g)$label = NA
V(g)$label = V(g)$name
plot(g, vertex.label.cex = 1, vertex.label.color = "black")
plot(g, layout=layout.fruchterman.reingold,
     vertex.size=22,
     vertex.label.cex=0.8)
?plot()
png("output/figures/funcomics/liana_res.net.png", 2000,2000)
plot(g,
     vertex.label.cex = 5,
     vertex.size= 15,
     vertex.label.color = "black",
     edge.width= 15,
     edge.arrow.size =5,
     edge.color= "black")
dev.off()






ct.sankey= comb.df2 %>%
  #filter(ligand %in% macs.candidates & receptor %in% fibs.candidates)%>%
  filter(ct_prio>0.95)%>%
  dplyr::select(ligand, receptor, ct_prio)%>%
  dplyr::rename(source= ligand,
                target= receptor,
                value= ct_prio)%>%
  get_sankey(.)
conserved.sankey= comb.df2 %>%
  #filter(ligand %in% macs.candidates & receptor %in% fibs.candidates)%>%
  filter(both_prio>0.99)%>%
  dplyr::select(ligand, receptor, both_prio)%>%
  dplyr::rename(source= ligand,
                target= receptor,
                value= both_prio)%>%
  get_sankey(.)
#
#     left_join(gex.fib%>%
#                 rownames_to_column("receptor")%>%
#                 rename(fc.fib= avg_log2FC), by= "receptor")%>%
#     left_join(gex.mac%>%
#                 rownames_to_column("ligand")%>%
#                 rename(fc.mac= avg_log2FC)%>%
#                 filter(ligand%in% macs.candidates), by= "ligand")%>%
#
#
#   nodes <- data.frame(
#     name=c(as.character(links$source),
#            as.character(links$target)),
#     expr=c(links$fc.mac, links$fc.fib)
#   )%>% distinct()
# intersect(links$source, links$target)
#   links$IDsource <- match(links$source, nodes$name)-1
#   links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network

my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#69b3a2", "steelblue"])'
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name",colourScale = my_color,
                   sinksRight=FALSE)

p


# t---------------------------------------------------------------------------------------------
p=   FeaturePlot(macs,cols= c("darkgrey","blue"),ncol = 2,
              features=c("Spp1", "Vegfa", "Gas6", "Pdgfa", "Tnf","Vegfb", "Gdf15",  "Plau"),
              combine= F, keep.scale = "all"
  )

p= map(p, function(x){
  x +theme(axis.title = element_text(size= 9),
           axis.text = element_blank(),
           axis.ticks= element_blank(),
           legend.position = "none",
           plot.title = element_text(size=10))


})
cowplot::plot_grid(plotlist = p, ncol = 2)

pls= FeaturePlot(macs,  cols= c("darkgrey","blue"),ncol = 2,
            features=c("Spp1", "Vegfa", "Tnf", "Gas6", "Pdgfa","Itga9", "Vegfb", "Gdf15",  "Plau"),combine = F, keep.scale = "all"  )


pls= map(pls, function(x){
 x+theme(axis.title = element_blank(),
         axis.text = element_blank())
})

plot_grid(plotlist = pls)

pdf("output/figures/main/Fig5/features.pdf",
    height= 6,
    width= 11)
plot_grid(plotlist = pls, ncol = 4)

dev.off()



Idents(macs) = "group"
p2= VlnPlot(macs,
#        idents = "cellstate",
            #cols= c("darkgrey","blue"),
            features=c("Spp1", "Vegfa", "Gas6", "Pdgfa", "Tnf","Vegfb", "Gdf15",  "Plau"))

pdf("output/figures/funcomics/ligand_macs_vln.pdf",
    height= 20,
    width= 10)
p2
dev.off()
# ripple heatmap -----------------------------------------------------------------------------


links= liana_res$mac.fibs %>%
  filter(ligand %in% macs.candidates & receptor %in% fibs.candidates)%>%
  left_join(gex.fib%>%
              rownames_to_column("receptor")%>%
              rename(fc.fib= avg_log2FC), by= "receptor")%>%
  left_join(gex.mac%>%
              rownames_to_column("ligand")%>%
              rename(fc.mac= avg_log2FC)%>%
              filter(ligand%in% macs.candidates), by= "ligand")%>%
  dplyr::select(ligand, receptor, hf_score, fc.mac, fc.fib)%>%
  mutate(label =  paste0(ligand, "->", receptor))

links %>% write.csv("output/liana_res.csv")
links= links%>% filter(hf_score>0.95)
h1= links %>% mutate(cell= "macrophages")%>%
  ggplot(., aes(x= cell, y= reorder(label,hf_score),  fill = fc.mac))+
  geom_tile()+
  scale_fill_gradient(low= "red", high = "darkred")+
  theme_minimal()+
  theme(#axis.text = element_blank(),
    axis.ticks = element_blank(),
    #legend.position = "none",
    axis.title = element_blank())
h1
h12= links%>% mutate(cell= "fibroblasts")%>%
  ggplot(., aes(x= cell, y= reorder(label,hf_score),  fill = fc.fib))+
  geom_tile()+
  scale_fill_gradient2(low= "blue", mid= "white", high = "red")+
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        #legend.position = "none",
        axis.title = element_blank())
h12
h2= links%>% mutate(cell= "macrophages")%>%
  ggplot(., aes(x= cell, y= reorder(label,hf_score), fill =hf_score))+
  geom_tile()+
  scale_fill_gradient(low= "white", high = "darkgreen")+
  theme_minimal()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        axis.title = element_blank())
h2

cowplot::plot_grid(h2, h1, h12, nrow = 1, rel_widths = c(0.7, 1,0.5))
# network plotting attempts -------------------------------------------------------------------
library(igraph)
links= hfpef.sankey
vertex= pivot_longer(links, names_to = "type", values_to="node", -value)%>%
  distinct(type, node)%>% column_to_rownames("node")


g= graph_from_data_frame(links)
V(g)$type= vertex[V(g)$name,1]
V(g)$type <- bipartite_mapping(g)$type
V(g)$color <- ifelse(V(g)$type, "lightblue", "salmon")

V(g)$size <- V(g)$expr
plot(g, vertex.label.cex = 0.5, vertex.size= 15,  vertex.label.color = "black")
plot(g, layout=layout.bipartite, vertex.size=30, vertex.label.cex=1)

png("output/figures/funcomics/liana_res.net.png", 2000,2000)
plot(g, vertex.label.cex = 5, vertex.size= 15,  vertex.label.color = "black",
     edge.width= 20, edge.arrow.size = 0)
dev.off()

library(visNetwork)
nodes <- data.frame(id = 1:3)
edges <- data.frame(from = c(1,2), to = c(1,3))
visNetwork(nodes = nodes%>% distinct() %>% dplyr::rename(id= name,
                                                         value= expr),
           edges =  links%>% dplyr::rename(to= source,
                                           from = target,
                                           strength= value), width = "100%")%>%
  visHierarchicalLayout()
##  to the netw
?visNetwork
<- graph.data.frame(links, directed=FALSE)
g <- make_bipartite_graph( rep(0:1,length=10), c(1:10))
print(g, v=TRUE)

fib_sigs= readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET.rds")

rs[rs %in% unlist(c(fib_sigs$unique$hfpef, fib_sigs$overlap))]

p.dot.receptor= DotPlot(subset(seu, idents = c( "fibroblasts")
                               #paste0("hfpef.",receiver)
), features = rs[1:20], cols = "RdYlBu", split.by = "group") +
  coord_flip()+
  theme(axis.text.x = element_text(angle= 60, hjust=1))

fibs= readRDS("output/cell_specific_obj/fibroblasts.rds")

FeaturePlot(fibs, features= rs[1:10])
pdf("output/figures/funcomics/liana_receptr_fib_expr.pdf" ,
    height= 25,
    width= 25)
VlnPlot(fibs, group.by = "opt_state" ,features = rs)
dev.off()



# add fib->mac --------------------------------------------------------------------------------

df1= liana_1 %>% filter(source=="fibroblasts",
                        target== "macrophages") %>%
  arrange((median_rank))%>%
  mutate(new_rank= rank(aggregate_rank, ties.method = "average"))
#print(n=100)

df2= liana_2 %>% filter(source=="fibroblasts",
                        target== "macrophages") %>%
  arrange((aggregate_rank))%>%
  mutate(new_rank= rank(aggregate_rank))

comb.df= df2 %>%
  inner_join(., df1, by= c("ligand", "receptor"))%>%
  mutate(new_score= (1-aggregate_rank.x)*aggregate_rank.y)%>%
  #mutate(new_score= (aggregate_rank.y)/(1-aggregate_rank.x))%>%
  arrange(desc(new_score))

# plot FIB-> Mac ------------------------------------------------------------------------------


p.aggregate.rank=
  comb.df%>%
  #filter(new_score>0.95)%>%
  ggplot(., aes(x=aggregate_rank.x, y= aggregate_rank.y,  col = new_score))+
  geom_point()+
  #geom_label_repel(aes(label= label), max.overlaps = 100)+
  labs(x= "ranks_hfpef",
       y= "ranks_ct",
       col= "hfpef specificity")+
  ggtitle("Liana LR pairs | fibroblasts->macrophages")+
  theme_minimal()
p.aggregate.rank

## plot head

p.head=
  comb.df%>%filter(aggregate_rank.x<0.1)%>%
  ggplot(., aes(x=aggregate_rank.x, y= aggregate_rank.y, col = new_score))+
  geom_point()+
  #geom_label_repel(aes(label= label), max.overlaps = 100)+
  labs(x= "ranks_hfpef",
       y= "ranks_ct",
       col= "hfpef specificity")+
  ggtitle("Liana LR pairs_ fibroblasts->macrophages")+
  theme_minimal()
p.head

##


## plot head with labels:
top_LR = comb.df%>% filter(new_score >0.95)

L= top_LR %>% pull(ligand)%>% unique()
R= top_LR %>% pull(receptor)%>% unique()


gex.fib= FoldChange(subset(seu, idents = c("fibroblasts")),
                    ident.1= "hfpef",
                    ident.2= "ct",
                    group.by= "group")

library(ComplexHeatmap)

fibs.candidates=  gex.fib%>%
  rownames_to_column("gene")%>%
  filter(gene %in% L,
         pct.1>0.1,
         avg_log2FC>0.1)%>%
  pull(gene)

gex.mac= FoldChange(subset(seu, idents = c("macrophages")),
                    ident.1= "hfpef",
                    ident.2= "ct",
                    group.by= "group")

macs.candidates= gex.mac%>%
  rownames_to_column("gene")%>%
  filter(gene %in% R,
         pct.1>0.1)%>%
  pull(gene)

gex.mac %>%
  rownames_to_column("gene")%>%
  filter(gene %in% macs.candidates)%>%
  dplyr::select(gene, avg_log2FC)%>%
  column_to_rownames("gene")%>%#filter(pct.1>0.1) %>%
  Heatmap()
Idents(seu) = "celltype"
p.dot.ligands= DotPlot(subset(seu, idents = c("macrophages")),
                       features = unique(macs.candidates),
                       cols = "RdYlBu",
                       split.by = "group") +
  coord_flip()+
  theme(axis.text.x = element_text(angle= 60, hjust=1))
p.dot.ligands
p.dot.receptor=
  DotPlot(subset(seu, idents = c( "fibroblasts")
                 #paste0("hfpef.",receiver)
  ), features = fibs.candidates, cols = "RdYlBu", split.by = "group") +
  coord_flip()+
  theme(axis.text.x = element_text(angle= 60, hjust=1))

p.dot.receptor

cowplot::plot_grid(p.dot.receptor,p.dot.ligands )



hfpef.sankey=
  comb.df %>%
  filter(ligand %in% fibs.candidates & receptor %in% macs.candidates)%>%
  #filter(new_score>0.95)%>%
  dplyr::select(ligand, receptor, new_score)%>%
  dplyr::rename(source= ligand,
                target= receptor,
                value= new_score)%>%
  get_sankey(.)


## only CT


comb.df2= comb.df %>%
  mutate(ct_prio= (1-aggregate_rank.y)*aggregate_rank.x)%>%
  mutate(both_prio= (1-aggregate_rank.y)*(1-aggregate_rank.x))%>% arrange(both_prio)
ct.sankey=
  comb.df2 %>%
  filter(ct_prio>0.999)%>%
  #filter(new_score>0.95)%>%
  dplyr::select(ligand, receptor, ct_prio)%>%
  dplyr::rename(source= ligand,
                target= receptor,
                value= ct_prio)%>%
  get_sankey(.)

ct.sankey
