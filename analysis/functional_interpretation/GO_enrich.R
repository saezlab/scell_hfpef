## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-19
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## perform GO enrichment of different gene sets
## 1) Cellstate marker (Fibs and Macs)
## 2) DEG fibs

library(Seurat)
library(tidyverse)
library(WriteXLS)
library(ggpubr)
library("enrichR")

source("analysis/utils.R")

# prepare GO ----------------------------------------------------------------------------------

## prepare DBs

dbs <- listEnrichrDbs()

GO_dbs= dbs$libraryName[grepl("GO_", x = dbs$libraryName)]
GO_dbs= GO_dbs[grepl("2021", x = GO_dbs)]

## function to apply enrichr. for GO terms to a list of gene sets

enrich.wrap= function(sets){

  enrich.r= lapply(sets, function(x){
    y= enrichr(genes = x,
               databases =GO_dbs)

  })


  bio.process= map(names(enrich.r), function(x){
    enrich.r[[x]]$GO_Biological_Process_2021%>% mutate(cluster= x)
  })%>% do.call(rbind, .)%>% as_tibble()

  mol.func= map(names(enrich.r), function(x){
    enrich.r[[x]]$GO_Molecular_Function_2021%>% mutate(cluster= x)
  })%>% do.call(rbind, .)%>% as_tibble()

  return(list("mf"= mol.func,
              "bp"= bio.process))
}

## clean names of gene sets (remove GO ID)
clean_names= function(df, col= "Term"){
  df[[col]]=  gsub(pattern = "\\(.*",replacement = "",x = df[[col]])
  df
}

filter_terms= function(df, terms= c("matrix","inflam","fibro", "immune", "heat"), pos=T){

  hits= lapply(terms, function(x){
    df$Term[grepl(x, df$Term)]
  })%>% unlist()

  if(pos){
    df = df %>% filter(Term %in% hits)
  }else{
    df = df %>% filter(!Term %in% hits)
  }
  }
# cell state marker ---------------------------------------------------------------------------

##int state fib marker
int.fib.marker= readRDS( "output/fib_integration/marker_list/integrated_marker.rds")

cut.off.state= 100

int.fib.marker= int.fib.marker%>% group_by(cluster) %>%
  filter(p_val_adj <0.05,
         avg_log2FC>0)%>%
  slice(1:cut.off.state)

int.genesets= split(int.fib.marker$gene, int.fib.marker$cluster)

## prepare gene sets:

cellstate_marker =readRDS("output/cell_state_marker.rds")

mac.marker= cellstate_marker$macrophages%>%
  group_by(cluster) %>%
  filter(p_val_adj<0.05,
         avg_log2FC>0)%>%
 top_n(100, p_val_adj)
# %>%
#   mutate(cluster= as.character(cluster),
#     cluster= ifelse(cluster==0, "Ccr2-/MHCII-/resident", cluster),
#          cluster= ifelse(cluster==1, "Ccr2-/MHCII-/Lyve1+", cluster),
#          cluster= ifelse(cluster==2, "Ccr2+/Ly6c+/monocyte", cluster),
#          cluster= ifelse(cluster==3, "Ccr2-/MHCII-/Cxcl2+", cluster),
#          cluster= ifelse(cluster==4, "Ccr2+/MHCII+", cluster))

macs.sets= split(x = mac.marker$gene, f= mac.marker$cluster)

fib.marker= cellstate_marker$fibroblasts%>%
  group_by(cluster) %>%
  filter(p_val_adj<0.05,
         avg_log2FC>0)#%>%
  top_n(100, p_val_adj)

fibs.sets= split(x = fib.marker$gene, f= fib.marker$cluster)


## perform enrich
macs.r= enrich.wrap(macs.sets)
fibs.r= enrich.wrap(fibs.sets)

int.r = enrich.wrap(int.genesets)

df=int.r$bp%>%
  filter(Adjusted.P.value<0.1)# %>%

split(df, f= df$cluster)%>%
  WriteXLS(., "output/GO_fib_cellstate_int.xls")

saveRDS(int.r,
        "output/GO_cellstates_int.rds")
#save as xls
df=fibs.r$bp%>%
  filter(Adjusted.P.value<0.05)# %>%
distinct(cluster, Term)

library(WriteXLS)
split(df, f= df$cluster)%>%
  WriteXLS(., "output/supp_tables/GO_fib_cellstate.xls")

df=macs.r$bp%>%
  filter(Adjusted.P.value<0.05) %>% filter(cluster==2)
split(df, f= df$cluster)%>%
  WriteXLS(., "output/supp_tables/GO_mac_cellstate.xls")

#save obj

saveRDS(list("fibs"= fibs.r, "macs"= macs.r),
          "output/funcomics_res/GO_cellstates.rds")



# DEGs of fibroblasts form differnt studies ----------------------------------------------------------------------------------------
signatures=readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
map(signatures$unique, length)

## perform enrich
fib.sigs=  enrich.wrap(signatures$total)
fib.sigs= lapply(fib.sigs, clean_names)


## BP
T_count= fib.sigs$bp%>%
  filter(Adjusted.P.value<0.05)%>%
  select(Term, cluster)%>% group_by(Term) %>% count()

common_T= T_count%>%
  filter(n==3)%>% pull(Term)

fib.sigs$bp = fib.sigs$bp %>% left_join(T_count)

p.bp.common= fib.sigs$bp%>%
  filter(Term %in% common_T)%>%
  ggplot(aes(x= cluster, y= reorder(Term,-log10(Adjusted.P.value)),  fill= -log10(Adjusted.P.value)))+
  geom_tile()+
  scale_fill_gradient(low = "darkred", high = "red")+
  labs(x="", y= "")+
  theme_minimal()

#pring HFpEF top
p.sigs= map(unique(fib.sigs$mf$cluster), function(x){
  df = fib.sigs$bp
  df= filter_terms(df, terms=c("eye","nigra", "COPI", "sulfur", "odonto"), pos = F)
  df%>%
    filter(!Term %in% common_T,
           cluster ==x,
           Adjusted.P.value<0.05)%>%
    mutate(n= as.factor(n))%>%
    group_by(cluster)%>%
    top_n(10, -Adjusted.P.value)%>%
    ggplot(aes(x= reorder(Term,-log10(Adjusted.P.value)),
               y=  -log10(Adjusted.P.value)
               #               fill = n
               )
           )+
    geom_col()+
    #facet_grid(vars(rows=cluster))+
    coord_flip()+
    labs(y="-log10(p.adj)", x= "")+
    theme_minimal()


})

p.bp= cowplot::plot_grid(p.sigs[[1]] + labs(y=""),
                         p.sigs[[2]] + labs(y=""),
                         p.sigs[[3]],
                         ncol = 1, align = "v")
p.bp

###MF

T_count= fib.sigs$mf%>%
  filter(Adjusted.P.value<0.05)%>%
  select(Term, cluster)%>% group_by(Term) %>% count()

common_T= T_count%>%
  filter(n==3)%>% pull(Term)

fib.sigs$mf = fib.sigs$mf %>% left_join(T_count)

p.mf.common= fib.sigs$mf%>%
  filter(Term %in% common_T)%>%
  ggplot(aes(x= cluster, y= reorder(Term,-log10(Adjusted.P.value)),  fill= -log10(Adjusted.P.value)))+
  geom_tile()+
  scale_fill_gradient(low = "darkred", high = "red")+
  labs(x="", y= "")+
  theme_minimal()

single_T= fib.sigs$mf%>%
  filter(Adjusted.P.value<0.05)%>%
  select(Term, cluster)%>% group_by(Term) %>% count()%>%
  filter(n==1)%>% pull(Term)

#pring HFpEF top
p.sigs= map(unique(fib.sigs$mf$cluster), function(x){
  df = fib.sigs$mf
  df= filter_terms(df, terms=c("eye","nigra", "COPI", "sulfur", "odonto"), pos = F)
  df%>%
    filter(!Term %in% common_T,
           cluster ==x,
           Adjusted.P.value<0.05)%>%
    mutate(n= as.factor(n))%>%
    top_n(15, -Adjusted.P.value)%>%
    ggplot(aes(x= reorder(Term,-log10(Adjusted.P.value)),
               y=  -log10(Adjusted.P.value)
               #               fill = n
    )
    )+
    geom_col()+
    coord_flip()+
    labs(y="-log10(p.adj)", x= "")+
    theme_minimal()+
    ggtitle(x)

})

p.mf= cowplot::plot_grid(plotlist = p.sigs, ncol = 1, align = "v")


p.mf
p.bp
p.bp.common
p.mf.common
pdf("output/figures/integration_studies/deg.hmaps/GO_common.bp",
    width = )

  # plot heatmap: -------------------------------------------------------------------------------

table(df$Term, df$cluster)

bio.process= filter_terms(int.r$bp, terms)

bio.process%>%
  #filter(Adjusted.P.value<0.01) %>%
  group_by(cluster)%>% top_n(., n = 3, wt = -Adjusted.P.value) %>%
  ggplot(., aes(y= Term, x= cluster, fill = -log10(Adjusted.P.value)))+
  geom_tile()+
  theme(axis.text.x = element_text(angle= 60, hjust = 1))

mol.func%>%
  #filter(Adjusted.P.value<0.01) %>%
  group_by(cluster)%>% top_n(., n = 3, wt = -Adjusted.P.value) %>%
  ggplot(., aes(y= Term, x= cluster, fill = -log10(Adjusted.P.value)))+
  geom_tile()

bio.process%>%filter(cluster== 0)

bio.process %>% ggplot(., aes(x= cluster, y= Term, fill = -log10(Adjusted.P.value)))+
  geom_tile()

# add sankey ----------------------------------------------------------------------------------

library(networkD3)
library(dplyr)
plot_sankey= function(df, p.val.alpha){
  links= df %>%
    group_by(cluster)%>%
    dplyr::filter(Adjusted.P.value<p.val.alpha) %>%
    #top_n(., n= 5, wt = Adjusted.P.value)%>%
    mutate(n= -log10(Adjusted.P.value))%>%
    dplyr::select(Term, cluster, n)%>%
    dplyr::rename(source= Term,
                  target= cluster,
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

}

# enrich the fibrosis signatures with GO terms---------------------------------------------------------------------------------------------

signature= readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled2.rds")

enrich.r= lapply(signature$total, function(x){
  y= enrichr(genes = x,
             databases =GO_dbs)

})

bio.process= map(names(enrich.r), function(x){
  enrich.r[[x]]$GO_Biological_Process_2018%>% mutate(cluster= x)
})%>% do.call(rbind, .)%>% as_tibble()


terms= c("matrix","inflam","fibro", "immune", "heat", "heart", "cardia", "cyto")
bio.process= filter_terms(bio.process, terms)

mol.func= map(names(enrich.r), function(x){
  enrich.r[[x]]$GO_Molecular_Function_2018%>% mutate(cluster= x)
})%>% do.call(rbind, .)%>% as_tibble()

cell.comp= map(names(enrich.r), function(x){
  enrich.r[[x]]$GO_Cellular_Component_2018%>% mutate(cluster= x)
})%>% do.call(rbind, .)%>% as_tibble()

plot_sankey(bio.process, 0.1)
plot_sankey(mol.func, 0.05)
plot_sankey(cell.comp, 0.05)
signatures.all = c(signature$unique, signature$overlap)


map(names(signature$unique ), function(x){
    enframe(signature$unique[[x]], name = "geneset", value = "gene")%>%
      select(gene, geneset)%>%
      mutate(geneset =x)
    })%>% do.call(rbind,.)%>%
    mutate(gene= toupper(gene)) %>%
    write.csv(., file= "output/fib_integration/marker_list/marker_for_ReHeaT_enrichment.csv", row.names = F)
??toupper
