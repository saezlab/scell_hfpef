## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-10-25
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## 1.) Identify markers of the integrated fibroatlas
## 2.) Perform ORA with interpretational genesets

library(Seurat)
library(tidyverse)
library(WriteXLS)

integrated_data = readRDS( "output/seu.objs/study_integrations/harmony_fib_filt.rds")
source("analysis/utils.R")

# find marker  --------------------------------------------------------------------------------
Idents(integrated_data)= "opt_clust_integrated"
DefaultAssay(integrated_data)= "RNA"
# Get markers per cluster:
marker_int= FindAllMarkers(integrated_data, assay= "RNA", slot = "data")

#save rds
saveRDS(marker_int, "output/fib_integration/marker_list/integrated_marker.rds")

#save as xlss
markers.list= split(x=marker_int, f = marker_int$cluster)
WriteXLS(markers.list,
           ExcelFileName = paste0("output/fib_integration/marker_list/statemarker.xlsx")
  )

int.fib.marker= readRDS( "output/fib_integration/marker_list/integrated_marker.rds")
marker.fibs= readRDS( "output/fib_integration/marker_list/local_fib_state_marker.rds")

cut.off.state= 75
marker.fibs.genes= marker.fibs%>% group_by(cluster) %>%
  filter(p_val_adj <0.05,
         avg_log2FC>0)%>%
  slice(1:cut.off.state)
local.genesets= split(marker.fibs.genes$gene, marker.fibs.genes$cluster)

int.fib.marker= int.fib.marker%>% group_by(cluster) %>%
  filter(p_val_adj <0.05,
         avg_log2FC>0)%>%
  slice(1:cut.off.state)
int.genesets= split(int.fib.marker$gene, int.fib.marker$cluster)


## 1.NABA ECM

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

# use cell state marker of fibroblasts to track myofibs in fib atlas --------------------------
add_module_for_genesets= function(genesets, seurat_obj, col.name){

  sets= names(genesets)
  gset= sets[1]
  for(gset in sets){
    features = list(gset = genesets[[gset]])
    ctrl_genes = sum(genesets[[gset]] %in% rownames(seurat_obj))
    #ctrl_genes = 35

    if(ctrl_genes>0){
      seurat_obj = AddModuleScore(object = seurat_obj,
                                  features = features,
                                  name = paste0(col.name,gset),
                                  ctrl = ctrl_genes,
                                  seed = 77)
    }
  }
  return(seurat_obj)
}

integrated_data= add_module_for_genesets(genesets, integrated_data, col.name= "hfpef.fib.clust")

##plot feature plots:
pls= map(paste0("hfpef.fib.clust", names(genesets), "1"), function(x){
  VlnPlot(integrated_data, x)
}
)
pls
integrated_data$hf
meta_marker = read_csv("~/R-projects/sc-exploration/output/marker/meta_marker.csv")

top_myofib_genes= meta_marker%>% top_n(., n= 50, wt= meanlogFC) %>% pull(gene)

gset_myo= list("myofib"= top_myofib_genes)
integrated_data= add_module_for_genesets(gset_myo, integrated_data, col.name = "")
integrated_data$myofib1
VlnPlot(integrated_data, "myofib1")
FeaturePlot(integrated_data, "myofibmyofib1")
DimPlot(integrated_data)

## Myofibroblasts are clust 2!

x= integrated_data@meta.data
saveRDS(x, "output/figures/integration_studies/meta.data.clustered.rds")
x%>% filter(opt_clust_integrated==2)


# ADD MODULE SCORES to aid with  cluster interpret -----------------------------------------------------
# we are going to add different gene set scores to help with cluster interpretation
#1) NABA
#2) DE genes from nature sc fibro atlas paper
#3) ReHeat disease score. (pending)

## 1) NABA

integrated_data= add_geneset_scores(integrated_data, NABA_SETS_mouse)

# Plot NABA Score overview:
meta= integrated_data@meta.data

meta= meta %>%
  gather(colnames(ECMscores), key= "NABA", value= "ECM_score") %>%
  rename(cluster= opt_clust_integrated)
p.naba.overview=  ggplot(data= meta, aes(x= cluster, y= ECM_score,fill = cluster))+
  geom_boxplot()+
  theme_minimal()+
  #scale_fill_manual(values = col.set)+
  #scale_color_manual(values = col.set)+
  facet_grid(cols= vars(NABA))

p.naba.select=  ggplot(data= meta %>% filter(!NABA %in% c("NABA_ECM_AFFILIATED","NABA_SECRETED_FACTORS", "NABA_MATRISOME_ASSOCIATED")),
                       aes(x= cluster, y= ECM_score,fill = cluster))+
  geom_boxplot()+
  theme_minimal()+
  #scale_fill_manual(values = col.set)+
  #scale_color_manual(values = col.set)+
  facet_grid(cols= vars(NABA))

p.naba.feat= FeaturePlot(integrated_data, features= colnames(ECMscores), keep.scale = "all")
p.naba.feat.noscale= FeaturePlot(integrated_data, features= colnames(ECMscores))


## 2) DE genes from nature sc fibro atlas paper
lst2= readRDS("data/prior_knowledge/fibro_atlas_deg/fibro_atlas_steady_state_marker.rds")

de.genes= map(lst2, function(x){
  x%>% arrange(p_val_adj)%>% slice(1:50)  %>% pull(Gene)
})

# add module scores and plot:
integrated_data= add_geneset_scores(integrated_data, de.genes)

meta= integrated_data@meta.data

meta= meta %>%
  gather(paste0(names(de.genes),"1"), key= "fib.cluster.steady.state", value= "fib_score") %>%
  rename(cluster= opt_clust_integrated)

p.fib.overview=  ggplot(data= meta, aes(x= cluster, y= fib_score,fill = cluster))+
  geom_boxplot()+
  theme_minimal()+
  #scale_fill_manual(values = col.set)+
  #scale_color_manual(values = col.set)+
  facet_grid(cols= vars(fib.cluster.steady.state))

p.fib.feature= FeaturePlot(integrated_data, features = paste0(names(de.genes),"1"), keep.scale= "all")
p.fib.feature.nocscale= FeaturePlot(integrated_data, features = paste0(names(de.genes),"1"))
pdf("output/fib_integration/fibroblast_module_scores.pdf",
    width= 15,
    height= 8)
p.naba.overview
p.naba.select
p.naba.feat
p.naba.feat.noscale
p.fib.overview
p.fib.feature
p.fib.feature.nocscale
dev.off()


# do ORA in addition to module score ----------------------------------------------------------
plot_ORA=function(Ora_res ){

  for (i in names(Ora_res)){
    Ora_res[[i]] = Ora_res[[i]] %>% mutate(cluster= i)
  }

  Ora_res2= as_tibble(do.call(rbind, Ora_res)) %>%
    #rename(cluster= gset) %>%
    mutate(corr_p_value = p.adjust(p_value, method= "BH"))%>%
    mutate(stars= ifelse(corr_p_value<0.05, "*", ""),
           stars= ifelse(corr_p_value<0.01, "**", stars),
           stars= ifelse(corr_p_value<0.001, "***",stars))


  p.p.ora= ggplot(Ora_res2, aes(x=  cluster,
                                y= gset, fill= -log10(corr_p_value)))+
    geom_tile()+
    scale_fill_gradient2(low= "white" , high= "red")+
    geom_text(mapping = aes(label= stars))+
    theme_minimal()+
    labs(fill= "-log10(q-value)")+
    theme(axis.text.x= element_text(angle=40, hjust= 1, size= 10),
          axis.title = element_blank(),
          axis.text= element_text(color ="black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )#+coord_flip()
  p.p.ora
  return(p.p.ora)
}

#NABA
Ora_res= lapply(int.genesets, function(x){
  GSE_analysis(x, NABA_select)

})

p.naba= plot_ORA(Ora_res)

pdf("output/figures/integration_studies/fibroblast_naba_ORA.pdf",
    width = 4.5,
    height= 2.5)

p.naba
dev.off()
## 2) DE genes from nature sc fibro atlas paper
steady.state= readRDS("data/prior_knowledge/fibro_atlas_deg/fibro_atlas_steady_state_marker.rds")
pert.state= readRDS("data/prior_knowledge/fibro_atlas_deg/fibro_atlas_perturbed_state_marker.rds")

steady.genes= map(steady.state, function(x){
  x%>% arrange(p_val_adj)%>%
    filter(avg_logFC>0,
           p_val_adj<0.05)%>%
    slice(1:75)  %>%
    pull(Gene)
})

pert.state= map(pert.state, function(x){
  x%>% arrange(p_val_adj)%>%
    filter(avg_logFC>0,
           p_val_adj<0.05)%>%
    slice(1:75)  %>%
    pull(Gene)
})

Ora_res.steady= lapply(int.genesets, function(x){
  GSE_analysis(x , steady.genes)

})

Ora_res.pert= lapply(int.genesets, function(x){
  GSE_analysis(x, pert.state)

})

p.pert= plot_ORA(Ora_res.pert)
p.stead= plot_ORA(Ora_res.steady)

pdf("output/figures/integration_studies/fibroblast_integration_state_ORA.pdf",
    width = 5,
    height= 5)
cowplot::plot_grid(p.pert+ggtitle("perturbed state marker"),
                   p.stead+ggtitle("steady state marker"),
                   #labels = "AUTO",
                   ncol=1)
dev.off()

##4 marker compar

Ora_res= lapply(int.genesets, function(x){
  GSE_analysis(x, local.genesets)

})


plot_ORA(Ora_res)

