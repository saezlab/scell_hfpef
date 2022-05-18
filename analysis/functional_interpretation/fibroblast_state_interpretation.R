## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-10-28
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
## Characterization of fibroblasts cell states.
## Use module scores
## Perform Overrepresentation analysis with different gene sets.

library(Seurat)
library(tidyverse)
library(WriteXLS)
library(rstatix)
library(ggpubr)

source("code/utils.R")


# load data and genesets ----------------------------------------------------------------------

### load sc data
seu.objs=  readRDS("output/seu.objs/cell.state.obj.list.rds")

seu.objs$fibroblasts@meta.data = seu.objs$fibroblasts@meta.data%>%
  mutate(cellstate= factor(cellstate, levels =  c("Col15a1+", "Igfbp3+",
                                                  "Pi16+","Cxcl1+","Cilp+","Wif1+")))

seu= seu.objs$fibroblasts
rm(seu.objs)



### Genesets

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
saveRDS(NABA_SETS_mouse, "data/prior_knowledge/NABA_mouse.rds")
NABA_SETS_mouse= readRDS( "data/prior_knowledge/NABA_mouse.rds")

## 2. Fibroblast atlas cell state marker (Buechler et. al 2021)

stead.marker= readRDS("data/prior_knowledge/fibro_atlas_deg/fibro_atlas_steady_state_marker.rds")
pert.marker= readRDS("data/prior_knowledge/fibro_atlas_deg/fibro_atlas_perturbed_state_marker.rds")

names(pert.marker) = paste0("pert.state", names(pert.marker))
names(pert.marker) = str_replace_all(names(pert.marker), pattern = "AllMarkers_", ".")
names(stead.marker) = paste0("steady.state", names(stead.marker))
names(stead.marker) = str_replace_all(names(stead.marker), pattern = "AllMarkers_", ".")

# define number of genes per state
cut.off.intmarker= 75

stead.genes= map(stead.marker, function(x){
  x%>% filter(p_val_adj<0.05,avg_logFC>0) %>% arrange(p_val_adj)%>% slice(1:cut.off.intmarker)  %>% pull(Gene)
})

pert.genes= map(pert.marker, function(x){
  x%>% filter(p_val_adj<0.05,avg_logFC>0) %>% arrange(p_val_adj)%>% slice(1:cut.off.intmarker)  %>% pull(Gene)
})


#save list with multiple cutoffs for other projects:
list.state.m= map(c(50,75,100,150), function(cut.off.intmarker){
  stead.genes= map(stead.marker, function(x){
    x%>% filter(p_val_adj<0.05,avg_logFC>0) %>% arrange(p_val_adj)%>% slice(1:cut.off.intmarker)  %>% pull(Gene)
  })

  pert.genes= map(pert.marker, function(x){
    x%>% filter(p_val_adj<0.05,avg_logFC>0) %>% arrange(p_val_adj)%>% slice(1:cut.off.intmarker)  %>% pull(Gene)
  })
  list("steady"= stead.genes,
       "pert"= pert.genes)

})
names(list.state.m)= paste0("top", c(50,75,100,150))
saveRDS(list.state.m,"data/prior_knowledge/fibro_atlas_deg/cross_organ_fibro_marker.rds")

target1= data.frame("states"= names(list.state.m$top50$steady),
                   "organ"= c("Pan", "Pan", "Spleen/Lymph", "Bone", "Lung", "Artery/Tendon", "Spleen/Liver", "Intestine", "Lung", "Intestine")
)
target2= data.frame("states"= names(list.state.m$top50$pert),
                   "organ"= c("Pan", "Pan", "Spleen/Lymph", "Bone", "Lung", "Artery/Tendon", "Intestine","Lung" ,"Muscle", "Panfibrosis" )
)
saveRDS(rbind(target1, target2), "data/prior_knowledge/fibro_atlas_deg/cross_organ_fibro_target.rds")

## 3.  DEG genes from contrast hf-ct in cell states:
intersect= readRDS("output/DEG_downsampled_intersect_hfpef.rds")

# de.list = readRDS(file = "output/cell_specific_obj/cell_type_based_DE/cells_DEA.rds")
de.subsampled= readRDS(file = "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")

# define number of DEGs
cutoff= 50
genes_de.fibs = de.subsampled$total$HFpEF #%>%
#   arrange((p_val_adj)) %>%
#   filter(avg_log2FC>0) %>%
#   slice(1:cutoff) %>%
#   pull(gene)

genes_de.fibs= intersect$fibroblasts$up

##4. IMPC cardiac phenotypes:
impc= readRDS("data/prior_knowledge/IMPC_processed.rds")

##5. cell state marker
Idents(seu)= "cellstate"
#seu=FindVariableFeatures(seu)
marker.fibs= FindAllMarkers(seu, assay= "RNA")

saveRDS(marker.fibs,  "output/fib_integration/marker_list/local_fib_state_marker.rds")

marker.fibs= readRDS( "output/fib_integration/marker_list/local_fib_state_marker.rds")

cut.off.state= 75
marker.fibs.genes= marker.fibs%>% group_by(cluster) %>%
  filter(p_val_adj <0.05,
         avg_log2FC>0)%>%
  slice(1:cut.off.state)

marker.fibs.genes= split(marker.fibs.genes$gene, marker.fibs.genes$cluster)

map(marker.fibs.genes, function(x){x[x%in% genes_de.fibs]})

# Module score --------------------------------------------------------------------------------
# Seurat's module score function is comparing average expression of a gene set to a control set of genes.

seu =add_geneset_scores(seu, de.subsampled$total)
seu =add_geneset_scores(seu,intersect$fibroblasts)

#map(paste0(names(de.subsampled$total), "1"), function(x) (fib.deg.feat= FeaturePlot(seu,x )))

meta= seu[[]]

p.deg.up=  ggplot(data= meta, aes(x= cellstate, y= up1,fill =group))+
  scale_fill_manual(values= hfpef_cols)+
  geom_boxplot()+
  theme_minimal()+
  #ggtitle(paste0("top ", cutoff, " genes from DEG" ))+
  ggtitle("upregulated DEG")+
  theme(axis.text.x = element_text(angle= 40, hjust=1),
        axis.text= element_text(color ="black"),
        axis.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
        plot.title = element_text(size = 12)
        )

p.deg.dn=  ggplot(data= meta, aes(x= cellstate, y= dn1,fill =group))+
  scale_fill_manual(values= hfpef_cols)+
  geom_boxplot()+
  theme_minimal()+
  ggtitle("downregulated DEG")+
  #ggtitle(paste0("top ", cutoff, " genes from DEG" ))+
  theme(axis.text.x = element_text(angle= 40, hjust=1),
        axis.title = element_blank(),
        axis.text= element_text(color ="black"),
        plot.title = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


p.module.scores= cowplot::plot_grid(p.deg.up, p.deg.dn, rel_widths = c(1, 1.25))
#FeaturePlot(seu,features = intersect$fibroblasts$up[21:25], split.by = "group", keep.scale = NULL)

df.test= meta%>%
  rownames_to_column(., "cellid") %>% as_tibble()%>%
  select(cellstate ,up1,cellid,dn1,   group)%>%
  group_by(group, cellstate)

## FC
df.test%>%
  summarise(m.up = median(up1),
            m.dn = median(dn1))%>%
  pivot_wider(names_from = group, values_from= c(m.up, m.dn))%>%
  mutate(fc.up = log2(m.up_hf/m.up_ct),
         fc.dn = log2(m.dn_hf/m.dn_ct))%>%
  mutate(fc.up =m.up_hf-m.up_ct,
         fc.dn = m.dn_hf-m.dn_ct)

df.test2= df.test %>%
  pivot_wider(names_from = group, values_from= c(up1, dn1))%>%
  group_by(cellstate)
df.test2%>%
  summarise(stat.up= wilcox.test(up1_ct,up1_hf)$statistic,
            p.up= wilcox.test(up1_ct,up1_hf)$p.value,
            stat.dn= wilcox.test(dn1_ct,dn1_hf)$statistic,
            p.dn= wilcox.test(dn1_ct,dn1_hf)$p.value)
df.test2%>%
  ggplot(aes(x= cellstate, y= stat.up))+
    geom_col()

#wilcox effect size

df.list= map(unique(meta$cellstate), function(x){
  print(x)
  y= meta%>% filter(cellstate== x)
  rbind(
  wilcox_effsize(y, formula = up1 ~ group) %>%
    mutate(cellstate= x,
           geneset= "up"),

  wilcox_effsize(y, formula = dn1 ~ group) %>%
    mutate(cellstate= x,
           geneset= "dn")
  )

})%>% do.call(rbind, .)


df.list %>% ggplot(aes(x= cellstate, y= effsize, fill= magnitude))+
  facet_grid(rows= vars(geneset))+
  geom_col()+
  scale_fill_manual(values= c("darkblue", "darkgreen", "darkorange"))+
  theme_bw()

## AUROC
library(pROC)
x= "Wif1+"
auc.df= sapply(unique(df.test$cellstate) ,function(x){
  df.test2= df.test %>% filter(cellstate==x)
  up.auc= auc( df.test2$group, df.test2$up1)
  dn.auc= auc( df.test2$group, df.test2$dn1, direction =  "<")
  c("up"= up.auc,"dn" =dn.auc)
})
colnames(auc.df)= unique(df.test$cellstate)
library(circlize)
col_fun = colorRamp2(c(0,0.5, 1), c("blue",  "white", "red"))
hmap.auroc= ComplexHeatmap::Heatmap(t(auc.df), col = col_fun, name= "AUROC", cluster_columns = F, cluster_rows = T)
hmap.auroc
p.cols= t(auc.df) %>% as.data.frame()%>% rownames_to_column("state")%>% pivot_longer(names_to= "set", values_to="auroc", -state)%>%
  ggplot(aes(state, auroc, fill= set))+
  geom_col(position = "dodge")+
  ylim(c(0.5, 1))

p.cols
pdf("output/figures/supp/auroc.deg.mod.fibstates.pdf",
    width= 2.5,
    height= 2)
hmap.auroc
dev.off()
## naba module scores:
seu= add_geneset_scores(seu, NABA_SETS_mouse)

meta= seu@meta.data
meta= meta %>%
  gather(paste0(names(NABA_SETS), "1"), key= "NABA", value= "ECM_score") %>%
  rename(cluster= opt_state)

p.naba.overview=  ggplot(data= meta, aes(x= NABA, y= ECM_score,fill = NABA))+
  geom_boxplot()+
  #theme_minimal()+
  #scale_fill_manual(values = col.set)+
  #scale_color_manual(values = col.set)+
  facet_grid(rows= vars(cluster))+
  coord_flip()

p.naba.select=  ggplot(data= meta %>%
                         filter(!NABA %in% c("NABA_ECM_AFFILIATED1","NABA_SECRETED_FACTORS1", "NABA_MATRISOME_ASSOCIATED1")),
                       aes(x= NABA, y= ECM_score,fill = NABA))+
  geom_boxplot()+
  #theme_minimal()+
  #scale_fill_manual(values = col.set)+
  #scale_color_manual(values = col.set)+
  facet_grid(rows= vars(cluster))+
  coord_flip()

### add fibro.atlas scores:

fib.marker= readRDS("data/prior_knowledge/fibro_atlas_deg/fibro_atlas_steady_state_marker.rds")

de.genes= map(fib.marker, function(x){
  x%>% arrange(p_val_adj)%>% slice(1:50)  %>% pull(Gene)
})

# add module scores and plot:
seu= add_geneset_scores(seu, de.genes)
meta= seu@meta.data
meta= meta %>%
  gather(paste0(names(de.genes),"1"), key= "fib.cluster.steady.state", value= "fib_score") %>%
  rename(cluster= opt_state)

p.fib.overview=  ggplot(data= meta, aes(x= cluster, y= fib_score,fill = cluster))+
  geom_boxplot()+
  theme_minimal()+
  #scale_fill_manual(values = col.set)+
  #scale_color_manual(values = col.set)+
  facet_grid(cols= vars(fib.cluster.steady.state))

p.fib.overview


p.fib.overview2=  ggplot(data= meta %>% group_by(cluster)%>% mutate(clust_mean= mean(fib_score)),
                         aes(x= fib.cluster.steady.state, y= fib_score,fill = fib.cluster.steady.state))+
  geom_boxplot()+
  #theme_minimal()+
  #scale_fill_manual(values = col.set)+
  #scale_color_manual(values = col.set)+
  facet_grid(rows= vars(cluster))+
  coord_flip()
p.fib.overview2

## impc?

impc= impc[lapply(impc, function(x) (length(x)>2))%>% unlist]

seu =add_geneset_scores(seu, impc)
names(impc)= str_replace_all(names(impc), " ", ".")
meta= seu@meta.data
meta= meta %>%
  gather(paste0(names(impc),"1"), key= "impc_set", value= "score") %>%
  rename(cluster= opt_state)

p.fib.overview=  ggplot(data= meta, aes(x= cluster, y= score,fill = cluster))+
  geom_boxplot()+
  theme_minimal()+
  #scale_fill_manual(values = col.set)+
  #scale_color_manual(values = col.set)+
  facet_grid(cols= vars(impc_set))
p.fib.overview
# ORA  --------------------------------------------------------------
## func for plotting:

plot_ORA=function(Ora_res , filter.names){

  for (i in names(Ora_res)){
    Ora_res[[i]] = Ora_res[[i]] %>% mutate(cluster= i)
  }

  Ora_res2= as_tibble(do.call(rbind, Ora_res)) %>%
    #rename(cluster= gset) %>%
    mutate(corr_p_value = p.adjust(p_value, method= "BH"))%>%
    mutate(stars= ifelse(corr_p_value<0.05, "*", ""),
           stars= ifelse(corr_p_value<0.01, "**", stars),
           stars= ifelse(corr_p_value<0.001, "***",stars))%>%
    mutate(cluster= factor(cluster, levels =  c("Col15a1+", "Igfbp3+",
                                                "Pi16+","Cxcl1+","Cilp+","Wif1+")))

  if (is.null(filter.names)){
    filter.names= Ora_res2$gset
  }

  Ora_res2 = Ora_res2 %>% filter(gset %in% filter.names)

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

#1. Cell state marker from fibroblast atlas

#combined.list= c("steady"=  stead.genes, "perturbed"= pert.genes)#, de.subsampled$unique)
combined.list= c(stead.genes, pert.genes)#, de.subsampled$unique)
p.list= lapply(marker.fibs.genes, function(x){
  Ora_res= GSE_analysis(geneList = x, Annotation_DB = combined.list)

})
ps= plot_ORA(p.list)

pdf("output/figures/fibroblast.state.atlas.pdf",
    height= 5,
    width= 5.0)
ps
dev.off()

#2. NABA

NABA_select= (NABA_SETS_mouse[!names(NABA_SETS_mouse) %in%
                                c("Secreted_factors","Matrisome", "Matrisome_associated", "Ecm_regulators", "Ecm_affiliated")])
Ora_res= lapply(marker.fibs.genes, function(x){
                                    GSE_analysis(geneList = x,Annotation_DB = NABA_SETS_mouse)

                                  })

p.ORA.naba= plot_ORA(Ora_res, filter.names = names(NABA_select))
saveRDS(p.ORA.naba,"output/figures/main/Fig2/ecm.ora.rds")

pdf("output/figures/fibroblast.NABA.pdf",
    height= 2.5,
    width= 4.5)
p.ORA.naba
dev.off()


#3. DEG
#get DEGs from hfpef vs ct in fibroblasts:

cutoff= 100

## get fib marker from atlas:
gene_signatures= readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled2.rds")
names(gene_signatures$total)= c("HFpEF", "MI", "AngII")
Ora_res= lapply(marker.fibs.genes, function(x){
  GSE_analysis(geneList = x,Annotation_DB = gene_signatures$total)

})

p.deg= plot_ORA(Ora_res)
pdf("output/figures/fibroblast.deg.hmap.pdf",
    height= 2.5,
    width= 4.5)
p.deg
dev.off()

# only hfpef gene set
Ora_res= lapply(marker.fibs.genes, function(x){
  GSE_analysis(geneList = x,Annotation_DB =gene_signatures$total)

})

p.deg= plot_ORA(Ora_res)

#4.

ora_impc= Ora_res= lapply(marker.fibs.genes, function(x){
    GSE_analysis(x, impc)
})



