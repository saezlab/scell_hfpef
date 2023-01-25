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

source("code/utils.R")
source("code/utils_funcomics.R")

integrated_data = readRDS( "output/seu.objs/study_integrations/harmony_fib_filt_2.rds")
msigDB= readRDS("data/prior_knowledge/msigDB.mouse_translated.rds")

# find marker  --------------------------------------------------------------------------------
Idents(integrated_data)= "opt_clust_integrated"
DefaultAssay(integrated_data)= "RNA"
# Get markers per cluster:
marker_int= FindAllMarkers(integrated_data, assay= "RNA", slot = "data")

#save rds
saveRDS(marker_int, "output/fib_integration/marker_list/integrated_marker.rds")
marker_int= readRDS("output/fib_integration/marker_list/integrated_marker.rds")

#save as xlss
markers.list= split(x=marker_int, f = marker_int$cluster)
WriteXLS(markers.list,
           ExcelFileName = paste0("output/fib_integration/marker_list/statemarker.xlsx")
  )

int.fib.marker= readRDS( "output/fib_integration/marker_list/integrated_marker.rds")
marker.fibs= readRDS( "output/fib_integration/marker_list/local_fib_state_marker.rds")

cut.off.state= 100
marker.fibs.genes= marker.fibs%>% group_by(cluster) %>%
  filter(p_val_adj <0.05,
         avg_log2FC>0)%>%
  slice(1:cut.off.state)
local.genesets= split(marker.fibs.genes$gene, marker.fibs.genes$cluster)

int.fib.marker= int.fib.marker%>% group_by(cluster) %>%
  filter(p_val_adj <0.05,
         avg_log2FC>0)%>%
  dplyr::slice(1:cut.off.state)
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

saveRDS(NABA_SETS_mouse, "data/prior_knowledge/")
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
plot_ORA=function(Ora_res , gsets= NULL){

  for (i in names(Ora_res)){
    Ora_res[[i]] = Ora_res[[i]] %>% mutate(cluster= i)
  }

  Ora_res2= as_tibble(do.call(rbind, Ora_res)) %>%
    #rename(cluster= gset) %>%
    mutate(corr_p_value = p.adjust(p_value, method= "BH"))%>%
    mutate(stars= ifelse(corr_p_value<0.05, "*", ""),
           stars= ifelse(corr_p_value<0.01, "**", stars),
           stars= ifelse(corr_p_value<0.001, "***",stars))

  if(!is.null(gsets))
    {Ora_res2 = Ora_res2%>% filter(gset %in% gsets)}


  p.p.ora= ggplot(Ora_res2, aes(x=  cluster,
                                y= gset, fill= -log10(corr_p_value)))+
    geom_tile(color= "black")+
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

NABA_select= NABA_SETS_mouse[!names(NABA_SETS_mouse) %in% c("Ecm_affiliated","Ecm_regulators", "Secreted_factors","Matrisome")]
Ora_res= lapply(int.genesets, function(x){
  GSE_analysis(x, NABA_SETS_mouse)

})

p.naba= plot_ORA(Ora_res, c("Basement_membranes"  ,
                            "Collagens",
                            "Core_matrisome",
                            "Ecm_glycoproteins"
                            )
                 )
unify_axis(p.naba)
pdf("output/figures/main/Fig2/fibroblast_naba_ORA.pdf",
    width = 4.9,
    height= 3.5)
unify_axis(p.naba)+
  theme(axis.title =  element_blank())+
  coord_equal()

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

names(pert.state)= str_replace_all(names(pert.state), pattern = "AllMarkers_", "")
names(steady.genes)= str_replace_all(names(steady.genes), pattern = "AllMarkers_", "")
Ora_res.steady= lapply(int.genesets, function(x){
  GSE_analysis(x , steady.genes)

})

Ora_res.pert= lapply(int.genesets, function(x){
  GSE_analysis(x, pert.state)

})

p.pert= plot_ORA(Ora_res.pert)
p.stead= plot_ORA(Ora_res.steady)

pdf("output/figures/supp/fibroblast_integration_state_ORA.pdf",
    width = 5,
    height= 5)
cowplot::plot_grid(p.pert+ggtitle("perturbed state marker"),
                   p.stead+ggtitle("steady state marker"),
                   #labels = "AUTO",
                   ncol=1)
dev.off()

##4 marker compar

Ora_res= lapply(int.genesets, function(x){
  GSE_analysis(x, c(local.genesets,gene_signatures$total) )

})

iont
marker.comp= plot_ORA(Ora_res)
saveRDS(marker.comp, "output/figures/main/Fig3/marker.comp.rds")

# MSIG_DB -------------------------------------------------------------------------------------
msigDB_m= msigDB

Ora_res= lapply(names(msigDB_m), function(y){
  lapply(names(int.genesets), function(x){
    GSE_analysis(int.genesets[[x]],Annotation_DB = msigDB_m[[y]])%>%
      mutate(msig_group= y,
             cluster= x)
    })%>% do.call(rbind,. )%>% as_tibble()
})%>% do.call(rbind,. )%>% as_tibble()


Ora_res.s= Ora_res %>% select(gset, GeneNames, p_value, corr_p_value, cluster, msig_group)

# merge with GO res:
GO.res= readRDS("output/GO_cellstates_int.rds")
GO.res$mf= GO.res$mf %>% rename(gset= Term,
                     #n.genes = Overlap,
                  p_value= P.value,
                  corr_p_value= Adjusted.P.value ,
                  GeneNames= Genes)%>%
  mutate(msig_group = "GO_MF")%>%
  select(gset,GeneNames, p_value,corr_p_value, cluster, msig_group)
GO.res$bp= GO.res$bp %>% rename(gset= Term,
                                #n.genes = Overlap,
                                p_value= P.value,
                                corr_p_value= Adjusted.P.value ,
                                GeneNames= Genes)%>%
  mutate(msig_group = "GO_BP")%>%
  select(gset,GeneNames, p_value,corr_p_value, cluster, msig_group)


Ora_res2= rbind(GO.res$mf, Ora_res.s, GO.res$bp)

## clean names:
#Ora_res$gset= clean_names(unlist(Ora_res$gset))

Ora_res2= Ora_res2 %>%
  distinct(gset, GeneNames, p_value, cluster,msig_group)%>%
  #rename(cluster= gset) %>%
  mutate(corr_p_value = p.adjust(p_value, method= "BH"))%>%
  mutate(stars= ifelse(corr_p_value<0.1, "*", ""),
         stars= ifelse(corr_p_value<0.05, "**", stars),
         stars= ifelse(corr_p_value<0.01, "***",stars))
Ora_res2$GeneNames= str_replace_all(Ora_res2$GeneNames, ";", ",")
vec= str_split(Ora_res2$GeneNames, pattern = ",")
vec= unlist(map(vec, length))
Ora_res2$n.genes= vec
Ora_res2

### order gene sets:
#1 get enriched genesets per cluster
#2 calculate jaccard distancs of enriched genes (overlap) between genesets
#3 hclust to define similar gene set groups
#4 plot each gene set group separately

Ora_res2 %>% group_by(cluster)%>% filter(corr_p_value<0.001)
Ora_res2 %>% filter(cluster==2)%>% arrange(corr_p_value) %>% print(n=100)

sigs= Ora_res2 %>% group_by(gset)%>%
  filter((corr_p_value<0.1))%>%
  arrange(corr_p_value)

listed.sigs= split(sigs, sigs$cluster)
WriteXLS(listed.sigs, "output/msigs_int_state_marker.xls")

sigs%>%
 ggplot(., aes(x=  cluster,
                              y= gset, fill= -log10(corr_p_value)))+
  geom_tile()+
  scale_fill_gradient2(low= "white" , high= "red")+
  geom_text(mapping = aes(label= stars))+
  facet_grid(rows=vars(msig_group))+
  theme_minimal()+
  labs(fill= "-log10(q-value)")+
  theme(axis.text.x= element_text(angle=40, hjust= 1, size= 10),
        axis.title = element_blank(),
        axis.text= element_text(color ="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )#+coord_flip()
p.p.ora


g.sets= c(#0  "REACTOME_ECM_PROTEOGLYCANS",
          #"KEGG_ECM_RECEPTOR_INTERACTION",
          "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
         # "REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES",
          "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES",
         #"REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
          "REACTOME_COLLAGEN_FORMATION",
         #"REACTOME_NON_INTEGRIN_MEMBRANE_ECM_INTERACTIONS",
         "REACTOME_NCAM1_INTERACTIONS",
        "cardiac epithelial to mesenchymal transition (GO:0060317)",
         "extracellular matrix disassembly (GO:0022617)",
         #"extracellular matrix organization (GO:0030198)",
         #"collagen fibril organization (GO:0030199)",
         "regulation of angiogenesis (GO:0045765)",
         "basement membrane organization (GO:0071711)",




          #1

          "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
         "REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS",
          #"PID_BMP_PATHWAY",
          #"NABA_ECM_GLYCOPROTEINS",
          #"BIOCARTA_GHRELIN_PATHWAY",
          "positive regulation of vascular endothelial growth factor production (GO:0010575)",
          "regulation of insulin-like growth factor receptor signaling pathway (GO:0043567)",
          "regulation of ossification (GO:0030278)",
          "regulated exocytosis (GO:0045055)",


         #2
         "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
         "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",
          "REACTOME_SIGNALING_BY_NOTCH4",
          "REACTOME_REGULATION_OF_APOPTOSIS",
          "REACTOME_REGULATION_OF_RUNX3_EXPRESSION_AND_ACTIVITY",
          "REACTOME_FOLDING_OF_ACTIN_BY_CCT_TRIC",
         "REACTOME_FORMATION_OF_TUBULIN_FOLDING_INTERMEDIATES_BY_CCT_TRIC",
          "antigen processing and presentation of exogenous peptide antigen via MHC class I (GO:0042590)",
         #"HALLMARK_HYPOXIA",

          #3
         "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
         "REACTOME_COLLAGEN_FORMATION",
         #"REACTOME_COLLAGEN_CHAIN_TRIMERIZATION",
         #"REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES",
         "REACTOME_ELASTIC_FIBRE_FORMATION",
         "REACTOME_SYNDECAN_INTERACTIONS",
         #"KEGG_ECM_RECEPTOR_INTERACTION",
         "KEGG_FOCAL_ADHESION",
         #"REACTOME_MOLECULES_ASSOCIATED_WITH_ELASTIC_FIBRES",
          "skeletal system development (GO:0001501)",
          "supramolecular fiber organization (GO:0097435)",
          "cell-matrix adhesion (GO:0007160)",
          "regulation of bone mineralization (GO:0030500)",

         #4
          "regulation of focal adhesion assembly (GO:0051893)",
          "cellular response to growth factor stimulus (GO:0071363)",


         #5
         "REACTOME_ELASTIC_FIBRE_FORMATION",
         #"HALLMARK_COAGULATION",
         #"REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION",
          "positive regulation of Wnt signaling pathway, planar cell polarity pathway (GO:2000096)",
         # "negative regulation of Wnt signaling pathway (GO:0030178)",
          "regulation of Wnt signaling pathway (GO:0030111)",
          "endochondral bone morphogenesis (GO:0060350)",

         #6
         "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        #"HALLMARK_MTORC1_SIGNALING",
         "HALLMARK_INFLAMMATORY_RESPONSE",
        # "REACTOME_SIGNALING_BY_INTERLEUKINS",
         "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",

          "cytokine-mediated signaling pathway (GO:0019221)",
          "positive regulation of cytokine production (GO:0001819)",
          "regulation of macrophage activation (GO:0043030)",
          "positive regulation of fibroblast proliferation (GO:0048146)",

         #7
         "HALLMARK_INTERFERON_GAMMA_RESPONSE",
         "REACTOME_INTERFERON_SIGNALING",
         "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
         "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
          "cytokine-mediated signaling pathway (GO:0019221)"
         # "cellular response to type I interferon (GO:0071357)"

         )

Ora_res3= Ora_res2%>%
  filter(gset %in% g.sets)%>%
  mutate(gset = factor(gset, levels= unique(g.sets)))

# # string cosmetic for plotting

clean_strings= function(g.sets){
  g.sets2=  gsub(pattern = "\\(.*",replacement = "",x =g.sets)
  g.sets2= str_remove_all(str_remove_all(str_remove_all(g.sets2, "KEGG_"), "REACTOME_"), "HALLMARK_")
  g.sets2= str_replace_all(g.sets2 , pattern = "_", " ")
  g.sets2= str_to_title(g.sets2)
  g.sets2

}
g.set2= clean_strings(unique(g.sets))

Ora_res3$gset= clean_strings(Ora_res3$gset)

p.p= Ora_res3%>%

  mutate(gset= factor(gset, levels= g.set2))%>%
  ggplot(., aes(x=  cluster,
                y= (gset), fill= -log10(corr_p_value)))+
  geom_tile(col = "grey")+
  scale_fill_gradient2(low= "white", mid = "red", high= "darkred", midpoint= 15)+
  #scale_fill_gradient(low= "white", high= "red")+
  geom_text(mapping = aes(label= stars))+
  theme_minimal()+
  labs(fill= "-log10(q-value)",
       x= "IFS", y= "")+
  theme(axis.text.x= element_text(angle=40, hjust= 1, size= 10),
        axis.title = element_blank(),
        axis.text= element_text(color ="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.y = element_text( colour = "black")
  )#+
#library(scales)
 #scale_y_discrete( labels = label_wrap(50))
# Ora_res3%>% filter(grepl("Elastic", gset), cluster== 4)

 p.p
pdf("output/figures/main/Fig2/gsea_statemarker.pdf",
    width= 15, height= 9)
unify_axis(p.p)+coord_equal()
dev.off()

pdf("output/figures/main/Fig2/gsea_statemarker.pdf",
    width=20, height= 6)
p.p+
  coord_flip()+
  theme(axis.text.x = element_text(angle=90, hjust= 1))


dev.off()



Ora_res2$gset
p.p.ora
 clean_names(df = Ora_res2, col = "gset")

df= sigs %>% filter(!grepl("NABA", gset))%>%
  group_by(cluster)%>%
  slice_max(order_by = -corr_p_value, n=10)%>% distinct(gset, cluster, corr_p_value, stars)

ggplot(df, aes(x= cluster, y= gset, fill = -log10(corr_p_value)))+
  geom_tile()

# gene set correlations -----------------------------------------------------------------------


# test an idea for geneset correlation
df=Ora_res2 %>% filter(cluster==2,
                       #           msig_group== "MSIGDB_HALLMARK",
                       corr_p_value<0.01)%>% arrange(corr_p_value)

df= df%>% group_by(gset)%>% mutate(corr_p_value= min(corr_p_value))%>% arrange(corr_p_value)%>% distinct(gset, .keep_all = T)
df$GeneNames
sets= str_split(df$GeneNames, ",")
names= df$gset

x= sets %>%
  qdapTools::mtabulate()%>%
  qdapTools::matrix2df("genes") %>%
  column_to_rownames("genes")

rownames(x)= names

jacc.dist.phecodes=
  dist((x), method= "binary") %>%
  as.matrix()

c.x= hclust(dist((x), method= "binary"), method = "ward.D2")
tree_labels = cutree(c.x, k = 5) %>% as_tibble() %>% mutate(gset= c.x$labels)
print(c.x)
tree_labels%>% arrange(value) %>% print(n=100)
map(unique(tree_labels$value), function(x){
  df%>%
    left_join(tree_labels)%>%
      filter(value==x)%>%
    ggplot(., aes(x=  cluster,
                  y= gset, fill= -log10(corr_p_value)))+
    geom_tile()+
    scale_fill_gradient2(low= "white" , high= "red")+
    geom_text(mapping = aes(label= stars))+
    #facet_grid(rows=vars(value))+
    theme_minimal()+
    labs(fill= "-log10(q-value)")+
    theme(axis.text.x= element_text(angle=40, hjust= 1, size= 10),
          axis.title = element_blank(),
          axis.text= element_text(color ="black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )#+coord_flip()
})
library(ComplexHeatmap)
p.disgenet.jacc= Heatmap(as.matrix(jacc.dist.phecodes), show_row_names = F, show_column_names = F, name= "Jaccard dist")
p.disgenet.jacc
