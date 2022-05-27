## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-03-11
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Interpretation of disease signatures

library(tidyverse)
library(decoupleR)
library(ggrepel)
library(ggplot2)
library(cowplot)
library(rlang)

library(ComplexHeatmap)
library(progeny)



source("code/utils_funcomics.R")
source("code/utils.R")

msigDB= readRDS("/home/jan/R-projects/sc_hfpef/data/prior_knowledge/Genesets_Dec19.rds")
gene_signatures= readRDS( "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
gene_translate= readRDS( "data/prior_knowledge/gene_translate.rds")


# data prep -----------------------------------------------------------------------------------

clean_names= function(string2){
  clean.names= lapply(strsplit(string2, "_" ),function(x){
    #print(x)
    xy= x[2:length(x)]
    #print(xy)
    paste(unlist(xy), collapse = "_")
  } )
  return(unlist(clean.names))
}

#select Msig DB genes of interest and translate genesets to mouse
names(msigDB)
msigDB_m= lapply(msigDB[c("MSIGDB_HMARKS",
                          "MSIGDB_REACTOME",
                          "MSIGDB_TF" ,
                          "MSIGDB_CANONICAL",
                          "MSIGDB_KEGG", "MSIGDB_BIOCARTA" )], function(y){
  lapply(y, function(x){
    gene_translate%>% filter(Gene.name %in% x)%>% pull(MGI.symbol)
  })

})


msigDB_m %>% saveRDS(., "data/prior_knowledge/msigDB.mouse_translated.rds")
msigDB_m <-readRDS( "data/prior_knowledge/msigDB.mouse_translated.rds")
# msigDB_H_m= lapply(msigDB$MSIGDB_HMARKS, function(x){
#   gene_translate%>% filter(Gene.name %in% x)%>% pull(MGI.symbol)
# })
#
# msigDB_H_m= lapply(msigDB$MSIGDB_REACTOME, function(x){
#   gene_translate%>% filter(Gene.name %in% x)%>% pull(MGI.symbol)
# })



# perform ORA ---------------------------------------------------------------------------------

## ora with both sets:
msig_sigs= map(names(gene_signatures$unique), function(x){
  lapply(names(msigDB_m), function(y){
    GSE_analysis(gene_signatures$total[[x]],msigDB_m[[y]])%>%
      mutate(study=x,
             collection = y)

  })%>% do.call(rbind,. )%>% as_tibble()
}) %>% do.call(rbind, .) %>% as_tibble()

## clean names:
#msig_sigs$gset= clean_names(unlist(msig_sigs$gset))

plot_set= clean_names(plot_set)

hfpef_hits=msig_sigs %>%
  filter( study== "HFpEF", corr_p_value<0.01)%>%
  arrange(corr_p_value)%>% pull(gset)%>% unique()

msig_sigs %>%
  filter( study== "HFpEF", corr_p_value<0.01)%>%
  arrange(corr_p_value)%>%
  print(n=100)

angii_hits=msig_sigs %>%
  filter( study== "AngII", corr_p_value<0.01)%>%
  arrange(corr_p_value)%>%  pull(gset)%>% unique()

# hfpef_hits= msig_sigs%>% filter(grepl("INTEGRIN", gset))%>% filter( study== "HFpEF", corr_p_value<0.1)%>%
#   arrange(corr_p_value)%>% pull(gset)%>% unique()

msig_sigs %>% filter(gset %in% hfpef_hits)%>%
  mutate(gset= factor(gset, levels= rev(hfpef_hits)))%>%
  mutate(label = paste0(GenesInList,  "/", GenesInPathway))%>%
  ggplot(., aes(x= study, y= gset, fill = -log10(corr_p_value)))+
  geom_tile()+
  scale_fill_gradient(low= "white", high= "darkred")+
  geom_text(aes(label= label))+
  theme_minimal()+
  theme(axis.text = element_text(colour = "black"))

# select gsets with biggest difference to ANgII
hfpef_hits= msig_sigs%>%
  pivot_wider(id_cols = c(study, corr_p_value, gset), names_from= study, values_from= corr_p_value) %>%
  unnest() %>%
  mutate(diff= -log10(HFpEF)+log10(AngII))%>% arrange(desc(diff)) %>% filter(gset %in% hfpef_hits)%>% pull(gset)%>% unique()




# final plot set manually curated from results above
plot_set= c("HALLMARK_ANGIOGENESIS",
            "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION" ,
            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" ,
            "NABA_MATRISOME",
            "REACTOME_ELASTIC_FIBRE_FORMATION" ,
            "REACTOME_COLLAGEN_DEGRADATION",
            "REACTOME_CROSSLINKING_OF_COLLAGEN_FIBRILS",
            "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES",
            "REACTOME_SIGNALING_BY_PDGF",
            "REACTOME_COLLAGEN_FORMATION",
            "NABA_ECM_REGULATORS",
            "REACTOME_LAMININ_INTERACTIONS",
            "PID_INTEGRIN1_PATHWAY",
            "NABA_BASEMENT_MEMBRANES" ,
            "REACTOME_O_GLYCOSYLATION_OF_TSR_DOMAIN_CONTAINING_PROTEINS",
            "REACTOME_FOCAL_ADHESION",
            "REACTOME_MOLECULES_ASSOCIATED_WITH_ELASTIC_FIBRES"
            #"REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
            #"CANONICAL_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX"
            )

plot_set= clean_names(plot_set)


p.msig= msig_sigs %>%
   mutate(corr_p_value = p.adjust(p_value, method = "BH"))%>% #we correct again as all test are done now
   filter(gset %in% plot_set)%>%
   mutate(#gset= factor(gset, levels= rev(plot_set)),
         sig= ifelse(corr_p_value<1e-2, "*", ""),
       sig= ifelse(corr_p_value<1e-5, "**", sig),
       sig= ifelse(corr_p_value<1e-10, "***", sig))%>%
  mutate(label = paste0(GenesInList,  "/", GenesInPathway),
         gset= str_to_title(gset))%>%
  ggplot(., aes(x= study, y= reorder(gset, -log10(corr_p_value)), fill = -log10(corr_p_value)))+
  geom_tile()+
  scale_fill_gradient(low= "white", high= "darkred")+
  geom_text(aes(label= sig))+
  theme_minimal()+
  theme(axis.text = element_text(colour = "black"))+
  labs(y= "", x= "", fill = "-log10(q-value)")
p.msig

pdf("output/figures/integration_studies/deg.hmaps/deg_msig.pdf",
    width= 6.5,
    height= 4)
p.msig
dev.off()


# TFs -----------------------------------------------------------------------------------------

library(dorothea)
data("dorothea_mm")

logFCs=readRDS(file = "output/fib_integration/marker_list/DEG_per_study_LOGFC.rds")
names(logFCs)=c("HFpEF", "MI", "AngII")
logFC.c= lapply(names(logFCs), function(x){
  rownames_to_column(logFCs[[x]], "gene") %>%
    mutate(study= x)
})%>% do.call(rbind, .)

net = dorothea_mm %>%
  filter(confidence %in% c("A", "B", "C")) %>% dplyr::rename(source= tf)%>% mutate(likelihood =1)%>% distinct(source, target, mor, likelihood)

mat= logFC.c %>% select(-pct.1, -pct.2)%>% pivot_wider(names_from= study, values_from = avg_log2FC) %>% drop_na
mat= as.matrix(column_to_rownames(mat, "gene"))

df = decouple(mat = mat, network = net, statistics = "ulm") %>%
  mutate(    p.adj = p.adjust(p_value, method = "BH"))

top_tf= df %>%
  filter(statistic== "ulm")%>%
  #group_by(condition)%>%
  filter(p.adj<0.05, condition== "HFpEF")%>%
  arrange(desc((score))) %>% #slice(1:50)%>%
  pull(source)%>% unique()

p.TFs= df %>% group_by(condition)%>%
  arrange(desc(score))%>% #%>% select(-p_value)%>% pivot_wider(values_from = score, names_from = condition)%>%
  mutate( sig= ifelse(p.adj<0.05, "*", ""),
    sig= ifelse(p.adj<0.01, "**", sig),
    sig= ifelse(p.adj<0.001, "***", sig),
    source= factor(source, levels=rev(top_tf)))%>%
  filter(source %in% top_tf)%>%
  filter(statistic== "ulm")

p.pTFS = p.TFs %>%distinct(condition, source, score, sig)%>%
  mutate(condition=factor(condition, levels= c("HFpEF", "MI", "AngII")) )%>%
  ggplot(., aes(x= condition, y= source, fill = score))+
  geom_tile()+
  geom_text(aes(label= sig))+
  scale_fill_gradient2(low= "blue", mid= "white", high= "red")+
  theme_minimal()+
  theme(axis.text= element_text(color= "black"),
        panel.border = element_rect(colour = "black",fill=NA, size=1),
        axis.text.x=element_text(angle= 40, hjust= 1))+
  labs(x= "", y= "", fill = "TF activity")+
  coord_equal()

p.pTFS
unify_axis(p.pTFS)
pdf("output/figures/main/Fig4/TF_study_wise.pdf",
    height= 5,
    width= 7)
unify_axis(p.pTFS)
dev.off()


df = map(names(logFC.c), function(x){
  mat= logFC.c %>% filter(study== x)%>%select(gene, avg_log2FC)
  mat= as.matrix(column_to_rownames(mat, "gene"))

  df = decouple(mat = mat, network = net, statistics = "mlm")
  df %>% mutate(study= x)
})
geneS= dorothea_mm%>% filter(tf=="Ppara")%>% pull(target)
logFC.c %>% filter(gene %in% geneS)

# progeny -------------------------------------------------------------------------------------
#progeny
source("code/utils_funcomics.R")
gex= logFC.c %>%
  select(-pct.1, -pct.2) %>%
  pivot_wider(names_from = study, values_from= avg_log2FC) %>%
  distinct(gene, HFpEF, MI, AngII)

gex= as.matrix(column_to_rownames(as.data.frame(gex), "gene"))

M.Progeny = run_progeny(gex, .label = colnames(gex))

hmap= Heatmap(t(M.Progeny), cluster_columns = F, name = "Progeny_score")
hmap= Heatmap(t(M.Progeny),cluster_rows = T,

        cluster_columns = F,
        name = "progeny \n score",
        column_names_rot = 40,
        border = T,
        #rect_gp = gpar(ol = "black", lty = 1),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10))

pdf("output/figures/main/Fig4/progeny_study_contrast.pdf",
    width= 2.6,
    height= 4)
print(hmap)
dev.off()


# gset plot -----------------------------------------------------------------------------------
genes= msigDB_m$MSIGDB_CANONICAL$PID_INTEGRIN1_PATHWAY
PID_INTEGRIN1_PATHWAY
genes

df= logFC.c%>% filter(grepl("Itg", gene))%>% select(-pct.1, -pct.2) %>%
  pivot_wider(names_from = study,values_from = avg_log2FC)%>% column_to_rownames("gene")
df= logFC.c%>% dplyr::filter(gene %in% genes)%>% select(-pct.1, -pct.2) %>%
  pivot_wider(names_from = study,values_from = avg_log2FC)%>% column_to_rownames("gene")

Heatmap(df,name = "avg_log2FC")



# calculate gset modules ----------------------------------------------------------------------

### order gene sets:
#1 get enriched genesets per cluster
#2 calculate jaccard distancs of enriched genes (overlap) between genesets
#3 hclust to define similar gene set groups
#4 plot each gene set group separately

# test an idea for geneset correlation
studies= unique(msig_sigs$study)
trees= map(studies, function(x){
  df=msig_sigs %>%
    filter( corr_p_value<0.05,
            study== x)%>%
    arrange(corr_p_value)%>%
    group_by(gset)%>%
    mutate(corr_p_value= min(corr_p_value))%>%
    arrange(corr_p_value)%>%
    distinct(gset, .keep_all = T)

  sets= str_split(df$GeneNames, ",")
  r.names= df$gset

  onehot.sets= sets %>%
    qdapTools::mtabulate()%>%
    qdapTools::matrix2df("genes") %>%
    column_to_rownames("genes")

  rownames(onehot.sets)= r.names

  jacc.dist= dist((onehot.sets), method= "binary")
  #ComplexHeatmap::Heatmap(jacc.dist.phecodes)
  c.x= hclust(jacc.dist, method = "ward.D2")

  res= sapply(c(1:10), function(y){
    enframe(cutree(c.x, k = y))$value
    # tree_labels = cutree(c.x, k = y)%>%
    #   as_tibble() %>%
    #   #mutate(gset= c.x$labels)%>%
    #   rename(!!(paste0("res.", y)) := value)
    #
  })
  colnames(res)= c(paste0("res.", c(1:10)))
  res.= res %>% as_tibble()%>%
    mutate(gset= c.x$labels)

  return(list(.tree= c.x,
              resolutions.= res.))
})
names(trees)= studies
plot(trees$HFpEF$.tree, labels = F)

hfpef= map(unique(trees$HFpEF$resolutions.$res.7), function(x){
   msig_sigs%>%
    left_join(trees$HFpEF$resolutions.)%>%
    filter(res.7==x)%>%
    mutate(#gset= factor(gset, levels= rev(plot_set)),
      sig= ifelse(corr_p_value<0.05, "*", ""),
      sig= ifelse(corr_p_value<0.01, "**", sig),
      sig= ifelse(corr_p_value<0.001, "***", sig),
      gset= clean_names(gset))%>%
    mutate(study= factor(study, levels= c("HFpEF", "MI", "AngII")))%>%
    ggplot(., aes(x=  study,
                  y= gset, fill= -log10(corr_p_value)))+
    geom_tile()+
    scale_fill_gradient2(low= "white" , high= "red" )+
    # scale_fill_gradientn(colors= c("white", "red"),
    #                       breaks=c(0,2,4,6,8,10),#labels=c("Minimum",0.5,"Maximum"),
    #                      limits=c(0,10))+
     geom_text(mapping = aes(label= sig))+
    #facet_grid(rows=vars(value))+
    theme_minimal()+
    labs(fill= "-log10(q-value)")+
    theme(axis.text.x= element_text(angle=40, hjust= 1, size= 10),
          axis.title = element_blank(),
          axis.text= element_text(color ="black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )#+coord_equal()
})


angii= map(unique(trees$AngII$resolutions.$res.7), function(x){
  msig_sigs%>%
    left_join(trees$AngII$resolutions.)%>%
    filter(res.7==x)%>%
    mutate(#gset= factor(gset, levels= rev(plot_set)),
      sig= ifelse(corr_p_value<0.05, "*", ""),
      sig= ifelse(corr_p_value<0.01, "**", sig),
      sig= ifelse(corr_p_value<0.001, "***", sig),
      gset= clean_names(gset))%>%
    mutate(study= factor(study, levels= c("HFpEF", "MI", "AngII")))%>%
    ggplot(., aes(x=  study,
                  y= gset, fill= -log10(corr_p_value)))+
    geom_tile()+
    scale_fill_gradient2(low= "white" , high= "red")+
    geom_text(mapping = aes(label= sig))+
    #facet_grid(rows=vars(value))+
    theme_minimal()+
    labs(fill= "-log10(q-value)")+
    theme(axis.text.x= element_text(angle=40, hjust= 1, size= 10),
          axis.title = element_blank(),
          axis.text= element_text(color ="black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )#+coord_equal()
})
legend_b <- get_legend(
  hfpef[[2]]
)
p.1= cowplot::plot_grid(#hfpef[[1]]+ theme(legend.position = "none"),
                    hfpef[[7]]+ theme(legend.position = "none",
                                      axis.text.x= element_blank()),
                   hfpef[[3]]+ theme(legend.position = "none",
                                     axis.text.x= element_blank()),
                   hfpef[[5]]+ theme(legend.position = "none",
                                     axis.text.x= element_blank()),
                   hfpef[[6]] + theme(legend.position = "none",
                                       axis.text.x= element_blank()),
                   angii[[4]] + theme(legend.position = "none",
                                      axis.text.x= element_blank()),
                   align= "v", axis= "lr", ncol = 1, rel_heights = c(1,1.3,0.9, 1, 0.8))
p.1

p.2= plot_grid(p.1, legend_b, rel_widths = c(1,0.3))
pdf("output/figures/main/Fig4/msig.hfpef_all.pdf",
    width= 6.5,
    height= 5)
p.2
dev.off()
pdf("output/figures/main/Fig4/msig.hfpef_all2.pdf",
    width= 5,
    height= 2)
hfpef[[1]]
dev.off()



hfpef[[3]]
cowplot::plot_grid(hfpef[[1]],
                   hfpef[[7]],
                   hfpef[[3]],
                   hfpef[[5]],
                   hfpef[[6]],
                   align="v", ncol = 1)

pdf("output/figures/supp/msig.angii.pdf")
angii
dev.off()
### calculate not enriched gene jaccards but rather full geneset overlap
df=msig_sigs %>%
    filter( corr_p_value<0.05)%>%
    pull(gset)%>%
  unique()
msig_sigs_filt = lapply(msigDB_m, function(x){
  x[names(x) %in% df]
})

x= unlist(msig_sigs_filt, recursive = F)
#unique(length(names(x))
names(x)
names(x) = str_replace_all(names(x), pattern = "^.*\\.", replacement = "")

onehot.sets = x[unique(names(x))] %>%
    qdapTools::mtabulate()%>%
    qdapTools::matrix2df("genes") %>%
    column_to_rownames("genes")

rownames(onehot.sets)= unique(names(x))

length(rownames(onehot.sets))

length(unique(names(x)))
jacc.dist= dist((onehot.sets), method= "binary")
  #ComplexHeatmap::Heatmap(jacc.dist.phecodes)
  c.x= hclust(jacc.dist, method = "ward.D")
  plot(c.x, labels = F)

  res= sapply(c(5:20), function(y){
    enframe(cutree(c.x, k = y))$value
    # tree_labels = cutree(c.x, k = y)%>%
    #   as_tibble() %>%
    #   #mutate(gset= c.x$labels)%>%
    #   rename(!!(paste0("res.", y)) := value)
    #
  })
  colnames(res)= c(paste0("res.", c(5:20)))
  res.= res %>% as_tibble()%>%
    mutate(gset= c.x$labels)%>%
    select(gset, everything())


class(res.$res.7)
  map(unique(res.$res.10), function(x){
    p. =msig_sigs%>%
      filter(gset %in% df)%>%
      left_join(res., by= "gset")#%>%

    p.%>%dplyr::filter(res.10 == x)%>%
      ggplot(., aes(x=  study,
                    y= gset, fill= -log10(corr_p_value)))+
      geom_tile()+
      scale_fill_gradient2(low= "white" , high= "red")+
      #geom_text(mapping = aes(label= stars))+
      #facet_grid(rows=vars(value))+
      theme_minimal()+
      labs(fill= "-log10(q-value)")+
      theme(axis.text.x= element_text(angle=40, hjust= 1, size= 10),
            axis.title = element_blank(),
            axis.text= element_text(color ="black"),
            panel.border = element_rect(colour = "black", fill=NA, size=1)
      )+coord_equal()
  })

