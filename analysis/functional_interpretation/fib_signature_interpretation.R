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

source("analysis/utils_funcomics.R")
source("analysis/utils.R")

msigDB= readRDS("/home/jan/R-projects/sc_hfpef/data/prior_knowledge/Genesets_Dec19.rds")
gene_signatures= readRDS( "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
gene_translate= readRDS( "data/prior_knowledge/gene_translate.rds")


# data prep -----------------------------------------------------------------------------------

clean_names= function(string){
  clean.names= lapply(strsplit(string, "_" ),function(x){
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
msig_sigs$gset= clean_names(unlist(msig_sigs$gset))

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


msigDB_m$MSIGDB_HMARKS[grepl("Angptl4", msigDB_m$MSIGDB_HMARKS)]
length(msigDB_m$MSIGDB_REACTOME)


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

top_tf= df %>% filter(p.adj<0.05, condition== "HFpEF")%>% arrange(desc((score))) %>% #slice(1:50)%>%
  pull(source)%>% unique()

p.TFs= df %>% group_by(condition)%>%
  arrange(desc(score))%>% #%>% select(-p_value)%>% pivot_wider(values_from = score, names_from = condition)%>%
  mutate( sig= ifelse(p.adj<0.05, "*", ""),
    sig= ifelse(p.adj<0.01, "**", sig),
    sig= ifelse(p.adj<0.001, "***", sig),
    source= factor(source, levels=rev(top_tf)))%>%
  filter(source %in% top_tf)

p.pTFS= p.TFs %>%
  ggplot(., aes(x= condition, y= source, fill = score))+
  geom_tile()+
  geom_text(aes(label= sig))+
  theme_minimal()+
  theme(axis.text= element_text(color= "black"),
        panel.border = element_rect(colour = "black",fill=NA, size=0.5))+
  scale_fill_gradient2(low= "blue", mid= "white", high= "darkred")+
  labs(x= "", y= "", fill = "TF activity")

p.TFs
df = map(names(logFC.c), function(x){
  mat= logFC.c %>% filter(study== x)%>%select(gene, avg_log2FC)
  mat= as.matrix(column_to_rownames(mat, "gene"))

  df = decouple(mat = mat, network = net, statistics = "mlm")
  df %>% mutate(study= x)
})
geneS= dorothea_mm%>% filter(tf=="Ppara")%>% pull(target)
logFC.c %>% filter(gene %in% geneS)

# progeny -------------------------------------------------------------------------------------
progeny
source("analysis/utils_funcomics.R")


library(progeny)
gex= logFC.c %>%
  select(-pct.1, -pct.2) %>%
  pivot_wider(names_from = study, values_from= avg_log2FC) %>%
  distinct(gene, HFpEF, MI, AngII)

gex= as.matrix(column_to_rownames(as.data.frame(gex), "gene"))

M.Progeny = run_progeny(gex, .label = colnames(gex))

Heatmap(t(M.Progeny), cluster_columns = F, name = "Progeny_score")


# gset plot -----------------------------------------------------------------------------------
genes= msigDB_m$MSIGDB_CANONICAL$PID_INTEGRIN1_PATHWAY
PID_INTEGRIN1_PATHWAY
genes

df= logFC.c%>% filter(grepl("Itg", gene))%>% select(-pct.1, -pct.2) %>%
  pivot_wider(names_from = study,values_from = avg_log2FC)%>% column_to_rownames("gene")
df= logFC.c%>% dplyr::filter(gene %in% genes)%>% select(-pct.1, -pct.2) %>%
  pivot_wider(names_from = study,values_from = avg_log2FC)%>% column_to_rownames("gene")

Heatmap(df,name = "avg_log2FC")
