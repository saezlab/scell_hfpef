## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-05
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## using the single study DEGs, we will analyze the integrated cell states to explore which
## cell state is induced, and express the disease relevant signatures

library(Seurat)
library(tidyverse)

source("code/utils.R")
source("code/utils_funcomics.R")
int.fibs= readRDS("output/seu.objs/study_integrations/harmony_fib_filt.rds")
int.fibs@meta.data %>% group_by(study, group) %>% count
gene_signatures= readRDS( "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
state_marker= readRDS("output/fib_integration/marker_list/integrated_marker.rds")

# Feature Plots -------------------------------------------------------------------------------

heat.hfpef= DoHeatmap(int.fibs, features= signatures$unique$hfpef)
heat.circ= DoHeatmap(int.fibs, features= signatures$unique$circ)
heat.all= DoHeatmap(int.fibs, features= signatures$overlap$all)

pdf("output/figures/integration_studies/fibroblast.study.DEG.heatmaps.pdf",
    height= 10,
    width = 10)
heat.hfpef
heat.circ
heat.all
dev.off()
# Module Score --------------------------------------------------------------------------------

# add module scores and plot:
get_module_score_for_genesetlist= function(seu, geneset.list){
  int.fibs= add_geneset_scores(int.fibs,geneset.list )

  meta= int.fibs@meta.data

  meta= meta %>%
    gather(paste0(names(geneset.list),"1"), key= "fib.cluster.steady.state", value= "fib_score") %>%
    rename(cluster= opt_clust_integrated)

  p.fib.overview=  ggplot(data= meta, aes(x= cluster, y= fib_score,fill = cluster))+
    geom_boxplot()+
    theme_minimal()+
    #scale_fill_manual(values = col.set)+
    #scale_color_manual(values = col.set)+
    facet_grid(cols= vars(fib.cluster.steady.state))+
    theme(axis.text = element_text(size= 13),
          strip.text.x = element_text(size = 14, colour = "black")

          )



}

p.fib.overview_total = get_module_score_for_genesetlist(int.fibs,signatures$total )
p.fib.overview = get_module_score_for_genesetlist(int.fibs,c(signatures$unique, signatures$overlap))
p.fib.overview_unique = get_module_score_for_genesetlist(int.fibs,red_geneS )
red_geneS= gene_signatures$unique
red_geneS$MI= red_geneS$MI[1:50]
pdf("output/fib_integration/fibroblast_module_scores_deg.pdf",
    width= 10,
    height= 7)
p.fib.overview
p.fib.overview_total
dev.off()

# ORA of cell state marker --------------------------------------------------------------------

#prepare_cell state marker:
marker = state_marker %>% group_by(cluster) %>% filter(p_val_adj<0.05,
                                                       avg_log2FC>0)%>%
  top_n(n = 100, avg_log2FC)
marker=split(marker$gene, marker$cluster)

# ADD NABA
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
    dplyr::pull(MGI.symbol)

})

names(NABA_SETS_mouse)= str_replace_all(names(NABA_SETS_mouse),"NABA_", "")
names(NABA_SETS_mouse)= str_to_title(names(NABA_SETS_mouse))

## run ORA NABA
Ora_res_NABA= lapply(NABA_SETS_mouse[!names(NABA_SETS_mouse) %in%
                                       c("Secreted_factors","Proteoglycans", "Ecm_regulators", "Ecm_affiliated")], function(x){
  GSE_analysis(x, marker)

})

## run ORA fib signatures:
Ora_res= lapply(c(gene_signatures$unique, gene_signatures$overlap), function(x){
  GSE_analysis(x, marker)

})

Ora_res= lapply(c(gene_signatures$total), function(x){
  GSE_analysis(x, marker)

})


p.DEG_ORA = plot_ORA(Ora_res)
p.int.NABA= plot_ORA(Ora_res_NABA)
cowplot::plot_grid(p.int.NABA+theme(legend.position = "bottom"),
          p.DEG_ORA+theme(legend.position = "bottom"))

pdf(file = "output/figures/main/fibroblast.study.DEG.ORA.pdf",
    width = 4,
    height = 2)
p.DEG_ORA

dev.off()

pdf(file = "output/figures/main/fibroblast.study.ECM.ORA.pdf",
    width = 5,
    height = 1.7)
p.int.NABA

dev.off()



# compare signatures within  cellstates -------------------------------------------------------

int.fibs= add_geneset_scores(int.fibs,gene_signatures$total )

meta= int.fibs@meta.data

studies= unique(meta$study)
map(names(gene_signatures$total), function(y){
  map(studies, function(x){

    meta%>%
      filter(study== x)%>%
      ggplot(., aes(x= opt_clust_integrated, y= !!as.symbol(paste0(y, "1")), fill = group ))+
      geom_boxplot()+
      ggtitle(x)
  })
})

## add auroc for better comparison
df.test= meta%>%
  rownames_to_column(., "cellid") %>% as_tibble()%>%
  dplyr::select(opt_clust_integrated ,HFpEF1,AngII1, MI1, cellid, study ,group)%>%
  group_by(group, opt_clust_integrated)


library(pROC)

auc.df= lapply(unique(df.test$study) ,function(y){
  df= sapply(sort(unique(df.test$opt_clust_integrated)) ,function(x){
    df.test2= df.test %>% dplyr::filter(opt_clust_integrated==x, study== y)
    hfpef.auc= auc( df.test2$group, df.test2$HFpEF1)
    ang.auc= auc( df.test2$group, df.test2$AngII1)
    mi.auc= auc( df.test2$group, df.test2$MI1)
    c("HFpEF"= hfpef.auc,"AngII" =ang.auc, "MI"= mi.auc)
    })
  colnames(df)= sort(unique(df.test$opt_clust_integrated))
  return(df)

  })

names(auc.df)= unique(df.test$study)
auc.df

library(circlize)
col_fun = colorRamp2(c(0.5, 1), c( "white", "red"))
hmaps= map(auc.df, function(x){
  hmap.auroc= ComplexHeatmap::Heatmap(t(x), col = col_fun, name= "AUROC", cluster_columns = F, cluster_rows = F)

})
p.cols= t(auc.df$circ) %>% as.data.frame()%>% rownames_to_column("state")%>% pivot_longer(names_to= "set", values_to="auroc", -state)%>%
  ggplot(aes(state, auroc, fill= set))+
  geom_col(position = "dodge")#+
  ylim(c(0.5, 1))

p.cols
pdf("output/figures/supp/auroc.deg.mod.fibstates_across.pdf",
    width= 1.5,
    height= 2)
print(hmaps[[1]])
print(hmaps[[2]])
print(hmaps[[3]])
dev.off()
