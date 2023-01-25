## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-02-07
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## reheat enrich


library(ggpubr)
library(Seurat)
library(tidyverse)
library(decoupleR)
library(ggrepel)

source("code/utils.R")
source("code/utils_funcomics.R")

## sc data
int.fibs= readRDS("output/seu.objs/study_integrations/harmony_fib_filt.rds")
int.fibs@meta.data %>% group_by(study, group) %>% count
gene_signatures= readRDS( "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
state_marker= readRDS("output/fib_integration/marker_list/integrated_marker.rds")
marker = state_marker %>% group_by(cluster) %>% filter(p_val_adj<0.05,
                                                       avg_log2FC>0)%>%
  top_n(n = 100, avg_log2FC)
marker=split(marker$gene, marker$cluster)

## load study specific effect sizes
logFCs=readRDS(file = "output/fib_integration/marker_list/DEG_per_study_LOGFC.rds")
names(logFCs)=c("HFpEF", "MI", "AngII")

## das et.al 2019
bulk= readRDS("data/_Hfpef_bulk_study/HF_pEF_processed.rds")

## reheat data can be downloaded https://zenodo.org/record/3797044#.XsQPMy2B2u5
reheat= readRDS("/home/jan/R-projects/HF_meta-analysis/data/shiny/directed_ranking.rds")
reheat.p= readRDS("/home/jan/R-projects/HF_meta-analysis/data/shiny/fisher_rank.rds")
## gene translate
gene_translate= readRDS( "data/prior_knowledge/gene_translate.rds")
#translate reheat
names(reheat)
reheat.translated= enframe(reheat, name= "Gene.name") %>% left_join(gene_translate) %>% drop_na
reheat = reheat.translated %>% pull(value)
names(reheat)= reheat.translated$MGI.symbol

# Feature Plots -------------------------------------------------------------------------------

# add module scores and plot:
reheat_scores= list(
 "reheat_50"= names(reheat[1:50]),
  "reheat_100"= names(reheat[1:100]),
  "reheat_250"= names(reheat[1:250]),
  "reheat_500"= names(reheat[1:500])
)

int.fibs= add_geneset_scores(int.fibs,reheat_scores )

meta= int.fibs@meta.data

meta= meta %>%
  gather(paste0(names(reheat_scores),"1"), key= "fib.cluster.steady.state", value= "fib_score") %>%
  rename(cluster= opt_clust_integrated)

my_comparisons= list( c("3","7"))
meta$reheat2501
s= map( paste0(names(reheat_scores),"1"), function(x){
  print(x)
  p.fib.overview=  ggplot(data= meta, aes(x= reorder(cluster,-{{x}}), y= {{x}},fill = cluster))+
    geom_violin()+
    geom_boxplot(fill = "white",width = 0.4)+
    theme_minimal()+
    scale_fill_manual(values = col_vector)+
    #scale_color_manual(values = col.set)+
    #facet_grid(cols= vars(fib.cluster.steady.state))+
    theme(axis.text = element_text(size= 13),
          strip.text.x = element_text(size = 14, colour = "black"),
          legend.position = "none",
          axis.title.x= element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    )+
    labs(y= "module score")+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

})
s
p.fib.overview=  ggplot(data= meta, aes(x= reorder(cluster,-reheat2501), y= reheat2501,fill = cluster))+
  geom_violin()+
  geom_boxplot(fill = "white",width = 0.4)+
  theme_minimal()+
  scale_fill_manual(values = col_vector)+
  #scale_color_manual(values = col.set)+
  #facet_grid(cols= vars(fib.cluster.steady.state))+
  theme(axis.text = element_text(size= 13),
        strip.text.x = element_text(size = 14, colour = "black"),
        legend.position = "none",
        axis.title.x= element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
        )+
  labs(y= "module score")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

p.fib.overview

pdf("output/figures/integration_studies/reheat_module_scores.pdf",
    width= 4.2,
    height= 4)
p.fib.overview
dev.off()

FeaturePlot(int.fibs, features = "reheat501", pt.size = 0)

VLNS = VlnPlot(int.fibs, features = paste0(names(reheat_scores),"1"), pt.size = 0,
               same.y.lims = T, ncol= 5)

# ORA with  ReHEAT ----------------------------------------------------------------

##prepare gene set input
# disease signatures:
net= enframe(gene_signatures$total, name = "source", value = "target") %>% unnest(target) %>%
  mutate(mor = 1,
         likelihood= 1, )
#cell state marker:
net2= enframe(marker,  name = "source", value = "target") %>% unnest(target) %>%
  mutate(mor = 1,
         likelihood= 1, )

net3= rbind(net, net2)
## prepare reheat:
re.= as.matrix(reheat)

# decouple(re., net, statistics = "wmean")
# run_ora(re., net,n_up = 200, n_bottom = 0)

#ora
#1. signatures
#map different reheat cutoffs to perform ora with state marker
p.= map(c(100,200, 500), function(x){
  res= run_ora(re., net3,n_up = x, n_bottom = 0)
  res %>% mutate(n.reheat= x)
})%>% do.call(rbind, .)

p.signatures.hfref = p.%>%
  mutate(sig= ifelse(p_value<0.05, "*", ""),
         sig= ifelse(p_value<0.01, "**", sig),
         sig= ifelse(p_value<0.001, "***", sig))%>%
  ggplot(., aes(x= source, y= as.factor(n.reheat), fill= score))+
  geom_tile()+
  scale_fill_gradient(low= "white", high = "darkred")+
  geom_text(aes(label=sig))+
  theme_linedraw()+
  labs(x= "diseae signature",
       y= "top n ReHeaT genes")

#2. ora with state marker:
p.2= map(c(100,200, 500), function(x){
  res= ora.int.state= run_ora(re., net2,n_up = x, n_bottom = 0)
  res %>% mutate(n.reheat= x)
})%>% do.call(rbind, .)

p.states.hfref= p.2 %>%
  mutate(sig= ifelse(p_value<0.05, "*", ""),
         sig= ifelse(p_value<0.01, "**", sig),
         sig= ifelse(p_value<0.001, "***", sig))%>%
  ggplot(., aes(x= source, y= as.factor(n.reheat), fill= score))+
  geom_tile()+
  scale_fill_gradient(low= "white", high = "darkred")+
  geom_text(aes(label=sig))+
  theme_linedraw()+
  labs(x= "integrated fibroblast state",
       y= "top n ReHeaT genes")

# ORA with hfpef bulk study Das et al 2019 -----------------------------------------------------------------

#translate and prepare
dea= bulk$DEA%>%
  left_join(gene_translate%>% dplyr::rename(gene= Gene.name))%>%
  arrange(desc(t))

mat= as.matrix(dea$t)
rownames(mat)= dea$MGI.symbol

#wmean
wmean.res= decouple(mat, net, statistics = "wmean")
wmean.res%>%
  filter(statistic== "norm_wmean")%>%
ggplot(., aes(x= source, y= statistic, fill = score))+
  geom_tile()

# ora
# 1. disease signatures
#map different alpha levels for cut offs.
p.3= map(c(0.05, 0.01, 0.001), function(x){
  n.genes= length(dea%>% filter(adj.P.Val<x, t>0)%>% pull(MGI.symbol)%>% unique())
  print(n.genes)
  res= run_ora(mat, net3,n_up = n.genes, n_bottom = 0)
  #res=  run_ora(re., net,n_up = x, n_bottom = 0)
  res %>% mutate(alpha= x)
})%>% do.call(rbind, .)

p.signatures.hfpef= p.3 %>%
  mutate(sig= ifelse(p_value<0.05, "*", ""),
         sig= ifelse(p_value<0.01, "**", sig),
         sig= ifelse(p_value<0.001, "***", sig))%>%
  ggplot(., aes(x= source, y= as.factor(alpha), fill= score))+
  geom_tile()+
  scale_fill_gradient(low= "white", high = "darkred")+
  geom_text(aes(label=sig))+
  theme_linedraw()+
  labs(x= "disease signatures",
       y= "alpha")


# 2. cell state marker
#map different alpha levels for cut offs.
p.4= map(c(0.05, 0.01, 0.001), function(x){
  n.genes= length(dea%>% filter(adj.P.Val<x, t>0)%>% pull(MGI.symbol) %>% unique())
  print(n.genes)
  res= run_ora(mat, net2,n_up = n.genes, n_bottom = 0)
  #res=  run_ora(re., net,n_up = x, n_bottom = 0)
  res %>% mutate(alpha= x)
})%>% do.call(rbind, .)

p.marker.hfpef= p.4 %>%
  mutate(sig= ifelse(p_value<0.05, "*", ""),
         sig= ifelse(p_value<0.01, "**", sig),
         sig= ifelse(p_value<0.001, "***", sig))%>%
  ggplot(., aes(x= source, y= as.factor(alpha), fill= score))+
  geom_tile()+
  scale_fill_gradient(low= "white", high = "darkred")+
  geom_text(aes(label=sig))+
  theme_linedraw()+
  labs(x= "integrated fibroblast state",
       y= "alpha")+
  coord_equal()



# plot together -------------------------------------------------------------------------------
plot.df= rbind(p.%>% filter(n.reheat==500)%>% select(-n.reheat)%>% mutate(Human.ref= "HFrEF"),
      p.2%>% filter(n.reheat==500)%>% select(-n.reheat)%>% mutate(Human.ref= "HFrEF"),
      p.3 %>% filter(alpha== 0.05)%>% select(-alpha)%>% mutate(Human.ref= "HFpEF"),
      p.4 %>% filter(alpha== 0.05)%>% select(-alpha)%>% mutate(Human.ref= "HFpEF")
      )%>% mutate(q_value= p.adjust(p_value),
                  sig= ifelse(p_value<0.05, "*", ""),
                         sig= ifelse(p_value<0.01, "**", sig),
                         sig= ifelse(p_value<0.001, "***", sig))

max.p= max(-log10(plot.df$q_value))
p_d= plot.df %>% filter(source %in% c("MI_early", "MI_late", "HFpEF", "AngII"))%>%
  mutate(source= factor(source,levels= c(c("HFpEF","AngII","MI_early", "MI_late" ))))%>%
  ggplot(., aes(x= source, y= Human.ref, fill= score))+
  geom_tile(color= "darkgrey")+
  scale_fill_gradientn(colors= c("white", "red"),
                       breaks=c(0,2,4,6,8),#labels=c("Minimum",0.5,"Maximum"),
                       limits=c(0,10))+
  geom_text(aes(label=sig))+
  theme_minimal()+
  theme(axis.text= element_text(color= "black"),
        panel.border = element_rect(colour = "black",fill=NA, size=1),
        axis.text.x=element_text(angle= 40, hjust= 1))+
  labs(x= "Murine disease signature", y= "Human bulk transcriptome", fill = "-log10(q-value)")


p_d= unify_axis(p_d)

p_s= plot.df %>% filter(!source %in% c("MI", "MI_early", "MI_late", "HFpEF", "AngII"))%>%
  ggplot(., aes(x= source, y= Human.ref, fill= score))+
  geom_tile(color = "darkgrey")+
  scale_fill_gradientn(colors= c("white", "red"),
                                              breaks=c(0,2,4,6,8,10),#labels=c("Minimum",0.5,"Maximum"),
                                        limits=c(0,10))+
  #scale_fill_gradient(low= "white", high = "darkred")+
  geom_text(aes(label=sig))+
  theme_minimal()+
  theme(axis.text= element_text(color= "black"),
        panel.border = element_rect(colour = "black",fill=NA, size=1),
        axis.text.x=element_text(angle= 40, hjust= 1))+
  labs(x= "Murine IFS marker", y= "", fill = "-log10(q-value)")
p_s= unify_axis(p_s)

p.comb= cowplot::plot_grid(p_d+
                     theme( legend.position = "none",
                           plot.margin = unit(c(1, 0, 0, 0), "cm"),
                           axis.title.y = element_text(size= 9)),
                   p_s+theme( plot.margin = unit(c(1, 0, 0, 0), "cm"),
                              axis.text.y = element_blank(),
                              axis.title.x= element_text(margin = margin(t = 17))
                              ),
                   ncol =2,
                   align= "hv",
                   axis= "tb",
                   rel_widths   = c(1,1.8))
p.comb

pdf("output/figures/main/Fig4/bulk_heat.pdf",
        width= 6,
        height= 2.)
p.comb
dev.off()
# plot ORA results ----------------------------------------------------------------------------

p.complete= cowplot::plot_grid(p.signatures.hfref+
                                 coord_equal(), p.states.hfref+
                                 coord_equal(),
                   p.signatures.hfpef+
                     coord_equal(), p.marker.hfpef+
                     coord_equal(), ncol = 2,
                   rel_widths = c(1, 1.5)
                   )

pdf("output/figures/main/ora_res_reheat_hfpef.pdf",
    width= 7.5,
    height= 5)
p.complete
dev.off()



# check ranking -------------------------------------------------------------------------------

df.reheat= enframe(reheat) %>% mutate(r.v= rank(desc(value)))
map(gene_signatures$total, function(x){
  df.= df.reheat %>% filter(name %in%x)
  c(mean(df.$r.v), median(df.$r.v))

  })

# plot DE genes together -----------------------------------------------------------------

#bring logFCs into single df:
logFC.c= lapply(names(logFCs), function(x){
  rownames_to_column(logFCs[[x]], "gene") %>%
    mutate(study= x)
})%>% do.call(rbind, .)


## 1. plot HFPEF statistic against study wise logFC.

p.fc= logFC.c%>%
  left_join(dea%>% dplyr::select(-gene) %>% dplyr::rename(gene = MGI.symbol) , by= "gene")%>%
  as_tibble()
p.fc

de_genes= p.fc %>% filter(adj.P.Val< 0.05, t>0)%>% pull(gene) %>% unique()

p1= p.fc %>%
  mutate(topr= ifelse(gene %in% de_genes, "top500", "" ),
         label= ifelse(topr== "top500", gene, NA))%>%
  filter(study== "AngII" & gene %in% gene_signatures$total$AngII) %>%
  #drop_na %>%
  ggplot(., aes( x= logFC, y= avg_log2FC, col= topr))+
  geom_point()+
  geom_label_repel(aes(label = label))+
  theme_bw()+
  scale_color_manual(values= hfpef_cols)+
  theme(legend.position = "none")+
  labs(x= "log2 FC (Das et. al)",
       y= "avg_log2 FC (fibroblasts AngII")
p1
p2= p.fc %>%
  mutate(topr= ifelse(gene %in% de_genes, "top500", "" ),
         label= ifelse(topr== "top500", gene, NA))%>%
  filter(study== "HFpEF" & gene %in% gene_signatures$total$HFpEF) %>%
  #drop_na %>%
  ggplot(., aes( x= logFC, y= avg_log2FC, col= topr))+
  geom_point()+
  geom_label_repel(aes(label = label))+
  theme_bw()+
  scale_color_manual(values= hfpef_cols)+
  theme(legend.position = "none")+
  labs(x= "log2 FC (Das et. al)",
       y= "avg_log2 FC (fibroblasts HFpEF")
p2

p3= p.fc %>%
  mutate(topr= ifelse(gene %in% de_genes, "top500", "" ),
         label= ifelse(topr== "top500", gene, NA))%>%
  filter(study== "MI" & gene %in% gene_signatures$total$MI) %>%
  #drop_na %>%
  ggplot(., aes( x= logFC, y= avg_log2FC, col= topr))+
  geom_point()+
  #geom_label_repel(aes(label = label), max.overlaps = 100)+
  theme_bw()+
  scale_color_manual(values= hfpef_cols)+
  theme(legend.position = "none")+
  labs(x= "log2 FC (Das et. al)",
       y= "avg_log2 FC (fibroblasts MI")
p3

hfpef_consensus_plot= cowplot::plot_grid(p1,p2, p3,
                   nrow = 3,
                   labels= "AUTO")

pdf("output/figures/integration_studies/bulk_validation/gene_plots_hfpef.pdf",
    width= 5,
    height= 8)
hfpef_consensus_plot
dev.off()

## 2. plot ReheaT
re.= re. %>%as.data.frame() %>%  rownames_to_column("gene")

p.fc.r= logFC.c %>% left_join(re., by= "gene")

p1= p.fc.r %>%
  mutate(AngII= ifelse(gene %in% gene_signatures$total$AngII,"AngII",""),
         topr= ifelse(gene %in% re.[1:500,1], "top500", "" ),
         label= ifelse(AngII == "AngII" & topr== "top500", gene, NA))%>%
  filter(study== "AngII" & AngII == "AngII") %>%
  #drop_na %>%
  ggplot(., aes( x= V1, y= avg_log2FC, col= topr))+
  geom_point()+
  geom_label_repel(aes(label = label))+
  theme_bw()+
  scale_color_manual(values= hfpef_cols)+
  theme(legend.position = "none")+
  labs(x= "ReHeaT consensus statistic (directed)")
p1

intersect(de_genes,marker$`0`)
# plot volcanos -------------------------------------------------------------------------------

dea%>%
  #mutate(sig.hfpef= ifelse(MGI.symbol %in% marker$`0`, "y","no" ))%>%
  mutate(sig.hfpef= ifelse(MGI.symbol %in% gene_signatures$total$AngII, "y","no" ))%>%
  mutate(label= ifelse(sig.hfpef=="y", MGI.symbol, ""))%>%
  #filter(sig.hfpef=="y")%>%
  ggplot(.,(aes(x=logFC , y= -log10(adj.P.Val), color = sig.hfpef)))+
  geom_point(aes(alpha= 0.7))+
  ggrepel::geom_label_repel(mapping=aes(label=label),max.overlaps = 150)

p.fc= logFC.c%>%
  left_join(dea%>% dplyr::select(-gene) %>% dplyr::rename(gene = MGI.symbol) , by= "gene")%>%
  as_tibble()
p.fc

de_genes= p.fc %>% filter(adj.P.Val< 0.05, t>0)%>% pull(gene) %>% unique()

p1= p.fc %>%
  mutate(topr= ifelse(gene %in% de_genes, "top500", "" ),
         label= ifelse(topr== "top500", gene, NA))%>%
  filter(study== "AngII" & gene %in% gene_signatures$total$AngII) %>%
  #drop_na %>%
  ggplot(., aes( x= logFC, y= avg_log2FC, col= topr))+
  geom_point()+
  geom_label_repel(aes(label = label))+
  theme_bw()+
  scale_color_manual(values= hfpef_cols)+
  theme(legend.position = "none")+
  labs(x= "log2 FC (Das et. al)",
       y= "avg_log2 FC (fibroblasts AngII")
p1

# check correlation ---------------------------------------------------------------------------

x = p.fc%>%
  filter(study== "MI" & gene %in% gene_signatures$total$MI)
df1= cor.test(x$avg_log2FC, x$logFC)

x = p.fc%>%
  filter(study== "HFpEF" & gene %in% gene_signatures$total$HFpEF)
df2= cor.test(x$avg_log2FC, x$logFC)

x = p.fc%>%
  filter(study== "AngII" & gene %in% gene_signatures$total$AngII)
df3= cor.test(x$avg_log2FC, x$logFC)

df= data.frame(cor= c(df1$estimate, df2$estimate, df3$estimate),
               pval= c( df1$p.value, df2$p.value, df3$p.value),
               study= c("MI","HFpEF", "AngII"),
               bulk= rep("HFpEF", 3))

p.fc= logFC.c %>% left_join(re., by= "gene")

x = p.fc%>%
  filter(study== "MI" & gene %in% gene_signatures$total$MI)

df1= cor.test(x$avg_log2FC, x$V1)

x = p.fc%>%
  filter(study== "HFpEF" & gene %in% gene_signatures$total$HFpEF)
df2= cor.test(x$avg_log2FC, x$V1)

x = p.fc%>%
  filter(study== "AngII" & gene %in% gene_signatures$total$AngII)
df3= cor.test(x$avg_log2FC, x$V1)

df.2= data.frame(cor= c(df1$estimate, df2$estimate, df3$estimate),
               pval= c( df1$p.value, df2$p.value, df3$p.value),
               study= c("MI","HFpEF", "AngII"),
               bulk = rep("HFrEF", 3))
#df.2= column_to_rownames(df, "study")

rbind(df, df.2) %>%
  mutate(sig = ifelse(pval<0.05, "*", ""))%>%
  rename(rho= cor)%>%
  ggplot(., aes(x= study, y= bulk, fill = rho))+
  geom_tile()+
  scale_fill_gradient2(low= "darkblue", mid= "white", high = "red")+
  #geom_text(aes(label = sig))+
  geom_text(aes(label = paste0("p ",round(pval, 2))), size =3)+
  theme_bw()



# TF in bulk ----------------------------------------------------------------------------------

bulk$vooom$E
library(dorothea)
data("dorothea_hs")
net = dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C")) %>% dplyr::rename(source= tf)%>% mutate(likelihood =1)%>% distinct(source, target, mor, likelihood)

mat= bulk$DEA
mat= mat %>% select(gene, t)
mat= column_to_rownames(mat, "gene")

#mat= bulk$vooom$E
df = decouple(mat = mat, network = net, statistics = "ulm") %>%
  mutate(    p.adj = p.adjust(p_value, method = "BH"))

df %>% filter(source== "PPARA")%>% print(n=100)
top_tf= df %>% filter(p.adj<0.05)%>% arrange(desc((score)))%>% print(n=199)# %>% #slice(1:50)%>%
  pull(source)%>% unique()




# progeny in bulk -----------------------------------------------------------------------------
library(progeny)
library(ComplexHeatmap)
  gex= bulk$vooom$E
?run_progeny
M.Progeny = run_progeny(gex, .label = colnames(gex), organism = "Human")
hmap= Heatmap(t(M.Progeny), cluster_columns = T, name = "Progeny_score")
t.vec= bulk$DEA%>% select(t)%>% as.matrix()
rownames(t.vec)= bulk$DEA$gene
t.progeny = run_progeny(t.vec, .label = "t.s", organism = "Human")
hmap= Heatmap(t.progeny, cluster_columns = T, name = "Progeny_score")

naba= readRDS("data/prior_knowledge/NABA_mouse.rds")

naba$Basement_membranes
