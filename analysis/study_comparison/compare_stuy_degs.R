## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-10-27
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## compare gene signature sets considering the main direction of regulation within a dataset


library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(ggvenn)
#
# hfpef.z= readRDS("output/figures/integration_studies/hfpef.z.transform.rds")
# forte.z= readRDS("output/figures/integration_studies/forte.z.transform.rds")
# circ.z= readRDS("output/figures/integration_studies/circ.z.transform.rds")
objs= readRDS( "output/fib_integration/marker_list/intstate07_singlestudyDEG.rds")
objs= readRDS( "output/fib_integration/marker_list/intstate0_singlestudyDEG.rds")

gene_signatures= readRDS( "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled2.rds")
gene_signatures = readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
gene_signatures$total= gene_signatures$total[names(gene_signatures$total)!="MI" ]

logFCs=readRDS(file = "output/fib_integration/marker_list/DEG_per_study_LOGFC.rds")
#names(logFCs)=c("HFpEF", "MI", "AngII")

# 1st) plot DE genes together -----------------------------------------------------------------
#x= logFCs$hfpef

logFC.c= lapply(names(logFCs), function(x){
  rownames_to_column(logFCs[[x]], "gene") %>%
    mutate(study= x)
})%>% do.call(rbind, .)

#prepare for saving
df= logFC.c %>%
  filter(gene %in% unlist(gene_signatures$total)) %>%
  select(-pct.1, -pct.2)%>%
  pivot_wider(names_from = study, values_from= avg_log2FC)%>%
  rename_at(vars(-contains("gene")), .funs =  ~paste0(., "_logFC"))%>%
  mutate(HFpEF = ifelse(gene %in% gene_signatures$total$HFpEF, "HFpEF", ""),
         MI = ifelse(gene %in% gene_signatures$total$MI, "MI", ""),
         AngII = ifelse(gene %in% gene_signatures$total$AngII, "AngII", ""))

WriteXLS(df, ExcelFileName = "output/fib_integration/marker_list/deg_fibs_study_wise.xlsx")

#x= gene_signatures$unique$hfpef

library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("darkblue", "white", "darkred"))
col_fun(seq(-3, 3))

Hmaps.U= map(gene_signatures$unique, function(x){
  df= logFC.c %>% filter(gene %in% x) %>% select(-pct.1, -pct.2)%>%
    pivot_wider(names_from = study, values_from= avg_log2FC)
  df= column_to_rownames(df, "gene")
  Heatmap(df[,], col= col_fun, name= "avg_log2FC",
          cluster_columns = F,
          column_names_rot = 40,
          border = T,
          #rect_gp = gpar(ol = "black", lty = 1),
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10))
})

Hmaps.O= map(gene_signatures$overlap, function(x){
  df= logFC.c %>% filter(gene %in% x) %>% select(-pct.1, -pct.2)%>%
    pivot_wider(names_from = study, values_from= avg_log2FC)
  df= column_to_rownames(df, "gene")
  Heatmap(df, col= col_fun, name= "avg_log2FC",   cluster_columns = F,
          column_names_rot = 40,
          border = T,
          #rect_gp = gpar(ol = "black", lty = 1),
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10))
})

pdf("output/figures/main/Fig3/HFPEF.sig.pdf",
    width= 3,
    height= 8)
Hmaps.U[1]
dev.off()

pdf("output/figures/main/Fig3/Ang.sig.pdf",
    width= 3,
    height= 4.6)
Hmaps.U[3]
dev.off()

pdf("output/figures/main/Fig3/Fibro_all.pdf",
    width= 3,
    height= 3)
Hmaps.O[1]
dev.off()

x= grid.grabExpr(draw(Hmaps.U[[1]], padding = unit(c(0, 0, 0, 0), "mm")))
x1= grid.grabExpr(draw(Hmaps.U[[3]], padding = unit(c(0, 0, 0, 0), "mm")))
x2= grid.grabExpr(draw(Hmaps.O[[1]], padding = unit(c(0, 0, 0, 0), "mm")))
plot_grid(x,x1,x2, nrow= 1)
# save heatmaps -------------------------------------------------------------------------------


DEAs= map(names(objs), function(x) (objs[[x]]%>% mutate(study=x)%>% rownames_to_column("gene")))
names(DEAs)= names(objs)
DEAs$hfpef

get_genes= map(DEAs, function(x) (x %>%filter(p_val_adj<0.01,
                                              avg_log2FC>0.4) %>% pull(gene) ))




upr.venn= ggvenn(get_genes)

df.comb= do.call(rbind, DEAs)

unique_hfpef = get_genes$hfpef[!get_genes$hfpef %in% c(get_genes$forte, get_genes$circ)]

whole_fib_dea= readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET.rds")

#Q1 is the upregulation of the hfpfef sig in clust0 specific?
df.comb%>%
  mutate(star= ifelse(p_val_adj<0.01, "*", ""))%>%
  group_by(study) %>%
  filter(gene %in% unlist(whole_fib_dea$unique$hfpef))%>%
  ggplot(aes(x= study,y= reorder(gene, avg_log2FC), fill= avg_log2FC))+
  geom_tile()+
  scale_fill_gradient2(low="blue", mid="white", high= "red")+
  theme(axis.text = element_text(size= 13))+
  ggtitle("HFpEF Fibro genes")+
  geom_text(aes(label=star ))

df.comb%>%
  mutate(star= ifelse(p_val_adj<0.01, "*", ""))%>%
  group_by(study) %>%
  filter(gene %in% sample(unlist(whole_fib_dea$unique$forte), 50))%>%
  ggplot(aes(x= study,y= reorder(gene, avg_log2FC), fill= avg_log2FC))+
  geom_tile()+
  scale_fill_gradient2(low="blue", mid="white", high= "red")+
  theme(axis.text = element_text(size= 11))+
  ggtitle("circ fibrosis genes")+
  geom_text(aes(label=star ))


df.comb%>%
  mutate(star= ifelse(p_val_adj<0.01, "*", ""))%>%
  group_by(study) %>%
  filter(gene %in% unlist(whole_fib_dea$overlap$all))%>%
  ggplot(aes(x= study,y= reorder(gene, avg_log2FC), fill= avg_log2FC))+
  geom_tile()+
  scale_fill_gradient2(low="blue", mid="white", high= "red")+
  theme(axis.text = element_text(size= 13))+
  ggtitle("Common fibrosis signature")+
  geom_text(aes(label=star ))


all.hfpef.genes= df.comb%>% group_by(study) %>%
  filter(gene %in% unlist(get_genes$hfpef))%>%
  ggplot(aes(x= study,y= reorder(gene, avg_log2FC), fill= avg_log2FC))+
  geom_tile()+
  scale_fill_gradient2(low="blue", mid="white", high= "red")+
  theme(axis.text = element_text(size= 13))

unique.hfpef.genes= df.comb%>%
  group_by(study) %>%
  filter(gene %in% unique_hfpef)%>%
  mutate(labels = ifelse(gene %in% whole_fib_dea$unique$hfpef, "darkgreen", "black"))%>%
  ggplot(aes(x= study,y= reorder(gene, avg_log2FC), fill= avg_log2FC))+
  geom_tile()+
  scale_fill_gradient2(low="blue", mid="white", high= "red")+
  ggtitle("unique upregulated HFpEF")+
  theme(axis.text = element_text(size= 13))


unique.hfpef.genes.cell.level= df.comb%>%
  group_by(study) %>%
  filter(gene %in% unique_hfpef,
         !gene %in% whole_fib_dea$unique$hfpef)%>%
  ggplot(aes(x= study,y= reorder(gene, avg_log2FC), fill= avg_log2FC))+
  geom_tile()+
  scale_fill_gradient2(low="blue", mid="white", high= "red")+
  ggtitle("unique upregulated HFpEF")+
  theme(axis.text = element_text(size= 13))

unique.hfpef.genes.cell.level

hfpef.focus= df.comb%>%
  filter(p_val_adj <0.0001,
         study=="hfpef",
         abs(avg_log2FC)>0.3)%>%
  pull(gene)

df.comb %>%
  filter(gene %in% hfpef.focus)%>%
  ggplot(aes(x= study,y= reorder(gene, avg_log2FC), fill= avg_log2FC))+
  geom_tile()+
  scale_fill_gradient2(low="blue", mid="white", high= "red")+
  ggtitle("unique upregulated HFpEF")+
  theme(axis.text = element_text(size= 13))




# use complex heatmap for simultaneous clustering:

 matrix.wide= df.comb%>% group_by(study) %>%
  filter(gene %in% unique_hfpef) %>%
  select(study, gene, avg_log2FC) %>% pivot_wider(names_from = study, values_from= avg_log2FC)%>%
   arrange(desc(hfpef))%>%
   as.data.frame() %>% column_to_rownames("gene")
 #remove NAs:
 matrix.wide=  replace(matrix.wide, is.na(matrix.wide), 0)


library(ComplexHeatmap)
library(pheatmap)
pheatmap(matrix.wide)
Heatmap(matrix.wide)
Heatmap(matrix.wide[rownames(matrix.wide) %in% whole_fib_dea$unique$hfpef,],clustering_method_rows = "ward.D2"
        )
x1= dist(matrix.wide)
c= hclust(x1, method = "ward.D2")
plot(c)
genes_oi= DEAs$hfpef%>% filter(p_val_adj<0.05) %>% arrange(desc(abs(avg_log2FC))) %>%
  #slice(1:50)%>%
  pull(gene)


df.comb%>% group_by(study) %>%
  filter(gene %in% unlist(genes_oi))%>% mutate(star= ifelse(p_val_adj<0.05, "*", ""))%>%
  ggplot(aes(x= study, y= reorder(gene, avg_log2FC), fill= avg_log2FC))+
  geom_tile()+
  scale_fill_gradient2(low="blue", mid="white", high= "red")+
  geom_text(mapping = aes(label=star))

### look for genes from fibs in hfpef:
fib_de=  readRDS("output/cell_specific_obj/cell_type_based_DE/cells_DEA.rds")
genes_oi= fib_de$fibroblasts %>% filter(p_val_adj<0.05)%>% pull(gene)
gene_DE= df.comb %>%
  filter(gene %in% genes_oi) %>% select(study, gene, avg_log2FC) %>% pivot_wider(names_from = study, values_from= avg_log2FC)
#%>%  filter(hfpef >0 & forte <=0 & circ <= 0)%>%pull(gene)

genes_hfpef.de= df.comb%>% group_by(study) %>%
  filter(gene %in% unlist(gene_DE))%>% mutate(star= ifelse(p_val_adj<0.05, "*", ""))%>%
  ggplot(aes(x= study, y= reorder(gene,avg_log2FC),  fill= avg_log2FC))+
  geom_tile()+
  scale_fill_gradient2(low="blue", mid="white", high= "red")+
  geom_text(mapping = aes(label=star))+
  ggtitle("DEG in all fibroblasts hfpef vs. ct _compared")

get_genes

##use DE_genes from fibroblast analysis

# run DE analysis on z-transformed data


# plot results: -------------------------------------------------------------------------------

pdf(file = "output/figures/integration_studies/fibroblast.hfpef.marker.pdf",
    height= 12,
    width= 8)
upr.venn
all.hfpef.genes
unique.hfpef.genes
genes_hfpef.de
dev.off()



# compare ppara regulons:  --------------------------------------------------------------------
library(dorothea)
data("dorothea_mm")

plot_regulon= function(tf, logFC.c){
  require(circlize)
  tf_reg= dorothea_mm%>%
    dplyr::filter(tf==!!tf) %>%
    pull(target)

  df= logFC.c %>% filter(gene %in% tf_reg) %>% select(-pct.1, -pct.2)%>%
    pivot_wider(names_from = study, values_from= avg_log2FC)

  df= column_to_rownames(df, "gene")
  df= df %>%arrange(desc(HFpEF))
  tf

  col_fun = colorRamp2(c(min(df, na.rm = T), 0, max(df, na.rm = T)), c("darkblue", "white", "darkred"))

  Heatmap(df, col= col_fun, name= "avg_log2FC", cluster_rows = T)

}

Smad3= plot_regulon("Smad3", logFC.c)
PParg = plot_regulon("Pparg", logFC.c)
PPara= plot_regulon("Ppara", logFC.c)
HIf= plot_regulon("Hif1a", logFC.c)

ppara=dorothea_mm%>% filter(tf== "Ppara") %>%pull(target)
print(s)

library(progeny)
?progeny
# compare single study_hfpef genes ------------------------------------------------------------

obsj=readRDS( file = "output/fib_integration/marker_list/DEG_per_study_in_fibs.rds")

genesoi=   map(objs, function(x){
  x=x%>% filter(p_val_adj<0.05, avg_log2FC>0)
  rownames(x)})

genes.all.fibs= genesoi$hfpef[!genesoi$hfpef %in% genesoi$forte]
genes.clust.fibs= unique_hfpef[!unique_hfpef %in% genesoi$forte]



# correlate direction of study DEG ------------------------------------------------------------

#merge logFCs to one df
logFCs2= map(names(logFCs), function(y){
  logFCs[[y]]%>% dplyr::rename(!!y := avg_log2FC)%>% rownames_to_column("gene")
})

fc.df= logFCs2[[1]] %>%
  full_join(logFCs2[[2]], by= "gene") %>%
  full_join(logFCs2[[3]], by= "gene") %>%
  full_join(logFCs2[[4]], by= "gene")

#heatmap

#create dfs to store correlation estimates
studies= names(gene_signatures$total)
comb.studies= c("MI.AngII", "HFpEF.AngII", "HFpEF.MI")

cor.df <- data.frame(matrix(ncol = 4, nrow = 4))
rownames(cor.df)=colnames(cor.df) = studies
#colnames(cor.df)=comb.studies
cor.p= cor.df

# loop to perform correlations for unique gene sets
corr.frame= lapply(studies, function(sig){

  df= fc.df %>% filter(gene %in% gene_signatures$total[[sig]])

  for(s.x in studies){
    for(s.y in studies){
      print(s.x)
      print(s.y)
      res= cor.test(df[[s.x]],df[[s.y]])
      cor.p[s.x, s.y] = res$p.value
      cor.df[s.x, s.y] = res$estimate
    }
  }
  return(list(cor.df, cor.p))
})

names(corr.frame)= studies

pls = map(corr.frame, function(x){

  # tidy the data for plotting
  cor.df= x[[1]] %>% rownames_to_column("signature") %>% pivot_longer(-signature,names_to= "comparison" ,
                                                                      values_to = "rho")

  cor.p= x[[2]] %>% rownames_to_column("signature") %>% pivot_longer(-signature,names_to= "comparison" ,
                                                                    values_to = "p.val")

  p.correlation.comp= cor.df %>% left_join(cor.p) %>%
    mutate(sig= ifelse(p.val<0.001, "*", ""),
           comparison = factor(comparison, levels= c("HFpEF", "AngII", "MI_early", "MI_late")),
           signature= factor(signature, levels= c("HFpEF", "AngII", "MI_early", "MI_late")))  %>%
    ggplot(., aes(x= signature, y = comparison, fill = rho))+
    geom_tile(color= "black")+
    scale_fill_gradient2(low= "darkblue", mid= "white", high= "darkred", midpoint = 0, limits= c(-1, 1))+
    #scale_fill_continuous(c(0,1))+
    theme_minimal()+
    geom_text(aes(label= sig))+
    theme(axis.text= element_text(colour = "black"))+
    labs(x= "",
         y= "")+
    coord_equal()


})

cowplot::plot_grid(plotlist = pls)

res= map(names(corr.frame), function(x){
  #x="HFpEF"
  # tidy the data for plotting
  cor.df= corr.frame[[x]][[1]] %>% rownames_to_column("signature") %>% pivot_longer(-signature,names_to= "comparison" ,
                                                                      values_to = "rho")%>%
    filter(signature== x)

  cor.p=corr.frame[[x]][[2]] %>% rownames_to_column("signature") %>% pivot_longer(-signature,names_to= "comparison" ,
                                                                     values_to = "p.val")%>%
    filter(signature== x)

  cor.df %>% left_join(cor.p)
})%>% do.call(rbind, .)

res

p.correlation.comp= res %>%
  mutate(p.val= p.adjust(p.val, "BH"))%>%
  mutate(sig= ifelse(p.val<0.001, "**",
                     ifelse(p.val<0.01, "**",
                            ifelse(p.val<0.05, "*", ""))),
         sig = ifelse(signature== comparison, "", sig),
         rho = ifelse(signature== comparison, NA, rho),
         comparison = factor(comparison, levels= c("HFpEF", "AngII", "MI_early", "MI_late")),
         signature= factor(signature, levels= c("HFpEF", "AngII", "MI_early", "MI_late"))) %>%
  ggplot(., aes(x= signature, y = comparison, fill = rho))+
  geom_tile(color= "black")+
  scale_fill_gradient(low= "white", high= "darkred")+
  theme_minimal()+
  geom_text(aes(label= sig))+
  theme(axis.text= element_text(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  coord_equal()

pdf(file = "output/figures/main/Fig2/deg.corr.pdf",
    width= 3,
    height= 3)
unify_axis(p.correlation.comp)+
  theme(axis.text.x = element_text(angle= 45, hjust= 1))
dev.off()
#mi.ang
  res= cor.test(df$MI_late,df$AngII )
  cor.df[sig, 1]= res$estimate
  cor.p[sig, 1]= res$p.value

  #hfpef.ang
  res= cor.test(df$HFpEF,df$AngII )
  cor.df[sig, 2]= res$estimate
  cor.p[sig, 2]= res$p.value

  #hfpef, mi
  res= cor.test(df$HFpEF,df$MI )
  cor.df[sig, 3]= res$estimate
  cor.p[sig, 3]= res$p.value


}

# tidy the data for plotting
cor.df= cor.df %>% rownames_to_column("signature") %>% pivot_longer(-signature,names_to= "comparison" ,
                                                            values_to = "rho")

cor.p= cor.p %>% rownames_to_column("signature") %>% pivot_longer(-signature,names_to= "comparison" ,
                                                            values_to = "p.val")

p.correlation.comp= cor.df %>% left_join(cor.p) %>%
  mutate(sig= ifelse(p.val<0.05, "*", ""),
         signature= factor(signature, levels= c("HFpEF", "AngII", "MI")))  %>%
  ggplot(., aes(x= signature, y = comparison, fill = rho))+
  geom_tile()+
  scale_fill_gradient2(low= "darkblue", mid= "white", high= "darkred")+
  theme_minimal()+
  geom_text(aes(label= sig))+
  theme(axis.text= element_text(colour = "black"))

pdf("output/figures/integration_studies/deg.hmaps/deg.correlation.pdf",
    width = 3.3 ,
    height= 2)
p.correlation.comp
dev.off()

# SANDBOX -------------------------------------------------------------------------------------


gene.means= map(names(z.seu), function(y){

  gex1= GetAssayData(z.seu[[y]], assay = "RNA", slot= "count")
  gene.mean.z = apply(gex1, 1, function(x){mean(x, na.rm= T)})
  enframe(gene.mean.z, name= "gene", value= "mean.z") %>% mutate(study= y)

})



##plot mean z-scores
z.seu= list("forte"= forte.z$seu,
            "circ"= circ.z$seu,
            "hfpef"= hfpef.z$seu)

gene.means= map(names(z.seu), function(y){

  gex1= GetAssayData(z.seu[[y]], assay = "RNA", slot= "count")
  gene.mean.z = apply(gex1, 1, function(x){mean(x, na.rm= T)})
  enframe(gene.mean.z, name= "gene", value= "mean.z") %>% mutate(study= y)

})

gene.means

df= do.call(rbind, gene.means)

p.wide= df %>% pivot_wider( values_from = mean.z, names_from= study)
remove_inf <- function(x) (replace(x, is.infinite(x), max(x)) )

p.wide= p.wide %>% mutate_at(c("circ", "hfpef", "forte"), remove_inf)

p.wide %>% arrange(desc(hfpef))

p.wide= p.wide %>% mutate(sums= abs(forte)+abs(circ)+abs(hfpef)) %>% arrange(desc(sums))
