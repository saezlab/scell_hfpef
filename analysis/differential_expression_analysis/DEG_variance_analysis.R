## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-05-14
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## test for set of DEGs whether rather cell state marker or general response.
#

library(tidyverse)
library(Seurat)
library(lsr)
library(ComplexHeatmap)
library(ggrepel)
library(cowplot)
library(pROC)
library(rstatix)
library(ggpubr)


source("code/utils.R")

fib = readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
fib$total= fib$total[names(fib$total)!="MI" ]
fib$total= fib$total[c(1,2,4,3)]
degs= unique(unlist(fib$total))
degs= degs[!str_detect(degs, "-")]
# seu= readRDS("output/seu.objs/cell.state.obj.list.rds")
# seu= seu$fibroblasts
# degs= fib$total$hfpef
#
# circ= readRDS( file = "../sc-exploration/output/circ_obj/fibroblasts_seu.rds")
#
# DefaultAssay(circ)= "RNA"
# circ@meta.data= circ@meta.data %>% mutate(group = ifelse(group == "angII", "hf", "ct" ))
# circ@meta.data%>% select(treatment, group) %>% table()
# seu= circ

## load fib atlas

int.seu = readRDS("output/seu.objs/study_integrations/harmony_fib_filt_2.rds")

# funcs ---------------------------------------------------------------------------------------

prep_input= function(seu, group, clust, degs){

  #get normalized expr data
  log.data= GetAssayData(seu, slot = "data", assay = "RNA")
  # get cell state labels:
  Idents(seu)= clust
  labels.state= CellsByIdentities(seu, )
  labels.state= enframe(labels.state, name="state", value= "cellid") %>% unnest(cellid)

  Idents(seu) = group
  labels.group= CellsByIdentities(seu)
  labels.group= enframe(labels.group, name="group", value= "cellid" )%>% unnest(cellid)

  #dim(as.matrix(log.data[,]))
  degs= degs[degs %in% rownames(log.data)]

  #prepare dataframe to run anovas on
  model.df= (log.data[degs,])%>%as.data.frame()%>%
    t() %>%
    as.data.frame()%>%
    rownames_to_column("cellid")%>%
    as_tibble()%>%
    left_join(labels.state, by= "cellid")%>%
    left_join(labels.group, by= "cellid")
}

do.aov = function(model.df, degs){
  df.aov= sapply(degs, function(x){
    print(x)

    #rename "-"

    x= str_replace_all(x, "-", "_")
    colnames(model.df)= str_replace_all(colnames(model.df), "-", "_")

    if(!x %in% colnames(model.df)){
      return(NULL)
    }

    aov.1=aov(as.formula(paste0(x ,"~ state")), data = model.df)
    aov.2=aov(as.formula(paste0(x ,"~ group")), data = model.df)
    eta1= etaSquared(aov.1)[1,1]
    eta2= etaSquared(aov.2)[1,1]

    c(eta1, eta2, eta1*eta2)


  })
  rownames(df.aov)= c("eta.state", "eta.group", "sg.ratio")
  df.aov= t(df.aov)
}



# angII ---------------------------------------------------------------------------------------

cl= fib$total$AngII[fib$total$AngII != "Hbb-bs"]
unique(int.seu$study)

circ= subset(int.seu, subset = study=="circ")

colnames(circ@meta.data)

df.circ= prep_input(subset(int.seu, subset = study=="circ"),
                    group =  "group",
                    clust = "opt_clust_integrated",
                    degs = degs)
rm(circ)
aov.circ= do.aov(df.circ, degs)


Heatmap(aov.circ[,1])

p1= fgsea::plotEnrichment(pathway= fib$overlap$all, stats= sort(aov.circ[,1]))
x= fgsea::fgseaSimple(pathway= list("gset "= fib$overlap$all), stats= sort(aov.circ[,1]), nperm= 100)

top_exp = sort(aov.circ[,1])

top_comp = rev(sort(aov.circ[,1]))

VlnPlot(circ, features= names(top_exp[1:10]), split.by = "group")
VlnPlot(circ, features= names(top_comp[1:10]),  split.by = "group",assay= "RNA")

VlnPlot(circ, features= c("Col1a1", "Meox1"), split.by = "group")
VlnPlot(circ, features= names(top_comp[1:10]), assay= "RNA")
VlnPlot(circ, features= names(top_exp[1:6]), assay= "RNA")

plot(aov.circ[,1], y= aov.circ[,2])


# mi ------------------------------------------------------------------------------------------

##late
mi_l= subset(int.seu, subset = study=="forte" & (time== 14 | group == "ct"))

df.mi_l= prep_input(mi_l ,group =  "group", clust = "opt_clust_integrated", degs =degs)
rm(mi_l)
aov.mi_l= do.aov(df.mi_l, degs)

p1= fgsea::plotEnrichment(pathway= fib$overlap$all, stats= aov.mi_l[,1])

top_exp = sort(aov.mi_l[,1])
top_comp = rev(sort(aov.mi_l[,1]))

VlnPlot(mi_l, features= names(top_comp[1:10]), assay= "RNA", split.by = "group")
VlnPlot(mi_l, features= names(top_exp[1:10]), assay= "RNA", split.by = "group")


##early
mi_e= subset(int.seu, subset = study=="forte" & time!= 14)
df.mi= prep_input(mi_e ,group =  "group", clust = "opt_clust_integrated", degs = degs)
rm(mi_e)
aov.mi= do.aov(df.mi, degs)

p1= fgsea::plotEnrichment(pathway= fib$overlap$all, stats= aov.mi[,1])

top_exp = sort(aov.mi[,1])
top_comp = rev(sort(aov.mi[,1]))

VlnPlot(mi_e, features= names(top_comp[1:10]), assay= "RNA", split.by = "group")
VlnPlot(mi_e, features= names(top_exp[1:10]), assay= "RNA", split.by = "group")




# hfpef----------------------------------------------------------

hfpef= subset(int.seu, subset = study=="hfpef")

df.hfpef= prep_input(hfpef ,group =  "group", clust = "opt_clust_integrated", degs = degs)
rm(hfpef)
aov.hfpef= do.aov(df.hfpef,degs)


# union ---------------------------------------------------------------------------------------

df.= prep_input(int.seu ,group =  "group", clust = "opt_clust_integrated", degs = degs)

aov.= do.aov(df.,degs)
aov.=cbind(aov., "sg_ratio"= aov.[,"eta.state"] / aov.[, "eta.group"])

Heatmap(aov.[fib$total$HFpEF,1])

aov. %>%
  as.data.frame() %>%
  rownames_to_column("gene")%>%
  as_tibble()%>% arrange(eta.state)%>%
  filter(gene %in% fib$total$AngII)%>%
  print(n=100)


aov.tibble= aov. %>%
  as.data.frame() %>%
  rownames_to_column("gene")%>%
  as_tibble()%>% arrange(eta.state)%>%
  filter(gene %in% fib$total$MI_early)%>%
  print(n=100)



# get AUROCs for single genes -----------------------------------------------------------------

genes= c("Col1a1", "Col1a2", "Col4a1", "Col4a2", "Postn", "Angptl4")
genes = fib$total$HFpEF
genes = degs

cell.lists= readRDS("output/cell_names_study.rds")
df= int.seu@assays$RNA@data[genes,]
res=
  lapply(genes, function(gene){
  sapply(names(cell.lists), function(x){
    df= map(unique(int.seu@meta.data$opt_clust_integrated), function(y){
      # print(gene)
      # print(y)
      # print(x)
      cells= cell.lists[[x]]
      meta.x=  int.seu@meta.data[cells, ] %>% filter(  opt_clust_integrated==y)
      resp=meta.x%>% pull(group)
      cells.oi= meta.x %>% rownames()
      d= auc(response= resp,predictor =df[gene,cells.oi])
      return(as.numeric(d[1]))
    })
    names(df)= unique(int.seu@meta.data$opt_clust_integrated)
    return(df)
  })

} )
names(res) = genes

saveRDS(list("aurocs"= res,
             "aov"= aov.tibble),
        "output/aov.and.aurocs.all.rds")

res2 = readRDS("output/aov.and.aurocs.all.rds")
res= res2$aurocs
aov.tibble= res2$aov
c("HFpEF", "AngII", "MI_early", "MI_late")
#as.matrix(res[[gene]])
#pheatmap::pheatmap(as.matrix(res[[gene]]))
pls= map(genes, function(gene){

  df= as.data.frame(res[[gene]])%>% mutate_all(as.numeric)
  colnames(df)= c("HFpEF", "AngII", "MI_early", "MI_late")
  df= df%>%
    rownames_to_column("state")%>%
    pivot_longer(-state, names_to= "study",
                 values_to= "auroc")#%>%
  stat.test= df  %>%
    mutate(study= factor(study)) %>%
  wilcox_test(auroc ~ study, paired = FALSE)
  val= 0.1/(dim(stat.test)[1] -1)
  p= ggboxplot(df, x = "study", y = "auroc") +
    stat_pvalue_manual(
      stat.test, label = "p.adj.signif",

      y.position = seq(0.9:1, by = val)
    )+
    geom_point(aes(color= state))+
    scale_color_manual(values=col_vector)+
    ggtitle(gene)+
    theme(axis.text.x = element_text(angle= 45, hjust= 1))
  # ggplot(., aes(x= name, y= value))+
  #   geom_boxplot()+
  #   scale_color_manual(values = col_vector)+
  #   geom_point(aes(color= state))+
  #   ggtitle(gene)
  # #%>% Heatmap(., cluster_rows =  F, cluster_columns = F,
   #                                                               name = gene)
  #Heatmap(as.matrix(res[[gene]]), col = col_vector)
unify_axis(p)
  })

p1= plot_grid(plotlist = pls, nrow= 1)

pdf("output/figures/general_hfpef_sig.pdf",
    width= 15,
    height= 5)
p1
dev.off()
### calculate median AUROC for each gene
names(fib$total)= c("hfpef" ,"ang"  , "mi_e"  ,"mi_l" )

pls= map(colnames(res$Angptl4), function(study){
  med.aurocs= map(res, function(x){
    median(unlist(x[,study]))
  }) %>% enframe(name = "gene",
                 value= "median_AUROC")%>%
    mutate(median_AUROC= unlist(median_AUROC))

  df.joined= aov.tibble %>%
    right_join(med.aurocs)

  df.joined %>%
    #filter(gene %in% fib$total[[study]])%>%
    mutate(stig= gene %in% fib$total[[study]])%>%
    ggplot(aes(x= eta.state, y= median_AUROC))+
   #geom_point(color= "darkgrey")+
    geom_point(aes(color=stig))+
  #geom_label_repel(mapping=aes(label = gene),color = "black", size = 3, max.overlaps = 30)+
  theme_minimal()+
  theme(legend.position = "none")+
    ylim(c(0.4, 0.9))+
    xlim(c(0, 0.45))+
    ggtitle(study)
})

plot_grid(plotlist = pls)
pls[1]
pls[2]
# combine with state dependecy



# compare genes  ------------------------------------------------------------------------------
## save the lists:

mat= cbind("circ"= do.circ, "mi_e"= do.mi_e,"mi_l"=  do.mi_l,"hfpef"= do.hfpef)
boxplot(mat)
x= list(aov.hfpef,aov.circ, aov.mi, aov.mi_l)
names(x)= c("HFpEF", "AngII", "MI_e", "MI_l")

saveRDS(x, file = "output/aov.df.rds")
x= readRDS( file = "output/aov.df.rds")

## add s / g ratios
ratios= map(x, function(y){
  print(head(y))
  #table(fib$overlap$all %in% rownames(y))
  y=cbind(y, "sg_ratio"= y[,"eta.state"] / y[, "eta.group"])
  y
})

### plot eta state for hfpef:
vec= ratios$HFpEF[fib$total$HFpEF,1]
pdf(file = "output/figures/main/Fig4/eta_state_hfpef.pdf",
    width= 2.3, height= 11)
Heatmap( vec, name= "eta² state", color= "black")
dev.off()


y= map(names(ratios), function(x){
  as.data.frame(ratios[[x]])%>% rownames_to_column("gene")%>%
    mutate(study= x)%>%
    as_tibble()
})

# plot s/g ratio for all genes

p.all.genes= do.call(rbind, y)%>%
  mutate(study= factor(study, levels= c("HFpEF", "AngII", "MI_early", "MI_late")))%>%
  ggplot(., aes(x= study, y= sg_ratio))+
  geom_hline(yintercept = 1)+
  geom_violin()+
  geom_boxplot(width = 0.3, outlier.colour = NA)+
  scale_y_log10()+
  theme_bw()+
  theme(panel.border = element_rect(size= 1) )+
  labs(y= "eta state / eta group")+
  stat_compare_means(comparisons = list(c(1,2), c(2,3), c(1,3) , c(3,4)))

y2= map(names(ratios), function(x){
  as.data.frame(ratios[[x]])%>%
    rownames_to_column("gene")%>%
    filter(gene %in% gene_signatures$total[[x]])%>%
    mutate(study= x)%>%
    as_tibble()
})
library(ggpubr)
p2= do.call(rbind, y2)%>%
  mutate(study= factor(study, levels= c("HFpEF", "AngII", "MI_early", "MI_late")))%>%
  ggplot(., aes(x= study, y= sg_ratio))+
  geom_hline(yintercept = 1)+
  geom_violin()+
  geom_boxplot(width = 0.3, outlier.colour = NA)+
  scale_y_log10()+
  theme_bw()+
  theme(panel.border = element_rect(size= 1) )+
  labs(y= "eta state / eta group")+
  stat_compare_means(comparisons = list(c(1,2), c(2,3), c(1,3) , c(3,4)))
p2

pdf("output/figures/main/Fig4/s_g_ratio_quantified.pdf",
    width = 4,
    height= 4)
unify_axis(p.all.genes)
unify_axis(p2)
dev.off()

#only for  signatures



###### compare dropout rates ( 0 elements)
calc_drop_outs= function(df, genes){
  vec= c()
  for (gene in genes){
    t= prop.table(table(df[[gene]]==0))
    vec= c(vec, t[2])
  }
  names(vec)= genes
  return(vec)
}
do.circ= calc_drop_outs(df.circ, degs)
do.mi_l= calc_drop_outs(df.mi, degs)
do.mi_e= calc_drop_outs(df.mi_l, degs)
do.hfpef= calc_drop_outs(df.hfpef, degs)




####### plot eta squared in comparison between studies
plot_points= function(df.aov, gset= c() ,titles){
  plot.df= rownames_to_column(as.data.frame(df.aov), "gene")
  p1= plot.df %>%
    mutate(signature= gene %in% fib$total[[titles]])%>%
    mutate(col= ifelse(gene %in% gset, "y", "n"))%>%
    mutate(label= ifelse(gene %in% gset, gene, ""))%>%
    ggplot(., aes(x= eta.state, eta.group, label = label, col = signature))+
    geom_abline(slope= 1, intercept= 0, color = "darkgrey", lty = 3)+
    geom_point()+
    scale_color_manual(values= hfpef_cols)+
    geom_label_repel(mapping=aes(label = label),color = "black", size = 3, max.overlaps = 100)+
    theme_minimal()+
    theme(legend.position = "none",
          panel.border = element_rect(colour = "black", fill=NA, size=1))+
    ggtitle(titles)
  unify_axis(p1)

}

degs = unique(rownames((x$HFpEF)))

gset= c("Angptl4", "Postn", "Cilp", "Meox1", "Sparc", "Col1a1", "Col1a2", "Col4a1", "Col4a2", "Pi16",
        "Dkk3", "Cthrc1", "Actb", "Acta2", "Mmp1", "Mmp2", "Timp1", "Loxl1", "Loxl2")

names(fib$total)= c( "HFpEF"   , "AngII"  ,  "MI_early" ,"MI_late" )
names(ratios)=  c( "HFpEF"   , "AngII"  ,  "MI_early" ,"MI_late" )

pls= map2(names(ratios), ratios, function(x, df.aov){
  plot_points(df.aov, gset, x)
})

points= cowplot::plot_grid(plotlist = pls)

points

pdf("output/figures/supp/anova.sigs.pdf",
    width = 7.5,
    height= 7.5)
points
dev.off()



#gset= Reduce(intersect, fib$total)


#look up genes

#check collagen expr:
p= sapply(ratios, function(x){
  #x[c("Col1a1", "Col1a2"),4]
  x[fib$overlap$all,4]
  #x[c("Loxl1", "Loxl2"),4]
  #x[gset,4]
})

p.collagens= t(p)%>%as.data.frame()%>%
  rownames_to_column(.,"disease_model")%>%
  pivot_longer(names_to="gene", values_to = "sg.ratio", -disease_model)%>%
  mutate(disease_model= factor(disease_model, levels=c("HFpEF", "AngII", "MI_early", "MI_late")))%>%
  ggplot(aes(color= gene, y= sg.ratio, x= disease_model, group=gene))+
  scale_y_log10()+
  #geom_boxplot(alpha= 0.2, color= "black")+
  #scale_color_manual(values = (col_vector))+
  geom_hline(yintercept = 1, color= "black")+
  geom_point(size=3, alpha= 0.9)+
  geom_line(alpha= 0.5)+
  theme_minimal()+
  labs(x= "", y= "eta² state / eta² group")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle= 40, hjust= 1))
p.collagens


pdf("output/figures/main/Fig4/collagen_compare_sg.pdf",
    width= 4,
    height= 4)
unify_axis(p.collagens)
dev.off()

    comp.top = names(sort(ratios$HFpEF[fib$total$HFpEF, 1], decreasing = T)[1:10])
exp.top = names(sort(ratios$HFpEF[fib$total$HFpEF, 1], decreasing = F)[1:10])
p= sapply(ratios, function(x){
  x[c("Col1a1", "Col1a2", "Col4a1", "Col4a2"),4]
  #x[c("Angptl4", "Ech1","Hsbp1", ),4]
  #x[exp.top,4]
  #x[gset,4]
})

p.collagens4= t(p)%>%as.data.frame()%>%
  rownames_to_column(.,"disease_model")%>%
  pivot_longer(names_to="gene", values_to = "sg.ratio", -disease_model)%>%
  mutate(disease_model= factor(disease_model, levels=c("HFpEF", "AngII", "MI_early", "MI_late")))%>%
  ggplot(aes(color= gene, y= sg.ratio, x= disease_model))+
  geom_point(size=4)+
  scale_y_log10()+
  scale_color_manual(values = (col_vector))+
  geom_hline(yintercept = 1, color = "darkgrey")+
  theme_minimal()+
  labs(x= "", y= "eta² state / eta² group")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle= 40, hjust= 1))
unify_axis(p.collagens4)

pdf("output/figures/main/Fig4/collagen4_compare_sg.pdf",
    width= 3,
    height= 4)
unify_axis(p.collagens4)
dev.off()


## check ratio difference for hfpef and angii
ratios$HFpEF[fib$total$HFpEF,]%>% as.data.frame()%>%
  rownames_to_column("gene")%>% as_tibble()%>% arrange(sg_ratio)
new.vec= ratios$HFpEF[fib$total$HFpEF,4] / ratios$AngII[fib$total$HFpEF,4]
Heatmap(sort(log10(new.vec)))
  #check collagen expr:
p= sapply(ratios, function(x){
  x[rownames(x)%in% naba$Collagens,4]
  #x[naba$Collagens,4]
})

t(p)%>%
  as.data.frame()%>% rownames_to_column(.,"disease_model")%>%
  pivot_longer(names_to="gene", values_to = "sg.ratio", -disease_model)%>%
  ggplot(aes(x= gene, y= log10(sg.ratio), color= disease_model))+
  geom_point()+
  geom_hline(yintercept = 0)+
  theme_minimal()
# plot for the dis sigs:



# plot for the dis sigs:
names(ratios)
names(fib$total)
fib$total= fib$total[c(1,2,4,3)]
# x= ratios$HFpEF
# y= fib$total$HFpEF

p_full_sig= map2(ratios, fib$total, function(x,y){
  print(head(x))
    return(x[rownames(x) %in% y,4])
})%>% enframe()%>% #
  unnest()%>%
  mutate(name= factor(name, levels = c("HFpEF", "AngII", "MI_e", "MI_l")))%>%
  #filter(value<200)%>%
  ggplot(aes(x= name, y= value))+
  geom_jitter(alpha= 0.4)+
  geom_boxplot(width= 0.4)+
  #geom_hline(yintercept = 1, color= "darkgrey")+
  scale_y_log10()+
  labs(x= "",
       y= "log10(eta.state/eta.group)")+
  theme_minimal()

p.sigs= unify_axis(p_full_sig)
p.sigs

ratios$HFpEF%>% as.data.frame() %>% rownames_to_column("gene")%>% as_tibble()%>% arrange(desc(eta.state))
pdf("output/figures/main/Fig4/hfpef.vec.pdf",
    height= 11,
    width= 2.2)
print(Heatmap(sort(ratios$HFpEF[fib$total$HFpEF,1], decreasing = T), cluster_rows=F ,
        name = "Eta² cell states"
        )
      )
dev.off()



# plot single genes exp -----------------------------------------------------------------------

sort(ratios$HFpEF[fib$total$HFpEF,1], decreasing = T)[1:5]
p1= VlnPlot(hfpef, features = names(sort(ratios$HFpEF[fib$total$HFpEF,1], decreasing = T))[1:5],
            ncol = 5,
            cols = col_vector,
            combine = F)

remove_axis_labs= function(pls){

  pls= map(pls, function(p){
    p+
    labs(x="", y= "")+
    theme(legend.position = "none")
    #unify_axis(p)
    })
  cowplot::plot_grid(plotlist = pls, nrow= 1)
}

p1.m= remove_axis_labs(p1)

p2= VlnPlot(hfpef, features = names(sort(ratios$HFpEF[fib$total$HFpEF,1], decreasing = F))[1:5],
            ncol = 5, cols = col_vector,
            combine = F)#, split.by =  "group")

p2.m= remove_axis_labs(p2)

p3= cowplot::plot_grid(p1.m, p2.m, ncol = 2)

pdf("output/figures/main/Fig4/top_hfpef_genes.pdf",
    height= 3,
    width= 15)
p1.m
p2.m
dev.off()

#plot for the union
p_full= map(ratios, function(x){
  print(head(x))
  return(x[,4])
})%>% enframe()%>% #
  unnest()%>%
  mutate(name= factor(name, levels = c("HFpEF", "AngII", "MI_e", "MI_l")))%>%
  #filter(value<200)%>%
  ggplot(aes(x= name, y= value))+
  geom_jitter(alpha= 0.4)+
  geom_boxplot(width= 0.4)+
  geom_hline(yintercept = 1, color= "darkgrey")+
  scale_y_log10()+
  labs(x= "",
       y= "log10(eta.state/eta.group)")+
  theme_minimal()

p_full= unify_axis(p_full)

dis_sigs= plot_grid(p_full, p.sigs)

pdf("output/figures/supp/anova.sigs_boxpls.pdf",
    width = 7,
    height= 7)
dis_sigs
dev.off()

# plot naba boxplots --------------------------------------------------------------------------
fib$total
x= map(names(naba), function(x){

  cm = degs[degs %in% naba[[x]] ]

  p_Mat= map(ratios, function(x){
    x[cm,4]
  })%>% enframe()%>% #
    unnest()%>%
    #filter(value<200)%>%
    ggplot(aes(x= name, y= value))+
    geom_jitter(alpha= 0.4)+
    geom_boxplot(width= 0.4)+
    geom_hline(yintercept = 1, color= "darkgrey")+
    scale_y_log10()+
    labs(x= "",
         y= "log10(eta.state/eta.group)")+
    theme_minimal()+
    ggtitle(x)
  unify_axis(p_Mat)
})
p_nabas = cowplot::plot_grid(plotlist = x[c(2,4)])

pdf("output/figures/supp/anova.sigs_boxpls_naba.pdf",
    width = 11,
    height= 7)
p_nabas
dev.off()

#try to filter for degs
x= map(names(naba), function(x){

  cm = degs[degs %in% naba[[x]] ]

  p_Mat= map2(ratios,fib$total, function(x,y){
    cm = cm[cm %in% y ]
    x[cm,4]
  })%>% enframe()%>% #
    unnest()%>%
    #filter(value<200)%>%
    ggplot(aes(x= name, y= value))+
    geom_jitter(alpha= 0.4)+
    geom_boxplot(width= 0.4)+
    geom_hline(yintercept = 1, color= "darkgrey")+
    scale_y_log10()+
    labs(x= "",
         y= "log10(eta.state/eta.group)")+
    theme_minimal()+
    ggtitle(x)
  unify_axis(p_Mat)
})
p_nabas = cowplot::plot_grid(plotlist = x[c(2,7,5,1)], ncol = 1)
p_nabas

pdf("output/figures/supp/anova.sigs_boxpls_naba_select.pdf",
    width = 2.7,
    height= 7.5)
p_nabas
dev.off()

##plot as hmap


x= lapply(names(naba), function(x){

  cm = degs[degs %in% naba[[x]] ]

  p_Mat= map(ratios, function(x){
    median(x[cm,4])
  })%>% enframe()%>% #
    unnest()%>%
    mutate(naba =x)

})
do.call(rbind, x)%>%
  ggplot(., aes(x = name, y= naba, fill = log10(value)))+
  geom_tile()+
  scale_fill_gradient(low= "white", high = "red")
##


df= map(ratios, function(x){
  x[cm,4]
})

sort(df$HFpEF)
df%>%
  enframe()%>% #
  unnest()%>%
  #filter(value<200)%>%
  ggplot(aes(x= name, y= log10(value)))+
  geom_jitter(alpha= 0.4)+
  geom_boxplot(width= 0.4)+
  geom_hline(yintercept = log10(1))+
  #scale_y_log10()+
  labs(x= "Disease Model",
       y= "log10(state.e2/group.e2)")

x= map(df, function(x){
  res= t.test(x,mu = 1)
  res$statistic
})

enframe(x)%>% unnest(value)%>%
  ggplot(aes(x= name, y= value))+
  geom_col()
names(x)



naba= readRDS("data/prior_knowledge/NABA_mouse.rds")


# plot gene exp -------------------------------------------------------------------------------
features= c("Col1a1", "Col1a2")

p.gs= function(features, int.seu){
  p1= VlnPlot(subset(int.seu, study== "hfpef" ), features = features, split.by = "group", combine = F)%>%
    remove_axis_labs()

  p2= VlnPlot(subset(int.seu, study== "circ" ), features = features, split.by = "group", combine= F)%>%
    remove_axis_labs()
  p3=  VlnPlot(subset(int.seu, study== "forte" ), features = features, split.by = "group", combine= F)%>%
    remove_axis_labs()

  cowplot::plot_grid(plotlist= list(p1, p2,p3), nrow = 3)

}
p1= p.gs(features= c("Col1a1", "Col1a2"), int.seu)
p2= p.gs(features= c("Col4a1", "Col4a2"), int.fibs)
p3= p.gs(features= c("Mmp1", "Mmp2"), int.fibs)

pdf("output/figures/supp/collagen.plots.pdf",
    width= 4,
    height= 8)
p1
p2
dev.off()

map(names(naba), function(x){
  VlnPlot(subset(int.seu, study== "hfpef" ), features = features)
  VlnPlot(subset(int.seu, study== "circ" ), features =  naba$Collagens)
  VlnPlot(subset(int.seu, study== "forte" ), features =  naba$Collagens, group.by = "group")#, split.by= "study")
})

# old -----------------------------------------------------------------------------------------


library(fgsea)
##load genesets:

## use metabolic network from COSMOS to annotate genes involved in metabolism
load("data/prior_knowledge/meta_network.RData")
genes= unique(c(meta_network$source, meta_network$target))

View(meta_network)
gene_translate= readRDS("data/prior_knowledge/gene_translate.rds")

metabs <- meta_network[grepl("Metab__HMDB", meta_network$source) | grepl("Metab__HMDB", meta_network$target),]
metabs <- metabs[,-2]
metabs <- metabs[grepl("Gene", metabs$source) | grepl("Gene", metabs$target),]

metabs_reactant <- metabs[grepl("Metab__HMDB",metabs$source),]
metabs_products <- metabs[grepl("Metab__HMDB",metabs$target),]

names(metabs_reactant) <- c("metab","gene")
names(metabs_products) <- c("gene","metab")

metabs <- as.data.frame(rbind(metabs_reactant, metabs_products))
metabs$gene <- gsub("Gene.*__","",metabs$gene)
metabs$metab <- gsub("_[a-z]$","",metabs$metab)
metabs$metab <- gsub("Metab__","",metabs$metab)
str_to_title(metabs$gene)

m.genes= gene_translate %>% filter(Gene.name %in% metabs$gene)%>% pull(MGI.symbol)

fgseaSimple(pathways= list("fibrosis_core"= fib$overlap$all,
                           "BM"= naba$Basement_membranes) ,stats= df.aov[,1], nperm= 1000 )
p1= fgsea::plotEnrichment(pathway= fib$overlap$all, stats= df.aov[,1])
p2= fgsea::plotEnrichment(pathway= naba$Basement_membranes, stats= df.aov[,1])
p3= fgsea::plotEnrichment(pathway=str_to_title(metabs$gene), stats= df.aov[,1])
df.aov[rownames(df.aov) %in% str_to_title(metabs$gene),]

df.aov[m.genes,1]
cowplot::plot_grid(p1+ggtitle("shared fibrosis"),
                   p2+ggtitle("BM"))


topgroup
naba= readRDS("data/prior_knowledge/NABA_mouse.rds")
fgseaSimple(pathways= naba ,stats= df.aov[,1], nperm= 1000 )

msigDB= readRDS("data/prior_knowledge/msigDB.mouse_translated.rds")
msigDB$MSIGDB_KEGG

fgsea.res= fgseaSimple(pathways= list("glucose"= msigDB$MSIGDB_REACTOME$REACTOME_GLUCOSE_METABOLISM,
                                      "FA"= msigDB$MSIGDB_KEGG$KEGG_FATTY_ACID_METABOLISM),
                       stats= df.aov[,1], nperm= 2000 )

fgsea.res= fgseaSimple(pathways= msigDB$MSIGDB_REACTOME,
                       stats= df.aov[,1], nperm= 2000 )
fgsea.res%>% as_tibble()%>% arrange(NES, padj)



# metabolism annotation -----------------------------------------------------------------------


