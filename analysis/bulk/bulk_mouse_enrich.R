## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-12-07
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## enrich different sigs in samples#

library(tidyverse)
library(decoupleR)
library(cowplot)
library(edgeR)
# gex.objs
gex.obj= readRDS("output/bulk_mouse_obj.rds")

# separated
sep.obj= readRDS("output/bulk_mouse_fit.rds")


expr_plot = function(gene, gex.obj){
  enframe(gex.obj$exp[gene,], name="sample")%>%
    left_join(gex.obj$meta)%>%
    mutate(Group = factor(Group,
                          levels= c("Control",
                                    "10-weeks HFpEF",
                                    "15-weeks HFpEF")))%>%
    ggplot(., aes(x= Group, y= value))+
    #geom_jitter()+
    geom_boxplot(width= 0.5, alpha = 0.8)+
    #geom_jitter()+
    theme_bw()+
    theme(axis.text.x = element_text(angle= 60, hjust= 1),
          axis.text = element_text(size= 11, color ="black"),
          panel.border = element_rect(fill = NA, color ="black", size = 1))+
    ggtitle(gene)+
    labs(x= "", y= "normalized expression")
    #ggpubr::stat_compare_means(comparisons = list(c(1,2), c(2,3), c(1,3)))

  # enframe(gex.obj$exp[gene,], name="sample")%>%
  #   left_join(gex.obj$meta)%>%
  #   ggplot(., aes(x= Group2, y= value))+
  #   #geom_jitter()+
  #   geom_boxplot(width= 0.5, alpha = 0.8)+
  #   theme(axis.text.x = element_text(angle= 60, hjust= 1))+
  #   ggpubr::stat_compare_means(comparisons = list(c(1,2)))
  #
}
library(cowplot)
p1= plot_grid(
expr_plot(gene = "Angptl4", gex.obj),
expr_plot("Postn", gex.obj),
expr_plot(gene = "Col1a1", gex.obj),
expr_plot(gene = "Col1a2", gex.obj),
expr_plot(gene = "Col4a1", gex.obj),
expr_plot(gene = "Col4a2", gex.obj),
expr_plot(gene = "Col15a1", gex.obj),
expr_plot(gene = "Sparc", gex.obj),
expr_plot(gene = "Spp1", gex.obj)
)

p1

pdf("output/figures/supp/bulk_single_genes_exp.pdf",
    width= 4.5,
    height= 9)
p1
dev.off()


p2= plot_grid(
  expr_plot(gene = "Igfbp3", gex.obj),
  expr_plot("Pcsk6", gex.obj),
  expr_plot(gene = "Perp", gex.obj)
)

p2

genes= c("Angptl4", "Postn", "Col1a1", "Col4a1", "Col15a1", "Spp1")
wk10 = topTable(sep.obj$older$time,
                coef = "k10",
                adjust.method = "BH",number = Inf)
wk15 = topTable(sep.obj$older$time, coef = "wk15",adjust.method = "BH",number = Inf)
#wk15 = topTable(sep.obj$older$w15, coef = "wk15",adjust.method = "BH",number = Inf)
plot(wk10[genes,"t"] )
wk15[genes,]
wk10[genes,]
wk15[c("Igfbp3", "Pcsk6", "Perp"), ]
wk10[c("Igfbp3", "Pcsk6", "Perp"), ]

#genesets:
gene_sigs= readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
naba= readRDS("data/prior_knowledge/NABA_mouse.rds")
int.fib.marker= readRDS( "output/fib_integration/marker_list/integrated_marker.rds")

cut.off.state= 100
int.fib.marker= int.fib.marker%>% group_by(cluster) %>%
  filter(p_val_adj <0.05,
         avg_log2FC>0)%>%
  dplyr::slice(1:cut.off.state)
int.genesets= split(int.fib.marker$gene, int.fib.marker$cluster)


get.net= function(gset.list){
net= enframe(gset.list, name = "source", value= "target") %>%
  unnest(target) %>%
  mutate(mor= 1,
         likelihood= 1)

}

net= get.net(int.genesets)
net2= get.net(gene_sigs$total)
net3= get.net(naba)
# combine all genesets:
net= rbind(net,net2,net3)


# enrichment (decoupleR) -----------------------------------------------------------------------

## on the separated data by batch (exploration)
#on contrast log2fc:
mats= sep.obj$`5wk`

#names(mats)= gex.obj$de$gene
sep.decoupled= decouple(mat = as.matrix(mats), network = net, statistics = "fgsea")

p2= sep.decoupled%>%
  filter(statistic=="norm_fgsea")%>%
  ggplot(., aes(x=source , y= score, fill = -log10(p_value)))+
  geom_col()+
  coord_flip()+
  theme_bw()
p2


## on the t-stats from the 10 an 15 wk

mats= sep.obj$rest$hf$t
sep.decoupled= decouple(mat = as.matrix(mats), network = net, statistics = "fgsea")

p2= sep.decoupled%>%
  filter(statistic=="norm_fgsea")%>%
  ggplot(., aes(x=source , y= score, fill = -log10(p_value)))+
  geom_col()+
  coord_flip()+
  theme_bw()+
  ggtitle("hfpef_vs_ct")



mats= sep.obj$rest$time$t
sep.decoupled= decouple(mat = as.matrix(mats), network = net, statistics = "fgsea")

p2= sep.decoupled%>%
  filter(statistic=="norm_fgsea")%>%
  ggplot(., aes(x=source , y= score, shape = condition, color = -log10(p_value)))+
  geom_point()+
  #facet_grid(rows= vars(condition))+coord_flip()+
  theme_bw()+
  ggtitle("hfpef_vs_ct")+
  coord_flip()

p2

mats= cbind(sep.obj$older$w10$t, sep.obj$older$w15$t, sep.obj$older$hf$t)
fc= sep.obj$`5wk`[names(sep.obj$`5wk`) %in% rownames(mats)]
mats= mats[names(fc),]
mats= as.matrix(cbind(mats, fc))
colnames(mats)= c("week10_t", "week15_t", "week10+15_t", "week5_fc")

#enrich all
sep.decoupled= decouple(mat =mats,
                        network = net,
                        statistics = "ulm")

p2= sep.decoupled%>%
  filter(statistic=="norm_fgsea")%>%
  ggplot(., aes(x=source , y= score, shape = condition, color = -log10(p_value)))+
  geom_point()+
  #facet_grid(rows= vars(condition))+coord_flip()+
  theme_bw()+
  ggtitle("hfpef_vs_ct")+
  coord_flip()

#plot separate
names= list("NABA"= names(naba),
     "Fib_sig"= names(gene_sigs$total),
     "Fib_state"= names(int.genesets)
     )


p_list= map(names, function(x){
  #print(x)
  sep.decoupled%>%
    filter(statistic=="ulm",
           source %in% x)%>%
    mutate(sigs= ifelse(p_value<0.01, "**",
                        ifelse(p_value<0.05, "*", "ns")))%>%
    ggplot(., aes(x=source , y= score, shape = condition,
                  #color = -log10(p_value),
                  color = sigs
                  ))+
    geom_point(size= 4)+
    #facet_grid(rows= vars(condition))+coord_flip()+
    theme_bw()+
    ggtitle("")+
    coord_flip()+
    geom_hline(yintercept = 0, linetype= 3)+
    theme(axis.text = element_text(size= 11, color= "black"))
})

p_geneset_scores= cowplot::plot_grid(plotlist= p_list)

pdf("output/figures/supp/bulk_enrich_sets.pdf",
    width = 11,
    heigh = 8)
p_geneset_scores
dev.off()



mats

p_list= map(gene_sigs$total, function(x){

  res= apply(mats, 2, function(y){
    y2= sort(y, decreasing = T)

   names(y2[names(y2) %in% x])
    })
  res

  })

p_list= map(gene_sigs$total, function(x){

  res= apply(mats, 2, function(y){
    y2= sort(y, decreasing = T)

    enframe(y2)%>%
      arrange(desc(value))%>%
      mutate(gset= name %in% x,
             ranks = rank(desc(value)))%>%
      filter(gset)%>%
      ggplot(., aes(x=reorder(name, value),label= ranks,   y= value))+
      geom_point(size= 1.5)+
      geom_col(width= 0.1, color= "black")+
      coord_flip()+
      geom_text(color= "black", position = position_nudge(y = 0.5))

  })
  res
})




library(cowplot)
map(p_list, function(x){ plot_grid(plotlist = x)})
pdf("output/figures/supp/bulk_leading_edges_angII_hfpef_fib_sigs.pdf",
    width =11,
    height= 16
    )
plot_grid(plotlist = p_list$HFpEF)
plot_grid(plotlist = p_list$AngII)

dev.off()
pdf("output/figures/supp/bulk_leading_edges_mi_fib_sigs.pdf",
    width =11,
    height= 35
)
plot_grid(plotlist = p_list$MI_early)
plot_grid(plotlist = p_list$MI_late)

dev.off()



# plot 10 and 15 weeks ------------------------------------------------------------------------
mats= sep.obj$older$time$t
colnames(mats)= c("10weeks_HFpEF", "15weeks_HFpEF")

#remove the MI signature (whihc is now mi early and late)
net= net %>% filter(source != "MI")
#enrich all
sep.decoupled= decouple(mat =mats,
                        network = net,
                        statistics = "ulm")


p2= sep.decoupled%>%
  filter(statistic=="ulm")%>%
  ggplot(., aes(x=source , y= score, shape = condition, color = -log10(p_value)))+
  geom_point()+
  #facet_grid(rows= vars(condition))+coord_flip()+
  theme_bw()+
  ggtitle("hfpef_vs_ct")+
  coord_flip()
p2
#plot separate
names= list("NABA"= names(naba),
            "Fib_sig"= names(gene_sigs$total),
            "Fib_state"= names(int.genesets)
)


p_list= map(names, function(x){
  #print(x)
 df= sep.decoupled%>%
    filter(statistic=="ulm",
           source %in% x)%>%
    mutate(sigs= ifelse(p_value<0.01, "**",
                        ifelse(p_value
                               <0.05, "*", "ns")))

 #plot points:
   df %>%ggplot(., aes(x=source , y= score, shape = condition,
                  #color = -log10(p_value),
                  color = sigs
    ))+
    geom_point(size= 4)+
    #facet_grid(rows= vars(condition))+coord_flip()+
    theme_bw()+
    ggtitle("")+
    coord_flip()+
    geom_hline(yintercept = 0, linetype= 3)+
    theme(axis.text = element_text(size= 11, color= "black"),
          axis.text.x = element_text(angle= 45, hjust=1))


  df %>% ggplot(., aes(x=source , y= score, fill = condition,
                #color = -log10(p_value),
                #alpha = rev(sigs)
  ))+
    geom_bar(stat="identity", position=position_dodge(), color ="black")+
    #facet_grid(rows= vars(condition))+coord_flip()+
    theme_bw()+
    ggtitle("")+
    coord_flip()+
    geom_hline(yintercept = 0, linetype= 3)+
    theme(axis.text = element_text(size= 11, color= "black"),
          axis.text.x = element_text(angle= 45, hjust=1),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
})
pdf("output/figures/supp/bulk_fib_sigs1.pdf",
    width =5,
    height= 4
)
p_list[1]
dev.off()

pdf("output/figures/supp/bulk_fib_sigs2.pdf",
    width =5,
    height= 2
)
p_list[2]
dev.off()

pdf("output/figures/supp/bulk_fib_sigs3.pdf",
    width =5,
    height= 3
)
p_list[3]
dev.off()

### compare t-stats in bulk with eta.sq
df.aov2= readRDS( file = "output/aov.df.rds")

gset= gene_sigs$total$HFpEF

plot_xy= function(gset, df.aov){
  enframe(df.aov[rownames(df.aov) %in% gset,1]) %>%
    inner_join(enframe(mats[rownames(mats) %in% gset ,3]), by= "name")%>%
    ggplot(aes(x= value.x, y= value.y, label = name))+
    geom_point()+
    geom_text(position= position_nudge(y= 0.1))

}

plot_xy(gene_sigs$total$HFpEF, df.aov2$HFpEF)
plot_xy(gset = gene_sigs$total$AngII,df.aov =  df.aov2$AngII)
plot_xy(gset = gene_sigs$total$AngII,df.aov =  df.aov2$AngII)
plot_xy(gset = gene_sigs$total$MI_late,df.aov =  df.aov2$MI_l)
p_list$MI_late
sep.obj$rest$hf
gex.obj$de%>%
  filter(grepl("Ang", gene))



run_ulm_wrap= function(mats, net, target){
  mats= cpm(dge, log = T)
  x= decouple(mat = mats, network = net, statistics = "ulm")
  x= x%>%
    left_join(target%>%
                rename(condition= sample))
  p1= x%>%
    filter(statistic=="ulm",
           #source %in% net2$source)%>%
           source %in% names(int.genesets))%>%
    ggplot(., aes(x=Group , y= score, col = Group))+
    facet_grid(rows= vars(source), scales= "free")+
    #facet_grid(rows= vars(source), scales = "free")+
    geom_boxplot(outlier.fill = NA)+
    geom_jitter()
  p1
}

#mats= gex.obj$exp
#target= gex.obj$meta

mats= cpm(dge, log = T)
x= decouple(mat = mats, network = net, statistics = "ulm")
x= x%>%
  left_join(target%>%
              rename(condition= sample))
p1= x%>%
  filter(statistic=="ulm",
         #source %in% net2$source)%>%
         source %in% names(int.genesets))%>%
  ggplot(., aes(x=Group , y= score, col = Group))+
  facet_grid(rows= vars(source), scales= "free")+
  #facet_grid(rows= vars(source), scales = "free")+
  geom_boxplot(outlier.fill = NA)+
  geom_jitter()
p1

#on contrast:
mats= gex.obj$de$t
names(mats)= gex.obj$de$gene

x2= decouple(mat = as.matrix(mats), network = net, statistics = "ulm")
p2= x2%>%
  filter(statistic=="ulm")%>%
  ggplot(., aes(x=source , y= score, fill = -log10(p_value)))+
  geom_col()+
  coord_flip()+
  theme_bw()
p2
gex.obj$de%>%
  filter(grepl("Ang", gene))


fit2$t
mats= fit2$t
names(mats)= gex.obj$de$gene

x2= decouple(mat = as.matrix(mats), network = net, statistics = "ulm")
p2= x2%>%
  filter(statistic=="ulm",
         source %in% names(gene_sigs$total)
         )%>%
  ggplot(., aes(x=condition , y= score, fill = -log10(p_value)))+
  facet_grid(rows=vars(source))+
  geom_col()+
  coord_flip()+
  theme_bw()
p2
