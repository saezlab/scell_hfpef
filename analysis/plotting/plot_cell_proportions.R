## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-10-15
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Plot cell proportions

library(Seurat)
library(tidyverse)
library(ggsignif)
source("code/utils.R")


# CELL TYPES ----------------------------------------------------------------------------------


# full data set
seu= readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")


x= seu[[]]
x$celltype2= str_replace_all(".1","", string= x$celltype2 )
x$celltype2= str_replace_all(".2","", string= x$celltype2 )
unique(x$celltype2)
#x = x %>% as_tibble() %>% ungroup() %>%  mutate(celltype3= base::ifelse(grepl("Fib", celltype2), "Fibroblasts", celltype2))
p.by.sample= calc_props(x, cluster.col = "orig.ident", "celltype")
p.by.sample2= calc_props(x, cluster.col = "orig.ident", "celltype2")
p.by.group= calc_props(x, cluster.col = "group", "celltype")
p.by.group2= calc_props(x, cluster.col = "group", "celltype2")

df.mean= calc_mean_proportions_per_group(x, "celltype2", "group")

df. = df.mean$groupwise%>%
  mutate(celltype2  = factor(celltype2,levels = c("Fibroblasts",
                                                    "Endothelial",
                                                    "NK.cells",
                                                    "Macrophages",
                                                    "CD8.T.cells",
                                                "CD4.T.cells",
                                                    "B.cells",
                                                    "Granulocytes",
                                                    "SMC/Pericytes")))
p.mean= plot_mean_proportions( df., "celltype2", "Celltype Composition")

p.mean. = p.mean+
  geom_signif(
    y_position = c(0.6, 0.25, 0.1,0.18),
    xmin = c(0.8, 1.8, 3.8, 6.8),
    xmax = c(1.2, 2.2, 4.2, 7.2),
    annotation = rep("*", 4), tip_length = 0
  ) +
  ylim(c(0,0.7))+
  labs(y= "mean %",
       x= "",
       fill= "")+
  ggtitle("")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(color= "black"))

p.mean.
#save figure for compositional plottig
saveRDS(p.mean., "output/figures/main/proportions_obj.rds")


pdf(file = "output/figures/proportion_per_group_final.pdf",
    width= 6,
    height= 4)
p.mean.
dev.off()


pdf(file = "output/figures/proportion_per_group.pdf",
    width= 6,
    height= 5)
p.by.sample
p.by.sample2
p.by.group
p.by.group2
p.mean

dev.off()



# CELL STATES ---------------------------------------------------------------------------------

seu.objs=  readRDS("output/seu.objs/cell.state.obj.list.rds")

##fibs=
seu= seu.objs$fibroblasts

y= seu[[]]
p.mean= calc_mean_proportions_per_group(y, "cellstate", "group")
p.mean= plot_mean_proportions(p.mean$groupwise, "cellstate", "Cellstate Composition")

p.fibs=p.mean +
  labs(y= "mean %",
       x= "",
       fill= "")+
  geom_signif(
    y_position = c(0.4, 0.2, 0.15), xmin = c(0.8, 3.8, 5.8), xmax = c(1.2, 4.2, 6.2),
    annotation = rep("*", 3), tip_length = 0
  ) +
  ylim(c(0,0.5))+
  ggtitle("")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text = element_text(color= "black"))

saveRDS(p.fibs, "output/figures/main/Fig2/props.rds")

pdf("output/figures/main/cellstate_fibs_proportions.pdf",
    width = 5,
    height= 4)
p.fibs
dev.off()

##macs
seu= seu.objs$macrophages
y= seu[[]]

p.mean= calc_mean_proportions_per_group(y, "cellstate", "group")
p.mean= plot_mean_proportions(p.mean$groupwise, "cellstate", "Cellstate Composition")

pmacs= p.mean +
  labs(y= "mean %",
       x= "",
       fill= "")+
  # geom_signif(
  #   y_position = c(0.28, 0.28, 0.42, 0.4, 0.22),
  #   xmin = c(0.8, 1.8,2.8, 3.8,4.8),
  #   xmax = c(1.2,2.2, 3.2, 4.2, 5.2),
  #   annotation = rep("*",5), tip_length = 0
  # ) +
  ylim(c(0,0.5))+
  ggtitle("")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(color= "black"))
pmacs


pdf("output/figures/main/cellstate_macs_proportions.pdf",
    width = 5,
    height= 4)
pmacs
dev.off()


# plot heatmap of proportion change between studies -------------------------------------------

p.val.df=readRDS( "output/fib_integration/p.vals.proportional2.rds")
int.fibs= readRDS("output/seu.objs/study_integrations/harmony_fib_filt.rds")
meta_seu= int.fibs@meta.data
rm(int.fibs)
meta_seu = readRDS("output/seu.objs/study_integrations/meta_with_time_and_doubleCT.rds")
## prop fc calc

df= map(unique(meta_seu$study), function(x){
  print(x)
  p.mean= meta_seu %>% filter(study==x)

  p.mean= calc_mean_proportions_per_group(p.mean, "opt_clust_integrated", "group")
  return(p.mean)
})

df.fc= lapply(df, function(x){
x$groupwise %>%pivot_wider(names_from=group , values_from= mean.percent, -median.percent)%>%
  mutate(fc = hf-ct)

})

names(df.fc)= c("AngII", "HFpEF",  "MI_early", "MI_late" )

df.fc= map(names(df.fc), function(x){
  df.fc[[x]]%>% mutate(study= x)
})%>% do.call(rbind,.)

#add pval
p.val.df = map(names(p.val.df), function(x){
  p.val.df[[x]]%>% mutate(study= x)
})%>% do.call(rbind,.)%>%
  dplyr::rename(opt_clust_integrated= int_cellstate)
p.val.df= p.val.df%>%
  mutate(study= ifelse(study== "circ", "AngII",
                       ifelse(study== "forte", "MI",
                              ifelse(study== "hfpef", "HFpEF", study))),
         sig= ifelse(p.val<0.01, "**",
                     ifelse(p.val<0.1, "*","")))

p.fc.proportion= df.fc %>%
  left_join(p.val.df, by= c("study", "opt_clust_integrated"))%>%
mutate(study= factor(study, levels=c("HFpEF", "MI_early", "MI_late", "AngII")))%>%
  ggplot(., aes(x= study, y= opt_clust_integrated, fill= fc))+
  geom_tile()+
  #geom_text(mapping=aes(label= round(p.val,2))  )+
  geom_text(mapping=aes(label= sig))+
  scale_fill_gradient2(low= "blue", mid= "white",
                      high= "red", midpoint = 0)+
  theme_minimal()+
  theme(axis.text = element_text(color= "black"),
        axis.text.x=element_text(angle= 40, hjust= 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  coord_equal()+
    labs(x= "", y= "integrated cell state",
         fill= "composition \n difference ")
p.fc.proportion


pdf("output/figures/main/Fig3/proportion_change_fc.pdf",
    width =4.6,
    height= 3)
p.fc.proportion
dev.off()


ggplot(df.fc, aes(x= study, y= ))
