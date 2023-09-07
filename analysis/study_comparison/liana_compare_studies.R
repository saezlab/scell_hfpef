## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2023-01-27
##
## Copyright (c) Jan D. Lanzer, 2023
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## load liana results per study and compare

library(tidyverse)
library(ComplexHeatmap)


#
#hfpef
hfpef= readRDS("output/liana_combined_results.rds")
angii= readRDS( "output/liana_angii.rds")
mie= readRDS( "output/liana_mi_e.rds")
mil= readRDS( "output/liana_mi_l.rds")

gene_intersect= readRDS("output/gene_intersect.rds")

li= list("HFpEF"= hfpef$mac.fibs,
     "AngII"= angii$mac.fibs,
     "MI_e"= mie$mac.fibs,
      "MI_l"= mil$mac.fibs)

df= map(names(li), function(x) (li[[x]]%>% mutate(study= x)))%>% do.call(rbind, .)%>%
  mutate("LR"= paste0(ligand, "->", receptor))%>%
  filter(ligand %in% gene_intersect & receptor %in% gene_intersect)


x1= df %>%
  group_by(study)%>%
  filter(hf_score>0.95)%>%
  #slice_max(order_by = hf_score, n= 30)%>%
  pull(LR)

df %>%
  filter(LR %in% x1)%>%
  ggplot(., aes(x= study, y= LR, fill = hf_score ))+
  geom_tile()

df.matrix= df %>%
  filter(LR %in% x1)%>%
  dplyr::select(LR, hf_score, study)%>%
  pivot_wider(names_from = LR, values_from = hf_score)%>%
  as.data.frame()%>%
  column_to_rownames("study")

Heatmap(t(df.matrix))

pdf("output/figures/liana_meta_hmap.pdf",
    height= 15,
    widt= 5)
Heatmap(t(df.matrix))
dev.off()

VlnPlot(forte_early, features = c("Spp1", "Tnf", "Itga9"),
        group.by = "celltype", split.by = "group",split.plot = T,
        fill.by= "group")
p1 <- VlnPlot(forte_early, c("Spp1"), pt.size = 0.1, split.by = "group", group.by = "celltype") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5)

unique(forte_late$group)
df %>% filter(grepl("MI", study))%>%
  pivot_wider(names_from = study, values_from=hf_score)%>%
  ggplot(., aes(x= MI_e, y= MI_l))+
  geom_point()


df %>%
  filter(LR %in% x1, grepl("Spp1", ligand))




# investiagte SPP1 ----------------------------------------------------------------------------
mie$liana_ct$natmi %>% filter(ligand == "Spp1")

mie$liana_ct$natmi %>% filter(ligand == "Spp1")

mie$mac.fibs%>% filter(ligand== "Spp1")
