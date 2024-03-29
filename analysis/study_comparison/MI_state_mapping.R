## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-03-22
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
## use forte et al MI study to see at which time point the integrated cluster pop up.
##

library(tidyverse)

source("code/utils.R")
source("code/utils_celltype_atlas.R")

integrated_data= readRDS("output/seu.objs/study_integrations/harmony_fib_filt.rds")
meta= integrated_data[[]]
rm(integrated_data)
meta
forte= readRDS( file = "../sc-exploration/output/cell_mi_files/fibro_subset_integrated.rds")
forte@meta.data= forte@meta.data %>%
  mutate(group = ifelse(OP== "myocardial infarction", "hf", "ct"),
         study= "forte")

DefaultAssay(forte)= "RNA"

rename_cells= function(seu, string_add_to_cellid){
  seu = RenameCells(seu,  string_add_to_cellid)
}
forte= rename_cells(forte, "b")
meta.mi = forte[[]]

df1= meta.mi %>% rownames_to_column("cellid")%>% as_tibble()
df2= meta%>% rownames_to_column("cellid")%>% as_tibble()

df= df1 %>% inner_join(df2, by= "cellid")

df= df%>%
  dplyr::select(cellid, opt_clust_integrated, time , group,  OP, orig.ident.x)
unique(df$OP)

prop.= calc_mean_proportions_per_group(df %>% filter(OP == "myocardial infarction"),
                                cluster.col = "opt_clust_integrated",
                                sample.col = "orig.ident.x",
                                group.col = "time")

prop.$groupwise%>%
  ggplot(., aes(factor(time), opt_clust_integrated, fill = mean.percent))+
  geom_tile()+
  scale_fill_gradient(low= "white", high = "darkred")+
  labs(fill= "mean proportion", y= "integrated cell states")+
  theme_minimal()

prop.$samplewise%>%
  ggplot(.,aes( color = factor(time), x= opt_clust_integrated, y = props))+
  geom_boxplot()




# add timepoint data to integrated object -----------------------------------------------------

df= df1 %>%
  select(cellid, time)%>%
  right_join(df2, by= "cellid")

table(df$time, df$group)

table(colnames(integrated_data) == df_sort$cellid)
table(colnames(integrated_data) %in% df$cellid)

df_sort = df[match( colnames(integrated_data), df$cellid),]

saveRDS(df_sort, "output/seu.objs/study_integrations/meta_with_time.rds")

integrated_data= AddMetaData(integrated_data, df_sort$time, col.name= "time")

saveRDS(integrated_data, "output/seu.objs/study_integrations/harmony_fib_filt_2.rds")
