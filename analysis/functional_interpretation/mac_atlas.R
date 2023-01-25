## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-10-07
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## intepret mac atlas


mac.maker=  c("Lyz2", "Cd68", "Cx3cr1")
integrated_data = readRDS( "output/seu.objs/study_integrations/macs/harmony_macs.rds")

Idents(integrated_data)

state.markers= FindAllMarkers(integrated_data, assay = "RNA")
saveRDS(state.markers, "output/statemarker_macs_integrated.rds")

x= state.markers%>% as_tibble()%>% filter(cluster==10)%>% pull(gene)
boom = x[1:10]

interest= c("Tnf", "Spp1")
FeaturePlot(integrated_data, features = interest)


# check composition ---------------------------------------------------------------------------


xy = integrated_data[[]]
x= xy %>%
  rownames_to_column("cellid") %>%
  group_by(group, opt_clust_integrated, study)%>%
  count

seu.meta= integrated_data[[]]

fibs.counts= seu.meta%>% group_by(study, group) %>% count()

x= seu.meta%>% group_by(opt_clust_integrated, group, study) %>%
  summarise(n = n())
x=x %>%ungroup() %>% group_by(study)%>%
  left_join(fibs.counts %>%
              rename(n.all= n),
            by=c("study", "group")) %>%
  mutate(freq2= n/n.all)

x %>% filter(study== "hfpef")%>% arrange(freq2)
x%>%filter(opt_clust_integrated== 5)
p.int.fib.prop= ggplot(x, aes(x= study, y= freq2, fill =opt_clust_integrated))+
  geom_col()+
  scale_fill_manual(values = col_vector)+
  theme_minimal()

x2= x %>% mutate(freak= sum(freq2)) %>% ungroup() %>% mutate(sacaled.prop= freq2/freak)
x2 %>% print(n=100)
