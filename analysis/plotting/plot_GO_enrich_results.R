## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-02-28
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## plot GO term enrichments

library(tidyverse)

source("analysis/utils_funcomics.R")
go. = readRDS("output/funcomics_res/GO_cellstates.rds")

df= go.$fibs$bp

df= clean_names(df) %>%
  mutate(cluster= ifelse(cluster==0, "Col15a1+", cluster),
       cluster= ifelse(cluster==1, "Igfbp3+", cluster),
       cluster= ifelse(cluster==2, "Pi16+", cluster),
       cluster= ifelse(cluster==3, "Cxcl1+", cluster),
       cluster= ifelse(cluster==4, "Cilp+", cluster),
       cluster= ifelse(cluster==5, "Wif1+", cluster))


## barplot:
p.sigs= map(unique(df$cluster), function(x){
  p.df= df%>%
     #filter(#!Term %in% common_T,
    #        cluster ==,
    #        Adjusted.P.value<0.05)%>%
    # mutate(n= as.factor(n))%>%
      filter(cluster==x)%>%
     group_by(cluster)%>%
     top_n(10, -Adjusted.P.value)
  p.df%>%
    ggplot(aes(x= reorder(Term,-log10(Adjusted.P.value)),
               y=  -log10(Adjusted.P.value)
               #               fill = n
    )
    )+
    geom_col(fill= "#DF9B6D")+
    #facet_grid(vars(rows=cluster))+
    coord_flip()+
    labs(y="-log10(p.adj)", x= "")+
    theme_minimal()+
    geom_text(data = p.df,
              aes(label=Term),
              y = 0,
              hjust=0,
              vjust=0.75 )+
    theme(axis.text.y = element_blank())+
    ggtitle(label= x)
})


p.grid= cowplot::plot_grid(plotlist = p.sigs, nrow = 3)

pdf("output/figures/cell_type_analysis/fibros/cellstate_GO.pdf",
    width = 8,
    height= 9)
p.grid
dev.off()




# as heatmap ----------------------------------------------------------------------------------

##
T_count= df%>%
  filter(Adjusted.P.value<0.001)%>%
  select(Term, cluster)%>% group_by(Term) %>% count()

common_T= T_count%>%
  filter(n>4)%>% pull(Term)
specific_T= df %>%
  group_by(cluster)%>%
  filter(!Term %in% common_T)%>%
  top_n(5, wt = -Adjusted.P.value)%>%
  pull(Term)

df = df %>%
  left_join(T_count)


p.bp.common
p.bp.common= df%>%
  filter(Term %in% c( common_T, specific_T))%>%ungroup()%>%
  mutate(Term = factor(Term, levels = unique(c(common_T,specific_T))),
         cluster= factor(cluster, levels=c(df %>% group_by(cluster) %>% pull(cluster)%>% unique())
                         ),
         star = ifelse(Adjusted.P.value<0.001, "*", ""))%>%
  ggplot(aes(x= cluster, y= Term,  fill= -log10(Adjusted.P.value)))+
  geom_tile()+
  geom_text(aes(label= star))+
  scale_fill_gradient(low = "white", high = "darkred")+
  labs(x="", y= "", fill= "-log(q-value)")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle= 60, hjust = 1),
        axis.text=element_text(color= "black"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

p.bp.common

pdf("output/figures/funcomics/GO_fib_cellstate_heatmap.pdf",
    width = 6.5,
    height = 4.3)
p.bp.common
dev.off()
