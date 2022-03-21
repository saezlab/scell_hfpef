## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-10-15
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
## Differential cell proportion analysis
## Get null distributions by label randomization and calculate empirical p-values.

library(Seurat)
library(tidyverse)

# full data set
seu= readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")

source("analysis/utils.R")

meta_seu= seu[[]]


# function ------------------------------------------------------------------------------------

get_null= function(t= 500, x, clust_col= "cellstate", group_col= "orig.ident"){
  set.seed(4)

  res= lapply(c(1:t), function(i){
    colname2= paste0("rep",i)
    x1=  sample(x[[group_col]], replace = F)
    x[[colname2]]= x1
    prop=calc_mean_proportions_per_group(x, clust_col, colname2, colname2)
    prop$samplewise %>%
      mutate(group= str_replace_all(.data[[colname2]], "[:digit:]", ""))%>%
      group_by(.data[[clust_col]],group)%>%
      mutate(m.p = mean(props)) %>%
      distinct(.data[[clust_col]], group,m.p)%>%
      pivot_wider(names_from = group, values_from= m.p)%>%
      mutate(delta.p = hf-ct)
  }) %>%
    do.call(rbind, .)
}

get_empirical_p= function(real_props, nulls, clust_col= "cellstate"){

  cells= unique(real_props[[clust_col]])

  res= sapply(cells, function(x){

    print(x)
    null.cell= nulls%>% filter(.data[[clust_col]]==x)

    null.cell= null.cell %>%
      summarise(mean= mean(delta.p),
                sd= sd(delta.p))

    effect.size= real_props%>% filter(.data[[clust_col]]==x)%>% pull(delta.p)
    print(effect.size)
    #
    # p_emp_inc= length(null.cell$delta.p[null.cell$delta.p>effect.size]) /
    #   length((null.cell$delta.p))
    #
    # p_emp_dec= length(null.cell$delta.p[null.cell$delta.p<effect.size]) /
    #   length((null.cell$delta.p))
    #
    # p= c(p_emp_dec, p_emp_inc)
    # names(p)= c("dec", "inc")
    #
    # return(p)
    #
    # prop.table(table(null.cell$delta.p>effect.size))[2]
    #
    # p_emp_dec= prop.table(table(null.cell$delta.p<effect.size))

    p_1= pnorm(q = abs(effect.size),
              mean = null.cell$mean,
              sd= null.cell$sd,
              lower.tail = F
              )

    p_2= pnorm(q = abs(effect.size),
              mean = null.cell$mean,
              sd= null.cell$sd,
              lower.tail = T
              )
    print(c(p_1, p_2))
    return(min(c(p_1, p_2)))
  })

  res
  # names(res)= cells
  # enframe(res, name = "celltype", value= "p.val") %>%
  #   mutate(p.val= unlist(p.val))
}

get_null_visuals= function(real_prop, nulls){

  cells= unique(real_prop$celltype_mod)

  res= map(cells, function(x){

    null.cell= nulls%>% filter(celltype_mod==x)
    # null.cell= null.cell %>%
    #   summarise(mean= mean(delta.p),
    #             sd= sd(delta.p))
    #
    effect.size= real_prop%>% filter(celltype_mod==x)%>% pull(delta.p)

    ggplot(null.cell,aes(x= delta.p))+
      geom_histogram(bins=100)+
      geom_vline(xintercept = effect.size
      )+
      ggtitle(x)+
      theme_classic()

  })

  cowplot::plot_grid(plotlist = res)
}

test_normality= function(nulls){
  celltype = unique(nulls$celltype_mod)
  ps = map(celltype, function(x){
    df= nulls%>% filter(celltype_mod==x)
    test.r= shapiro.test(df$delta.p)
    test.r$p.value
  })
  if(all(ps>0.05)){
    message("all nulls passed")
  }else{
    message("some distributions are not normal")
    return(ps)
    }
}

# CELL TYPE ----------------------------------------
# keep T.cells separated (CD4, CD8) but test Fibroblasts together
meta_seu = meta_seu %>% mutate(celltype_mod= ifelse(grepl("Fib", celltype2), "Fibroblasts", celltype2))
nulls= get_null(t= 1000,x=  meta_seu[c("orig.ident", "celltype_mod")],clust_col =  "celltype_mod")

## get real observed proportions
real_prop= calc_mean_proportions_per_group(meta_seu, "celltype_mod", "group")

real_prop2 =real_prop$groupwise%>%
  distinct(celltype_mod, group, mean.percent)%>%
  pivot_wider(names_from = group, values_from=  mean.percent)%>%
  mutate(delta.p = hfpef-ct)

## calculate p-val per cell type
test_normality(nulls)
ps = get_empirical_p(real_prop = real_prop2 ,nulls =  nulls,clust_col =  "celltype_mod")

##plot =

cell_prop_ps= enframe(ps, name = "celltype", value = "p.val")%>%
  ggplot(., aes(x= celltype, y= -log10(p.val)))+
  geom_hline(yintercept = -log10(0.05), col= "grey")+
  geom_point(size= 2)+
  theme_minimal()+
  geom_col(width = 0.05, fill = "black")+
  coord_flip()+
  labs(x= "")+
  theme(axis.text.x = element_text(angle=00, hjust = 1, size= 10, color = "black"),
        axis.text.y = element_text(angle=00, hjust = 1, size= 10, color = "black")
        )
cell_prop_ps

p.mean= plot_mean_proportions(real_prop$groupwise, "celltype_mod", "Celltype Composition")
p.mean

p.nulls= get_null_visuals(real_prop2, nulls)

pdf("output/figures/proportions_test_nulls.pdf",
    width= 10)
  p.nulls
dev.off()

pdf("output/figures/proportions_test_ps.pdf",
    height= 2.5,
    width = 5)
cell_prop_ps
dev.off()




# CELL STATES ---------------------------------------------------------------------------------

### fibroblasts
wrap_cell_state= function(meta_seu, cell,  celltype_col= "celltype" , t= 1000){

  meta_seu2= meta_seu %>%
    filter(!!as.symbol(celltype_col)== cell)

  nulls= get_null(t= t,
                  x=  meta_seu2[c("orig.ident", "cellstate")],
                  clust.col = "cellstate")
  print(nulls)
  ## get real observed proportions
  real_prop= calc_mean_proportions_per_group(meta_seu2,
                                             "cellstate",
                                             "group")

  real_prop2 =real_prop$groupwise%>%
    distinct(cellstate, group, mean.percent)%>%
    pivot_wider(names_from = group, values_from=  mean.percent)%>%
    mutate(delta.p = hfpef-ct)
  print(real_prop2)

  ps = get_empirical_p(real_prop = real_prop2, nulls =  nulls, clust.col = "cellstate")
  print(ps)

  prop_ps= enframe(ps, name = "cellstate", value = "p.val")%>%
    ggplot(., aes(x= cellstate, y= -log10(p.val)))+
    geom_hline(yintercept = -log10(0.05), col= "grey")+
    geom_point()+
    theme_minimal()+
    geom_col(width = 0.05, fill = "black")+
    coord_flip()+
    labs(x= "")+
    theme(axis.text.x = element_text(angle=00, hjust = 1, size= 8),
          axis.text.y = element_text(angle=00, hjust = 1, size= 8))

  return(list(prop_ps, ps, real_prop2))

}

meta_seu = seu[[]]
p.macs= wrap_cell_state(meta_seu= meta_seu,
                        cell =  "macrophages",
                        t= 1000)
p.macs
meta_seu = seu[[]]

p.fibs= wrap_cell_state(meta_seu = meta_seu,
                        cell =  "Fibroblasts",
                        t= 1000)

pdf("output/figures/cell_type_analysis/cellstate_macs_fibs_p_vals.pdf",
    height= 1.3,
    width= 3.5)
p.macs
p.fibs
dev.off()


# integrated CELL STATES (by study) -----------------------------------------------------------

int.fibs= readRDS("output/seu.objs/study_integrations/harmony_fib_filt.rds")
meta_seu= int.fibs@meta.data
rm(int.fibs)


wrap_cell_state_int= function(meta_seu,
                              filter_val= "forte",
                              filer_col= "study" ,
                              group_col= "group",
                              t= 1000,
                              clust_col ="opt_clust_integrated",
                              ...){

  meta_seu2= meta_seu %>%
    filter(!!as.symbol(filer_col)== filter_val)

  nulls= get_null(t= t,
                  x=  meta_seu2[c(group_col, clust.col)],
                  clust_col = clust_col,group_col = group_col)
  test_normality(nulls)

  ## get real observed proportions
  real_prop= calc_mean_proportions_per_group(meta_seu2,
                                             clust.col,
                                             group_col)

  vals.contrast= real_prop$groupwise$group%>% unique()
  real_prop2 =real_prop$groupwise%>%
    distinct(!!as.symbol(clust_col), !!as.symbol(group_col), mean.percent)%>%
    pivot_wider(names_from = group, values_from=  mean.percent)%>%
    mutate(delta.p = hf-ct)
  print(real_prop2)

  ps = get_empirical_p(real_prop = real_prop2, nulls =  nulls, clust_col = clust_col)
  print(ps)

  prop_ps= enframe(ps, name = "opt_clust_integrated", value = "p.val")%>%
    ggplot(., aes(x= opt_clust_integrated, y= -log10(p.val)))+
    geom_hline(yintercept = -log10(0.05), col= "grey")+
    geom_point()+
    theme_minimal()+
    geom_col(width = 0.05, fill = "black")+
    coord_flip()+
    labs(x= "")+
    theme(axis.text.x = element_text(angle=00, hjust = 1, size= 8),
          axis.text.y = element_text(angle=00, hjust = 1, size= 8))
  return(( ps))

}


res= map(unique(meta_seu$study), function(x){
  wrap_cell_state_int(meta_seu,
                      filter_val= x,
                      filer_col= "study" ,
                      group_col= "group",
                      t= 1000,
                      clust_col ="opt_clust_integrated")

})


names(res)= unique(meta_seu$study)
p.val.df= map(res, function(x){
  names(x)= sort(unique(meta_seu$opt_clust_integrated))
  enframe(x, name= "int_cellstate", value= "p.val")
  })


saveRDS(p.val.df, "output/fib_integration/p.vals.proportional.rds")
p.val.df=readRDS( "output/fib_integration/p.vals.proportional.rds")
