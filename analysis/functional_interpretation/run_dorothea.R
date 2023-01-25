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
##


library(Seurat)
library(tidyverse)
library(dorothea)

library(devtools)
library(decoupleR)
library(ComplexHeatmap)
library(rstatix)
library(cowplot)
library(biomaRt)
source("analysis/utils.R")

seu= readRDS("output/seu.objs/cell.state.obj.list.rds")
logFCs = readRDS("output/fib_integration/marker_list/DEG_per_study_LOGFC.rds")

# Prepare reg.networks -------------------------------------

#### dorothea reg:
data("dorothea_mm")

regulons= dorothea_mm %>%
  filter(confidence %in% c("A", "B", "C")) %>%
  mutate(likelihood=1)%>%
  rename(source= tf)

tfs.corr= decoupleR::check_corr(network= regulons)%>% filter(correlation ==1)%>% pull(source)
regulons= regulons%>% filter(!source %in% tfs.corr)


x= merge(logFCs$HFpEF,logFCs$AngII,
         by="row.names",
         all.x=F,
         all.y= F,
         suffixes=c("HFpEF", "AngII"))

y= merge(logFCs$MI_late,logFCs$MI_early,
         by="row.names",
         all.x=F,
         all.y= F,
         suffixes=c("MI_late", "MI_early"))

xy= merge(x, y,by="Row.names",
      all.x=F,
      all.y= F )


df= xy[,grepl(x = colnames(xy),pattern = "avg_log2")]
rownames(df)= xy$Row.names

#remove NA

dec.res= decouple(mat = as.matrix(df),
            network = regulons,
            statistics = "mlm")

dec.res=
  dec.res %>% mutate(condition= str_replace_all(condition, "avg_log2FC", ""),
                     p_value= p.adjust(p_value,method= "BH"),
                     stars= ifelse(p_value<0.001, "***",
                                   ifelse(p_value<0.01, "**",
                                          ifelse(p_value<0.05, "*", "")
                                          )
                                   )
                     )

##plot sig TFs per study:
pls= map(unique(dec.res$condition), function(x){
  tfs= dec.res%>% filter(condition== x,
                         statistic== "mlm",
                    p_value <0.05)%>%
    arrange((score))%>%
    pull(source)

  p.df= dec.res %>%
    filter(source %in% tfs)%>%
    mutate(source= factor(source, levels = tfs),
           condition=factor(condition, levels= c("HFpEF", "AngII", "MI_early", "MI_late")))
  range(p.df$score)
  p.df%>% arrange(desc(score))
  p= p.df %>%
    ggplot(., aes(x= condition, y= source, fill = score,label =stars))+
    geom_tile(color ="darkgrey")+
    scale_fill_gradient2(low= "blue", mid= "white", high= "red", midpoint = 0, limits= c(-10,10))+
    theme_minimal()+
    geom_text(aes(label = stars))+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.text.x = element_text(angle= 45, hjust= 1))+
    coord_equal()+
    labs(y= "", x="", fill= "TF activity")

  unify_axis(p)

})

pdf("output/figures/main/Fig3/TF_acts_hfpef.pdf",
    width= 4,
    height= 6)
pls[[1]]
dev.off()


tfs= dec.res%>% filter(condition== "HFpEF",
                       statistic== "mlm",
                       p_value <0.05)%>%
  arrange((score))%>%
  pull(source)

p.df= dec.res %>%
  filter(source %in% tfs)%>%
  mutate(source= factor(source, levels = tfs),
         condition=factor(condition, levels= c("HFpEF", "AngII", "MI_early", "MI_late")))
range(p.df$score)
p.df%>% arrange(desc(score))

p=p.df %>%
  filter(statistic== "mlm")%>%
  ggplot(., aes(x= condition, y= source, fill = score,label =stars))+
  geom_tile(color ="darkgrey")+
  scale_fill_gradient2(low= "blue", mid= "white", high= "red", midpoint = 0)+
  theme_minimal()+
  geom_text(aes(label = stars))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle= 45, hjust= 1))+
  coord_equal()+
  labs(y= "", x="", fill= "TF activity")

unify_axis(p)

pdf(file = "output/figures/main/Fig3/TF_acts_hfpef2.pdf",
    width= 5,
    height= 7)
unify_axis(p)
dev.off()


tfdec.res %>% group_by(condition, source) %>% filter(any(p_value<0.05))%>% pull(source)

df= map(names(logFCs), function(gex1){

  gex= as.matrix(logFCs[[gex1]][,1])
  rownames(gex)= rownames(logFCs[[gex1]])

  x= decouple(mat = gex,
              network = regulons,
              statistics = "mlm")

  x %>% filter(statistic =="mlm") %>%
    arrange(desc(score))

  })


df$HFpEF%>% filter(p_value <0.05)%>%pull(source)


run_dorothea= function(gex.profile, regulons){

  dec.res = decouple(gex.profile,
                     network = regulons,
                     statistics ="wmean",
                     .source ="tf",
                     .target = "target",
                    consensus_score = F
                    )


  p1= dec.res  %>% filter(statistic== "corr_wmean")%>% ggplot( aes(x= condition, y= source, fill = score))+
    geom_tile()+
    scale_fill_gradient2(low= "blue", mid= "white", high= "red")

  print(p1)

  dec.res.matrix= dec.res %>%
    filter(statistic== "corr_wmean") %>%
    dplyr::select(score, source, condition) %>%
    pivot_wider(names_from = source, values_from = score)

  df= column_to_rownames(dec.res.matrix,"condition")
  df= scale(df)
  df= t(df)
  p2= Heatmap(df)

####  this part uses the raw corrected means to calculate fold changes between both groups and selects TF based on an absolute log2 FC > 1
  ## however we comparing the scale of the different samples here.
tf_fc2= t(column_to_rownames(dec.res.matrix, "condition")) %>% as.data.frame %>% rownames_to_column("tf") %>% as_tibble %>%
    mutate(sum1= (hf1+hf2),
           sum2= (ct1+ct2),
           fc= sum1/sum2,
           logfc = log2(fc)) %>%
    arrange(desc(logfc))

  df_match= as_tibble((df)) %>%
    mutate(tf= rownames(df))   %>%
    filter(sign(hf1)==sign(hf2),
           sign(ct1)== sign(ct2))

  df_match= column_to_rownames(df_match, "tf")

  p3= Heatmap(df_match)

  return(list("P"= list(p1, p2,p3),
              "tf"= rownames(df_match),
              "tf_log2fc"= tf_fc2,
              "df_scaled"= df,
              "df"= dec.res
              )
         )

  }

# run dorothea on single cell level for those that are interesting from the pseudobulk --------

## Now test for those TFs of interest:

seu_int= seu_int$fibroblasts

DefaultAssay(seu_int)= "RNA"
regulon = regulons

test_on_cell_level= function(cell_oi= "fibroblasts", tf_oi, seu_int, regulon){


  ## remove cells not of interest:
  sub.seu= subset(seu_int, celltype ==cell_oi)

  ## Subsampling before calculating TF statistics:

  cell.count= sub.seu@meta.data %>%
    rownames_to_column("cellid") %>%
    group_by(orig.ident) %>%
    filter(celltype== cell_oi) %>%
    count()

  # loop through every sample and sample the minimum number of cells per sample
  # taken from the cell.count:

  cell_ids= map(unique(sub.seu@meta.data$orig.ident), function(x){
    # all cells from sample x
    sub.seu2= subset(sub.seu, orig.ident==x)
    # sample minimum
    sub.seu2 <- sub.seu2[, sample(colnames(sub.seu2), size =min(cell.count$n), replace=F)]
    # return cell IDs
    colnames(sub.seu2)
  }
  )

  #subset cells to all those samplesd:
  sub.seu= sub.seu[, unlist(cell_ids)]

  # get normalized gex data:
  gex_fibs= GetAssayData(sub.seu, slot = "data")

  # subset regulons to provided tfs of interest
  regulons_filt= regulon %>% filter(tf %in% tf_oi)

  # run decoupler on the resulting ingredients
  dec.res = decouple(gex_fibs[,],
                     network = regulons_filt,
                     statistics ="wmean",
                     .source ="tf",
                     .target = "target",
                     consensus_score = F
                     )

  # Get cell ID to group mappings:
  groups= sub.seu@meta.data
  translate= cbind("condition"= rownames(groups), "group"= groups$group) %>% as_tibble

  # map cell_ids into the tf activity tibble:
  dec.res2= dec.res %>%
    filter(statistic== "corr_wmean")  %>%
    left_join(translate, by= "condition")%>%
    mutate(group= factor(group, levels = c("hf", "ct")))

  # calculate logfc
  fc.tf.activity= dec.res2 %>%
    group_by(source, group) %>%
    summarize(mean_= mean(score)) %>%
    pivot_wider(names_from= group, values_from= mean_)%>%
    mutate(fc =hf/ct,
           logfc= log2(fc))

  # perform t-test per tf per group
  dec.res2 = dec.res2 %>% group_by(source) %>% wilcox_test(score~group, ref.group	="hf") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")

  dec.res2= dec.res2 %>% left_join(fc.tf.activity, by= "source")

  return(list("proc"= dec.res2,
              "raw"= dec.res))
}

plot_enrichment_with_logfc = function(tf.res,
                                      pval.cutoff= 0.05, label.= "celltype",
                                      logfc.cutoff =0.3){

  tf.res %>%
    filter(p.adj<pval.cutoff,
           abs(logfc)>logfc.cutoff) %>%
    ggplot(aes(x= reorder(source,logfc), y= logfc, fill= -log10(p.adj)))+
    geom_col()+
    scale_fill_gradient(low = "red", high= "black")+
    #scale_fill_manual(values = col.set)+
    theme_classic()+
    coord_flip()+
    theme(axis.text = element_text( hjust = 1, size = 12))+
    labs(x= "", y= "log2fc")+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1)
    )

}



## DOROTHEA GRNS

#using dorothea GRN with only those tfs that showed consistent regulation on pseudobulk level:
tf_oi= fib.pb.dorothea$tf_log2fc %>% filter(abs(logfc)>0.5) %>% pull(tf)
tfs= unique(regulons$tf)[1:10]

fib.sc.dorothea= test_on_cell_level("fibroblasts",
                              tf_oi = tfs,
                              seu_int,
                              regulons)

fib.sc.dorothea$raw
p.fib.dorothea= plot_enrichment_with_logfc(fib.sc.dorothea,
                                           label. =  "fibroblasts_dorothea",logfc.cutoff = 0.25
                                           )

seu_fun= readRDS("output/seu.objs/integrated_funcomics.rds")
seu_fun_fib= subset(seu_fun, celltype=="fibroblasts")
seu_fun@assays$dorothea
DefaultAssay(seu_fun_fib)= "dorothea"
FeaturePlot(seu_fun_fib, features=c("Ppara", "Gli2", "Rfx5"),split.by = "group")


pdf("output/figures/funcomics/TF_fibroblasts.pdf",
    height= 4, width = 4
    )
p.fib.dorothea
dev.off()

