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
## compute funcomics per pseudobulked cell type.#
# to do: the sign thing is probably not the best way, find another way( log2FC  )#
# SUBSAMPLE BEFORE T-testing on single cells otherwise samples with less cells are under represented!!


library(Seurat)
library(tidyverse)
library(dorothea)
library(progeny)
library(devtools)
library(decoupleR)
library(ComplexHeatmap)
library(rstatix)
library(cowplot)
library(biomaRt)
source("analysis/utils.R")

pb = readRDS("output/cell_specific_obj/cell_type_based_DE/samplewise_cellwise_PB.rds")



# Prepare reg.networks -------------------------------------

#### dorothea reg:

regulons= dorothea_mm %>%
  filter(confidence %in% c("A", "B", "C")) %>%
  mutate(likelihood=1)

### fibs.reg from heart cell atlas human

fib.reg= read.csv("data/regulons_heart_atlas/Fib.txt", sep = '\t')

endo.reg= read.csv("data/regulons_heart_atlas/Endo.txt", sep = '\t')

# translate regulons to mouse with biomart:

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

ends.genes= convertHumanToMouse(endo.reg$target, human, mouse)

translate_regulons= function(regulon){

  genes= convertHumanToMouse(regulon$target, human, mouse)
  regulon.m= regulon %>% left_join(genes$full.genes %>% rename(target= HGNC.symbol)) %>% drop_na %>% as_tibble

  regulon.m= regulon.m %>%
    dplyr::select(-target, -n_reads, -max_reads) %>%
    rename(target= MGI.symbol,
           tf= source)

  regulon.m= regulon.m %>%
    distinct(tf, target, likelihood, mor) %>%
    group_by(tf, target) %>%
    summarise(likelihood= mean(likelihood),
              mor = mean(mor))

  return(regulon.m)
}

end.reg.m= translate_regulons(endo.reg)
fib.reg.m= translate_regulons(fib.reg)

# end.reg.m= end.reg.m %>%
#   dplyr::select(-target, -n_reads, -max_reads) %>%
#   rename(target= MGI.symbol,
#          tf = source)

regulons_heart_m= list("fib"= fib.reg.m,
                       "end"= end.reg.m)

saveRDS(regulons_heart_m, "data/regulons_heart_atlas/regulons_heart_m.rds")
regulons_heart_m= readRDS("data/regulons_heart_atlas/regulons_heart_m.rds")



# Run TF activity on pseudobulk: --------------------------------------------------------------
run_TMM_= function(pb, celltype){
  library(edgeR)
  library(limma)

  #subset pb object to desired celltype
  df= sapply(pb, function(x){
    x[, celltype]
  })


  #### Filter & Normalize
  #create DGE class object
  group <- c("ct","ct","hf","hf")
  print("make sure this aligns:")
  print(group)
  print(colnames(df))
  dge <- DGEList(counts=df, group=group)

  #detect and remove low expressed gene
  keep <- filterByExpr(dge,group = group )

  dge <- dge[keep,,keep.lib.sizes=FALSE]

  dge <- calcNormFactors(dge)

  # use limma voom to transform count data to log2 counts per million
  v <- voom(dge, plot=TRUE)

  ### save R-objects
  voom_count= v$E

  return(voom_count)
}

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

run_progeny= function(gex.profile, .label){
  prog.res= progeny(expr = gex.profile, organism = "Mouse", scale =T, perm= 1000,z_scores =T)
  p.progeny.z= Heatmap(t(prog.res), cluster_columns = F, name = .label)
}

## normalize pseudobulks:
fibs= run_TMM_(pb, "Fibroblasts")
endos= run_TMM_(pb, "Endothelial")
macs= run_TMM_(pb, "Macrophages")
## run TF activity on pseudobulks

#with dorothea GRNs
fib.pb.dorothea= run_dorothea(gex.profile = fibs, regulons = regulons)
end.pb.dorothea= run_dorothea(endos, regulons)

## cell type specific GRN
fib.pb.grn= run_dorothea(fibs,regulons_heart_m$fib)
end.pb.grn= run_dorothea(endos,regulons_heart_m$end)


#plot lfc without p-val

plot_pb_func= function(tf, .label){
  tf %>%
    filter(abs(logfc) >0.5) %>%
    ggplot(.,aes(x= reorder(tf,logfc), y= logfc))+
    geom_col()+
    scale_fill_gradient(low = "red", high= "black")+
    #scale_fill_manual(values = col.set)+
    theme_classic()+
    coord_flip()+
    theme(axis.text = element_text( hjust = 1, size = 12))+
    labs(x= "tf", y= "log2fc")+
    ggtitle(paste0("pseudbulked", .label))

}
p1 = plot_pb_func(fib.pb.dorothea$tf_log2fc, "_fibroblast_dorothea")
p2 = plot_pb_func(end.pb.dorothea$tf_log2fc, "_end_dorothea")
p3 = plot_pb_func(fib.pb.grn$tf_log2fc, "_fibroblast_grn")
p4 = plot_pb_func(end.pb.grn$tf_log2fc, "_end_grn")

plot_grid(p1,p2,p3,p4, ncol = 2)

#plot_progny
p5= run_progeny(fibs, "fibroblasts")
p6= run_progeny(endos, "endothelial")

pdf("output/figures/funcomics/pseudobulked_progeny_per_celltype.pdf")
p5
p6
dev.off()

DEgenes_cells= readRDS("output/cell_specific_obj/cell_type_based_DE/pseudobulk_DE.rds")


df= DEgenes_cells$fibroblasts

df= DEgenes_cells$Endothelial
df= as.data.frame(df$dea) %>% dplyr::select(t, gene) %>% column_to_rownames("gene")
x= run_progeny(as.matrix(df), "DE_genes")

prog.res= progeny(as.matrix(df), organism = "Mouse", z_scores = T, perm = 500)
prog.res

# run dorothea on single cell level for those that are interesting from the pseudobulk --------

## Now test for those TFs of interest:
seu_int= readRDS("output/seu.objs/integrated_funcomics.rds")

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
  fc.tf.activity= dec.res2 %>% group_by(source, group) %>% summarize(mean_= mean(score)) %>% pivot_wider(names_from= group, values_from= mean_)%>%
    mutate(fc =hf/ct,
           logfc= log2(fc))

  # perform t-test per tf per group
  dec.res2 = dec.res2 %>% group_by(source) %>% wilcox_test(score~group, ref.group	="hf") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  dec.res2= dec.res2 %>% left_join(fc.tf.activity, by= "source")

  return(dec.res2)
}

plot_enrichment_with_logfc = function(tf.res,
                                      pval.cutoff= 0.05, label.= "celltype"){

  tf.res %>% filter(p.adj<pval.cutoff) %>%
    ggplot(aes(x= reorder(source,logfc), y= logfc, fill= -log10(p.adj)))+
    geom_col()+
    scale_fill_gradient(low = "red", high= "black")+
    #scale_fill_manual(values = col.set)+
    theme_classic()+
    coord_flip()+
    theme(axis.text = element_text( hjust = 1, size = 12))+
    labs(x= "tf", y= "log2fc")+
    ggtitle(label.)

}



## DOROTHEA GRNS

#using dorothea GRN with only those tfs that showed consistent regulation on pseudobulk level:
tf_oi= fib.pb.dorothea$tf_log2fc %>% filter(abs(logfc)>0.5) %>% pull(tf)
tfs= unique(regulons$tf)

fib.sc.dorothea= test_on_cell_level("fibroblasts",
                              tf_oi = tfs,
                              seu_int,
                              regulons)

p.fib.dorothea= plot_enrichment_with_logfc(fib.sc.dorothea,label. =  "fibroblasts_dorothea")

# estimate whether pseudobulk fc agrees with sc level
sig.TFs= fib.sc.dorothea%>% filter(p.adj<0.1) %>% pull(source)
fib.pb.dorothea$tf_log2fc %>% filter(tf %in% sig.TFs)

# endos:
tf_oi= end.pb.dorothea$tf_log2fc %>% filter(abs(logfc)>0.5) %>% pull(tf)
end.sc.dorothea = test_on_cell_level(cell_oi = "Endothelial",
                               #tf_oi = endos.res$tf,
                               tf_oi= tf_oi,
                               seu_int,
                               regulons)

p.end.dorothea= plot_enrichment_with_logfc(end.sc.dorothea,label. =  "endothelial_dorothea")
#using dorotheaGRN with  ALL tfs [lots of FP] probably cannot be run locally:
# fibs.res2= test_on_cell_level("fibroblasts",
#                               regulons$tf,
#                               seu_int,
#                               regulon = regulons)
# end.res2= test_on_cell_level("Endothelial",
#                              regulons$tf,
#                              seu_int,
#                              regulon = regulons)

# CELL SPEC REGULONS:
tf_oi= fib.pb.grn$tf_log2fc %>% filter(abs(logfc)>0.3) %>% pull(tf)

fib.sc.grn= test_on_cell_level("fibroblasts",
                              tf_oi = tf_oi,
                              seu_int,
                              regulon = regulons_heart_m$fib)

p.fib.grn= plot_enrichment_with_logfc(fib.sc.grn,label. =  "fibroblast_grn")

tf_oi= end.pb.grn$tf_log2fc %>% filter(abs(logfc)>0.3) %>% pull(tf)

end.sc.grn= test_on_cell_level("Endothelial",
                            tf_oi = tf_oi,
                            seu_int,
                            regulon = regulons_heart_m$end)

p.end.grn= plot_enrichment_with_logfc(end.sc.grn,label. =  "endothelial_grn")


pdf("output/figures/funcomics/TF_overview_log2fc.pdf",
    height= 4, width = 8
    )
plot_grid(p.fib.grn, p.fib.dorothea)
plot_grid(p.end.grn, p.end.dorothea)
dev.off()


# Plot results:  ------------------------------------------------------------------------------

# we will plot TFs significant on sc level as scaled normalized pb :
#DOROTHEA:

tfs.target= endos.res2%>% filter(p.adj <0.05) %>% pull(source)

endos.res$tf_log2fc %>%
  filter(tf %in% tfs.target) %>%
  ggplot(.,aes(x= reorder(tf, logfc), y= logfc))+
  geom_col()

tfs.target= fibs.res2%>%
  filter(p.adj <0.05) %>%
  pull(source)

fibs.res$tf_log2fc %>%
  #filter(tf %in% tfs.target) %>%
  ggplot(.,aes(x= reorder(tf, logfc), y= logfc))+
  geom_col()

df.fibs.dorothea= fibs.res$df_scaled[fibs.res2%>% filter(p.adj <0.05) %>% pull(source), ]

Heatmap(df.dorothea)

vec= fibs.res2%>% filter(p.adj <0.05) %>% arrange(p.adj) %>% pull(source)
pdf("output/figures/funcomics/TF_dorothea_fibs_overview.pdf", height= 20, width=5)
fibs.res$P[3]
Heatmap(df.fibs.dorothea[vec,], cluster_rows = F)
dev.off()


df.end.dorothea= endos.res$df_scaled[endos.res2%>% filter(p.adj <0.05) %>% pull(source), ]
Heatmap(df.end.dorothea)

pdf("output/figures/funcomics/TF_dorothea_end_overview.pdf", height= 20, width=5)
fibs.res$P[3]
Heatmap(df.end.dorothea)
dev.off()

# CELL SPECIFIC GRNS:

df[]

df.end= end.end$df_scaled[
  tf.endo %>% filter(p.adj <0.05) %>% pull(source),]
Heatmap(df2)
end.end$P[3]

pdf("output/figures/funcomics/TF_endoGRN_overview.pdf", height= 20, width=5)
end.end$P[3]
Heatmap(df.end)
dev.off()


df.fibs= fibs.fibs$df_scaled[
  tf.fibs%>% filter(p.adj <0.05) %>% pull(source),]

pdf("output/figures/funcomics/TF_fibGRN_overview.pdf", height= 20, width=5)
  fibs.fibs$P[3]
  Heatmap(df.fibs)

dev.off()
#
# p.fibs = tf.fibs %>% filter(p.adj.signif != "ns") %>%
#   arrange(statistic) %>%
#   ggplot(., aes(x= .y., y= reorder(source,(statistic)),  fill = statistic))+
#   geom_tile()+
#   scale_fill_gradient2(low= "blue", mid= "white", high = "red")+
#   theme_classic()+
#   geom_label(mapping=aes(label= p.adj.signif),label.size = NA)+
#   labs(x= "fibroblasts_ HFpEF contrast")
#
#
# p.endo= tf.endo %>%
#   arrange(statistic) %>%
#   ggplot(., aes(x= .y., y= reorder(source,(statistic)),  fill = statistic))+
#   geom_tile()+
#   scale_fill_gradient2(low= "blue", mid= "white", high = "red")+
#   theme_classic()+
#   geom_label(mapping=aes(label= p.adj.signif),label.size = NA)+
#   labs(x= "endothelial cells_ HFpEF contrast")
#
#   endos.res$P[2]
#   endos.res$P[3]
#   p.endo
#
# pdf("output/figures/funcomics/TF_dorothea_overview_fibs_.pdf",
#     height=12, width= 6)
#   #fibs.res$P[2]
#   fibs.res$P[3]
#   p.fibs
#   #cowplot::plot_grid(p.fibs, p.endo)
# dev.off()
#
# pdf("output/figures/funcomics/TF_dorothea_overview_endo_.pdf",
#     height=12, width= 6)
# #fibs.res$P[2]
# endos.res$P[3]
# p.endo
# #cowplot::plot_grid(p.fibs, p.endo)
# dev.off()
#
