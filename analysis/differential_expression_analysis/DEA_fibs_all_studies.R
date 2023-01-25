## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-04
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## run a DEA for fibroblasts from all studies


library(Seurat)
library(tidyverse)
library(ggvenn)
source("code/utils.R")
# circ=
circ= readRDS( file = "../sc-exploration/output/circ_obj/fibroblasts_seu.rds")
DefaultAssay(circ)= "RNA"
circ@meta.data= circ@meta.data %>% mutate(group = ifelse(group == "angII", "hf", "ct" ))
circ@meta.data%>% select(treatment, group) %>% table()

#forte =
forte= readRDS( file = "../sc-exploration/output/cell_mi_files/fibro_subset_integrated.rds")
forte@meta.data= forte@meta.data %>% mutate(group = ifelse(OP == "myocardial infarction", "hf", "ct" ))
DefaultAssay(forte)= "RNA"

#hfpef=
hfpef= readRDS(file = "output/cell_specific_obj/fibroblasts.rds")
#hfpef@meta.data= hfpef@meta.data %>% mutate(group = ifelse(group == "hfpef", "hf", "ct" ))
DefaultAssay(hfpef)= "RNA"


# run DE --------------------------------------------------------------------------------------

DimPlot(circ)
DimPlot(forte)
DimPlot(hfpef)
Idents(circ)= "group"

circ.de= FindMarkers(circ, group.by = "group",ident.1 = "hf", ident.2 = "ct")
forte.de= FindMarkers(forte, group.by = "group",ident.1 = "hf", ident.2 = "ct")
hfpef.de= FindMarkers(hfpef, group.by = "group",ident.1 = "hf", ident.2 = "ct")



objs= list("hfpef"= hfpef.de,
           "forte"= forte.de,
           "circ"= circ.de)

saveRDS(objs, file = "output/fib_integration/marker_list/DEG_per_study_in_fibs.rds")
objs= readRDS( file = "output/fib_integration/marker_list/DEG_per_study_in_fibs.rds")

# downsample comparison2 ----------------------------------------------------------------------

get_subsampled_DEG= function(seu, seed, genes.to.test=NULL, prefilter= T){

  set.seed(seed)

  orig.ids= unique( seu@meta.data$orig.ident)

  cell.count= seu@meta.data %>%rownames_to_column("cellid") %>%
    group_by(orig.ident) %>%
    dplyr::count()

  random.cell.ids= map(orig.ids , function(x){
    df= seu@meta.data%>% filter(orig.ident== x)
    sample(rownames(df), min(cell.count$n), replace = F)
  })


  #print(cell.count)
  map((random.cell.ids), length)
  print(paste0("downsample to ", min(cell.count$n), " cells"))

  seu.sub= subset(seu, cells = unlist(random.cell.ids))
  DefaultAssay(seu.sub)= "RNA"
  Idents(seu.sub) = "group"
  if(prefilter==T){
    DEA= FindMarkers(seu.sub,
                     ident.1 = "hf",
                     ident.2 = "ct",
                     # pseudocount.use = 5,
                     #logfc.threshold = 0.1,
                     #features =intersect(rownames(seu.sub), genes.to.test),
                     #     min.pct = 0,
                     # min.cells.feature = 0,
                     # min.cells.group = 0,
                     # min.diff.pct = 0,
                     #min.pct	= 0,
                     pseudocount.use = 5,
                     logfc.threshold = 0.1)

  }else{
    DEA= FindMarkers(seu.sub,
                     ident.1 = "hf",
                     ident.2 = "ct",
                      pseudocount.use = 5,
                     logfc.threshold = 0,
                     #features =intersect(rownames(seu.sub), genes.to.test),
                          min.pct = 0,
                      min.cells.feature = 0,
                      #min.cells.group = 0,
                      #min.diff.pct = 0,
                     min.pct	= 0,
                     pseudocount.use = 5
                     )
  }
  print(dim(DEA))
  DEA = DEA %>%rownames_to_column("gene") %>% arrange(gene) %>% as_tibble()
}


wrap_deg_subsampled_with_fisher= function(seu,
                                          number_of_replicates= 10, ...){

  # for each repetition perform a downsampling with unique seed
  require(survcomp)

  deg.array=  sapply(seq(1:number_of_replicates), function(x){
    get_subsampled_DEG(seu, seed = x, genes.to.test = genes.to.test)
  }, simplify = F)

  #get all genes that were tested:
  tested.genes= map(deg.array, function(x){
    x$gene
  })

  tested.genes= sort(unique(unlist(tested.genes)))

 sig.genes= map(deg.array, function(x){
    x%>% filter(p_val_adj<0.05)%>% pull(gene)
  })

  count.genes= enframe(table(unlist(sig.genes)), name = "gene", value = "count")

  res = map(tested.genes, function(x){
    #print(x)
    gene_p_vec= lapply(seq(1:number_of_replicates), function(y){
      #print(deg.array[[y]])
      p.= deg.array[[y]] %>% filter(gene==x) %>% pull(p_val_adj)
      #print(p.)
    })
    #message("p-value ready")
    gene_fc_vec= lapply(seq(1:number_of_replicates), function(y){
      #print(deg.array[[y]])
      p.= deg.array[[y]] %>% filter(gene==x) %>% pull(avg_log2FC)
      #print(p.)
    })
    fish.p = combine.test(p= unlist(gene_p_vec), na.rm = T)
    median.fc= median(unlist(gene_fc_vec),na.rm = T)
    return(list(fish.p, median.fc))
  })

  names(res) = tested.genes
  head(res)
  df= enframe(res)# %>% mutate(value= unlist(value)) %>% rename(fisher.p= value) %>% arrange(fisher.p)
  df= df %>% unnest_wider(value )
  colnames(df)= c("gene", "fisher.p", "median.log2FC")
  df= left_join(df, count.genes)

  #optional filter for count

  df = df %>% filter(count>number_of_replicates*0.6)
  return(df%>% arrange(fisher.p))

}



# downsample for each study:: -----------------------------------------------------------------

## run for fibs
hfpef_f= wrap_deg_subsampled_with_fisher(seu = hfpef,number_of_replicates =  5)
forte_f= wrap_deg_subsampled_with_fisher(forte, 5)
circ_f= wrap_deg_subsampled_with_fisher(circ, 5)

#for the MI study we are adding the subset of early and late remodeling

forte_early= wrap_deg_subsampled_with_fisher(subset(forte,time != 14), 5)
forte_late= wrap_deg_subsampled_with_fisher(subset(forte,time %in% c(0,14) | group == "ct" ), 5)


#combine and make up and downregulated lists:
obj = list(hfpef_f, forte_f, circ_f)
#obj= list(forte_early,forte_late)
l.sets= map(obj, function(x) (x%>%filter(median.log2FC>0)%>% arrange(fisher.p, desc(abs(median.log2FC))))%>% #slice(1:100)
              pull(gene))
l.sets.dn= map(obj, function(x) (x%>%filter(median.log2FC<0)%>% arrange(fisher.p, desc(abs(median.log2FC))))%>% #slice(1:100)
              pull(gene))
names(l.sets)=names(l.sets.dn)= c("HFpEF", "MI", "AngII")

#plot venn of intersects:
source("analysis/utils.R") #load for col_vector
venn.up = ggvenn(l.sets, fill_color = col_vector, show_percentage = F,fill_alpha = 0.6, text_size = 5)
venn.dn = ggvenn(l.sets.dn, fill_color = col_vector, show_percentage = F,fill_alpha = 0.6, text_size = 5)
pdf("output/figures/integration_studies/deg.venns.pdf",
    width= 4,
    height=3)
venn.up
venn.dn
dev.off()


#find intersects and unions:

HFpEF_unique_sig= l.sets$HFpEF[!l.sets$HFpEF %in% c(l.sets$MI, l.sets$AngII)]
MI_unique_sig= l.sets$MI[!l.sets$MI %in% c(l.sets$HFpEF, l.sets$AngII)]
AngII_unique_sig=  l.sets$AngII[!l.sets$AngII %in% c(l.sets$HFpEF, l.sets$MI)]

overlap_all_sig= intersect(intersect(l.sets[[1]], l.sets[[2]]), l.sets[[3]])

HFpEF_MI_sig= intersect(l.sets$HFpEF, l.sets$MI)[!intersect(l.sets$HFpEF, l.sets$MI) %in% l.sets$AngII]
HFpEF_AngII_sig= intersect(l.sets$HFpEF, l.sets$AngII)[!intersect(l.sets$HFpEF, l.sets$AngII) %in% l.sets$MI]
MI_AngII_sig= intersect(l.sets$MI, l.sets$AngII)[!intersect(l.sets$MI, l.sets$AngII) %in% l.sets$HFpEF]

gene_signatures= list("unique"= list("HFpEF"= HFpEF_unique_sig,
                                     "MI"= MI_unique_sig,
                                     "AngII"= AngII_unique_sig),
                      "overlap"= list("all"= overlap_all_sig,
                                      "HFpEF_MI"= HFpEF_MI_sig,
                                      "HFpEF_AngII"= HFpEF_AngII_sig,
                                      "MI_AngII"= MI_AngII_sig),
                      "total"= l.sets
)

saveRDS(gene_signatures, "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
gene_signatures= readRDS( "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")

saveRDS(hfpef_f, "output/fib_integration/marker_list/DEG_hfpef.rds")


##expand with early and late for MI
obj= list(forte_early,forte_late)
l.sets= map(obj, function(x) (x%>%filter(median.log2FC>0)%>% arrange(fisher.p, desc(abs(median.log2FC))))%>% #slice(1:100)
              pull(gene))
gene_signatures$total$MI_late= l.sets[[2]]
gene_signatures$total$MI_early= l.sets[[1]]

saveRDS(gene_signatures, "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")

source("code/utils.R") #load for col_vector
venn.up = ggvenn(list("MI_late"= gene_signatures$total$MI_late,
                      "MI_early"= gene_signatures$total$MI_early,
                   "HFpEF"= gene_signatures$total$HFpEF,
                   "AngII" = gene_signatures$total$AngII), fill_color = c(col_vector[2], col_vector[1], col_vector[3],col_vector[4]), show_percentage = F,fill_alpha = 0.6, text_size = 4, set_name_size = 5)
venn.up

pdf("output/figures/main/Fig3/deg.venns.pdf",
    width= 4,
    height=3)
venn.up
dev.off()



# add jaccard ---------------------------------------------------------------------------------
gene_signatures$total= gene_signatures$total[names(gene_signatures$total)!="MI" ]
df= sapply(gene_signatures$total, function(x){
  sapply(gene_signatures$total, function(y){
        length(intersect(x,y))/length(union(x,y))
  })
})

p1= df %>% as.data.frame()%>% rownames_to_column("comparison")%>%
  pivot_longer(-comparison)%>%
  mutate(value = ifelse(name== comparison, NA, value),
       comparison = factor(comparison, levels= c("HFpEF", "AngII", "MI_early", "MI_late")),
       name= factor(name, levels= c("HFpEF", "AngII", "MI_early", "MI_late"))) %>%
  ggplot(., aes(x= name, y = comparison, fill = value))+
  geom_tile(color= "black")+
  scale_fill_gradient(low= "white", high= "darkred")+
  theme_minimal()+
  #geom_text(aes(label= sig))+
  theme(axis.text= element_text(colour = "black"),
        axis.text.x = element_text(angle= 45, hjust= 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  coord_equal()+
  labs(fill ="Jaccard \n Index",x="",  y= "")

pdf("output/figures/main/Fig3/jaccard_sigs2.pdf",
    height= 3,
    width= 3)
unify_axis(p1)
dev.off()


diag(df)= NA
df[upper.tri(df)]= NA

hmap= ComplexHeatmap::Heatmap(t(df),
                        cluster_rows = F, cluster_columns = F, name = "Jaccard Index",
                        rect_gp = gpar(col = "black", lwd = 1.5))

pdf("output/figures/main/Fig3/jaccard_sigs.pdf",
    height= 3,
    width= 3)
hmap
dev.off()

# Compare cell signatures and plot ------------------------------------------------------------

##1 = SET COMPARISON
# Take all DE genes upregulated and define unique genes for each study

#
# get_genes= map(gene_signatures$total, function(x) (rownames_to_column(x, "gene")  %>%filter(p_val_adj<0.01,
#                                               avg_log2FC>0) %>% pull(gene) ))
#
# pdf("output/figures/integration_studies/fibroblast.study.DEG.venn.pdf")
# ggvenn(gene_signatures$total)
# dev.off()
#
# hfpef_unique_sig= get_genes$hfpef[!get_genes$hfpef %in% c(get_genes$forte, get_genes$circ)]
# forte_unique_sig= get_genes$forte[!get_genes$forte %in% c(get_genes$hfpef, get_genes$circ)]
# circ_unique_sig=  get_genes$circ[!get_genes$circ %in% c(get_genes$hfpef, get_genes$forte)]
#
# overlap_all_sig= intersect(intersect(get_genes[[1]], get_genes[[2]]), get_genes[[3]])
#
# hfpef_forte_sig= intersect(get_genes$hfpef, get_genes$forte)[!intersect(get_genes$hfpef, get_genes$forte) %in% get_genes$circ]
# hfpef_circ_sig= intersect(get_genes$hfpef, get_genes$circ)[!intersect(get_genes$hfpef, get_genes$circ) %in% get_genes$forte]
# forte_circ_sig= intersect(get_genes$forte, get_genes$circ)[!intersect(get_genes$forte, get_genes$circ) %in% get_genes$hfpef]
#
# gene_signatures= list("unique"= list("hfpef"= hfpef_unique_sig,
#                                      "forte"= forte_unique_sig,
#                                      "circ"= circ_unique_sig),
#                       "overlap"= list("all"= overlap_all_sig,
#                                       "hfpef_forte"= hfpef_forte_sig,
#                                       "hfpef_circ"= hfpef_circ_sig,
#                                       "forte_circ"= forte_circ_sig)
#                       )
#
# saveRDS(gene_signatures, "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET.rds")
# gene_signatures= readRDS( "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET.rds")



# get full DEa results for logfc for all genes ------------------------------------------------


objs= list("hfpef"= hfpef,
           "forte"= forte,
           "circ"= circ)

df= map(objs, function(x){
  Idents(x)= "group"
  x_f= FoldChange(x, ident.1= "hf", ident.2= "ct")
})



saveRDS(df, file = "output/fib_integration/marker_list/DEG_per_study_LOGFC.rds")
df=readRDS( file = "output/fib_integration/marker_list/DEG_per_study_LOGFC.rds")
Idents(forte)= "group"
MI_early= FoldChange(subset(forte,time != 14), ident.1= "hf", ident.2= "ct")
MI_late= FoldChange(subset(forte,time %in% c(0,14)), ident.1= "hf", ident.2= "ct")
df$forte_late= MI_late
df$forte_early= MI_early
names(df)=c("HFpEF", "MI", "AngII", "MI_late", "MI_early")
df= df[names(df)!= "MI"]
saveRDS(df, file = "output/fib_integration/marker_list/DEG_per_study_LOGFC.rds")
#
# mean_exp= map(objs, function(x){
#   Idents(x)= "group"
#   x_f= AverageExpression(x, features = c(the_good, the_bad) ,assays = "RNA",  group.by= "group")
# })
#
# map(mean_expm , function(x){
#   x$RNA
# }
# mean_exp$hfpef$RNA



