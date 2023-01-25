## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-12-02
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## process murine hfpef data
## split the samples by batch

library(limma)
library(edgeR)
library(tidyverse)
library(ggrepel)
library(stringr)
library(sva)

target= read.csv("data/BulkRNAseq/Samples bulkSeq 11_2022.csv", skip = 1) %>%
  as_tibble()%>%
  select(Sample.ID, Group, X.4)%>%
  rename(factor1= X.4,
         sample= Sample.ID)%>%
  mutate(Group2= ifelse(grepl("HFpEF", Group), "HFpEF", "Control"))
target
counts= read.csv("data/BulkRNAseq/aligend_counts_time_course_hfpef.txt", sep = '\t')

colnames(counts)
df= column_to_rownames(counts, "X")
boxplot(df)
a = colnames(df)
res <- str_match(a, "_\\s*(.*?)\\s*_")[,2]
colnames(df)= res
colnames(df)== target$sample
df= df[,colnames(df)[match( target$sample, colnames(df))]]
colnames(df)== target$sample

target= target%>% mutate(factor1= ifelse(factor1 == "", 0, 1),
                         factor1= as.factor(factor1))
target$Group= map(target$Group, function(x){unlist(str_split(x, " ")[1])[1]}) %>% unlist()
target$Group= str_replace_all(target$Group, "-", "_")

target

## filter out
targets= split.data.frame(target, target$factor1)
df0= df[, targets$`0`$sample]
df1= df[, targets$`1`$sample]

# filter, TMM normalize and voom transform
get_dge= function(target, df){
  group= target$Group
  dge <- DGEList(counts=df, group = group)

  #detect and remove low expressed gene
  keep <- filterByExpr(dge, min.prop =0.5,min.count = 5 )
  #filterByExpr
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge$counts= dge$counts+5
  dge <- calcNormFactors(dge)
  # use limma voom to transform count data to log2 counts per million
  v <- voom(dge, plot=TRUE)

  ### save R-objects
  voom_count= v$E

  return(list("dge"= dge, "voom"= voom_count, "target"= target))
}



dge0 = get_dge(targets$`0`, df0)
dge1=  get_dge(targets$`1`, df1)

# PCA ------------------------------------------------------------------------------------

do_PCA= function(target, df){
  ##PCA
  PCA <- prcomp(t(df) ,center = TRUE, scale. = F)

  plot.pca = PCA$x %>%
    as.data.frame %>%
    rownames_to_column("sample") %>%
    as_tibble()%>%
    left_join(target)

  p.pca = ggplot(plot.pca,aes(x= PC1, y= PC2, col= Group, shape= factor1))+
    geom_point(size= 3)+
    theme_minimal()+
    labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
         y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
    ggtitle(paste0(""))+
    geom_text_repel(aes(label= sample),show.legend = FALSE)

  p.pca2 = ggplot(plot.pca,aes(x= PC3, y= PC4, col= Group, shape= factor1))+
    geom_point(size= 3)+
    theme_minimal()+
    labs(x= paste0("PC3 (",as.character(round(PCA$sdev[3]^2/sum(PCA$sdev^2)*100)),"%)"),
         y= paste("PC4 (",as.character(round(PCA$sdev[4]^2/sum(PCA$sdev^2)*100)),"%)"))+
    ggtitle(paste0(""))+geom_text_repel(aes(label= sample),show.legend = FALSE)
  #p.pca2

  cowplot::plot_grid(p.pca, p.pca2)

}

do_PCA(dge1$target, dge1$voom)
do_PCA(dge0$target, dge0$voom)


# DEA  ---------------------------------------------------------------------------------------------
## 5wks

dge1$target
boxplot(dge1$voom)

ct= dge1$voom[,c("H14")]
hf = dge1$voom[,colnames(dge1$voom) != "H14"]

# mean of hf group
hf_mean= apply(hf,1, mean)

# check if genes are orderd alike (sanity check) and calculate fold change (fc)
if(table(names(hf_mean)== names(ct))){
  fc= hf_mean/ct
}

# we cannot do more than logfc since strong batch in this lot and only one control

## 10, 15 wks do proper limma
# do contrast on disease (f1) or timepoints (f2)
f1 = factor(dge0$target$Group2, levels= c("Control", "HFpEF"))
f2 = factor(dge0$target$Group, levels= c("Control",
                                   "10_weeks",
                                   "15_weeks"
                                   ))


design = model.matrix(~f1)
cont.matrix = makeContrasts(hf= "f1HFpEF",
                              levels=design)
design.t = model.matrix(~f2)
cont.matrix.t = makeContrasts(k10= "f210_weeks",
                            wk15="f215_weeks",
                            levels=design.t)


get_fit= function(design,cont.matrix, voom_count){
  fit = lmFit(voom_count, design)
  #Empirical Bayes
  fit2 = contrasts.fit(fit, cont.matrix)
  fit2 = eBayes(fit2)
}

fit.hf = get_fit(design, cont.matrix, dge0$voom)
fit.t = get_fit(design.t, cont.matrix.t, dge0$voom)

## better to subset the df for the fit

w_10= dge0$target%>% filter(Group != "15_weeks")
f= factor(w_10$Group, levels= c("Control", "10_weeks"))
design10= model.matrix(~0+f)
cont.matrix.10= makeContrasts("wk10" = f10_weeks-fControl,
                              levels= design10)

w_15= dge0$target%>% filter(Group != "10_weeks")
f= factor(w_15$Group, levels= c("Control", "15_weeks"))
design15= model.matrix(~0+f)
cont.matrix.15= makeContrasts("wk15" = f15_weeks-fControl,
                              levels= design15)

fit10= get_fit(design10, cont.matrix.10, voom_count = dge0$voom[,w_10$sample])
fit15= get_fit(design15, cont.matrix.15, voom_count = dge0$voom[,w_15$sample])

saveRDS(list("5wk"= fc,
             "older"= list("hf"= fit.hf,
                          "time"= fit.t,
                          "w10"= fit10,
                          "w15"= fit15)),
        "output/bulk_mouse_fit.rds")




# ---------------------------------------------------------------------------------------------

get_toptable= function(g.type = "WT",
                       trtmnt = c("TAC", "Sham"),
                       contrast.name){

  meta= gex.obj$meta %>% filter(genotype %in% g.type,
                                treatment %in% trtmnt)

  #genotype= factor(meta$genotype, levels= c("WT", "S100a1", "Strit", "Strit_S100a1"))
  genotype= factor(meta$genotype, levels= g.type)

  #treatment= factor(meta$treatment, levels=c("Sham", "TAC"))
  treatment= factor(meta$treatment, levels= trtmnt)

  if(length(g.type)==1){
    design = model.matrix(~0+treatment)
  }else if(length(trtmnt)==1){
    design = model.matrix(~0+genotype)
  }
  print(design)
  #line= factor(meta$lane_id)

  fit = lmFit(gex.obj$exp[,meta$sample.names], design)
  #colnames(gex.obj$exp)

  # #Define contrasts
  # cont.matrix = makeContrasts(c.name= contrast.name,
  #                             levels=design)

  #makeContrast doesnt accept string, so weird work around here:
  cmd <- paste("cont.matrix <- makeContrasts(", contrast.name, ", levels = design)", sep =
                 '"')
  eval(parse(text = cmd))

  print(cont.matrix)
  #Empirical Bayes
  fit2 = contrasts.fit(fit, cont.matrix)
  fit2 = eBayes(fit2)

  #obtain differentially expressed genes
  DE_results2 = as.data.frame(topTable(fit2,adjust.method = "fdr",number = Inf)) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(desc(abs(t))) %>%
    as_tibble()


}



#obtain differentially expressed genes
DE_results = as.data.frame(topTable(fit2, coef = 2, adjust.method = "BH",number = Inf)) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(desc(abs(t))) %>%
  as_tibble()

DE_results%>%
  mutate(sig= ifelse(P.Value <0.05, "y", "n"),
         hfpef= ifelse(gene %in% gene_sigs$total$HFpEF, gene, ""),
         angii= ifelse(gene %in% gene_sigs$total$AngII, gene, ""))%>%
  ggplot(.,aes(x= logFC, y= -log10(P.Value), col =sig, label = angii))+
  geom_point()+
  geom_label_repel(aes(label = hfpef), max.overlaps = 1000, col = "black")


fit2$t
design2 = model.matrix(~0+f)

gene_sigs= readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
