## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-09-10
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## get pseudobulks of fibroblasts



library(Seurat)
library(tidyverse)
library(Matrix.utils)

seu= readRDS( file = "output/seu.objs/study_integrations/harmony_fib_filt.rds")

meta.integrated= seu@meta.data
#get pseudobulks of fibroblasts for each study

studies= unique(seu@meta.data$study)

# pseudobulk_fibroblasts ----------------------------------------------------------------------


x= studies[1]
pb.studies= lapply(studies, function(x){
  seu2= subset(seu, study==x)
  gex.count= GetAssayData(seu2, assay = "RNA", slot = "count")
  groups= seu2@meta.data[, c("orig.ident")]
  pb <- aggregate.Matrix(t(gex.count),
                         groupings = groups,
                         fun= "sum")



  })

rm(seu)

names(pb.studies) = studies

pb.studies2= map(pb.studies, t)
saveRDS(pb.studies2, "output/fib_integration/pseudobulked.per.study.sample.rds")

# pb= pb.studies2$circ
# x= meta.integrated %>% filter(study== "circ")
run_TMM_limma= function(pb, meta.integrated){
  library(edgeR)
  library(limma)
  library(ggrepel)

  #subset pb object to desired celltype
  df= pb
  target= meta.integrated %>% distinct(orig.ident, group)

  df= df[,target$orig.ident]
  colnames(df)== target$orig.ident
  group= target$group

  #### Filter & Normalize
  #create DGE class object

  # which(df==0)
  # rowSums(df)
  #
  #   df[which(df==0),]
  dge <- DGEList(counts=df, group=group)

  #detect and remove low expressed gene
  keep <- filterByExpr(dge,group = group, min.prop =1,min.count = 5 )
  #?filterByExpr
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge$counts= dge$counts + 5

  dge <- calcNormFactors(dge)
  # use limma voom to transform count data to log2 counts per million
  v <- voom(dge, plot=TRUE)

  ### save R-objects
  voom_count= v$E

  #Adjust a linear model to the expression data
  f = factor(group, levels= c("ct","hf"))
  design = model.matrix(~0+f)
  colnames(design) = c("ct","hf")
  ExpMat= voom_count
  fit = lmFit(ExpMat, design)

  #Define contrasts
  cont.matrix = makeContrasts(KO_sign = hf-ct,
                              levels=design)

  #Empirical Bayes
  fit2 = contrasts.fit(fit, cont.matrix)
  fit2 = eBayes(fit2)

  #obtain differentially expressed genes
  DE_results = as.data.frame(topTable(fit2,adjust.method = "BH",number = Inf)) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(desc(abs(t))) %>%
    as_tibble()

  p.pval.dist= DE_results %>% ggplot(., aes(x= P.Value))+geom_histogram(bins = 100)

  p.volcano= DE_results%>% ggplot(., aes(x= logFC, y= -log10(P.Value)))+
    geom_point()

  #PCA

  PCA <- prcomp(t(voom_count) ,center = TRUE, scale. = T)

  plot.pca = PCA$x %>%
    as.data.frame %>%
    rownames_to_column("orig.ident") %>%
    as_tibble()%>%
    left_join(target)

  p.pca = ggplot(plot.pca,aes(x= PC1, y= PC2,color = group))+
    geom_point(size= 3)+
    theme_minimal()+
    labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
         y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
    ggtitle(paste0("pseudobulked "))+
    geom_text_repel(aes(label= orig.ident),show.legend = FALSE)


  return(list("p" = p.volcano,
              "dea"= DE_results,
              "p2"= p.pca,
              "pca.df"= PCA,
              "pca.pca"= plot.pca
  )
  )

}


pb.de = lapply(studies, function(x){
  meta= meta.integrated %>% filter(study== x)
  pb.dea= run_TMM_limma(pb.studies2[[x]], meta)
})

names(pb.de)= studies

pb.de$circ
saveRDS(list("pseudobulkedGEX"= pb.studies2,
             "dea_res"= pb.de),
        "output/fib_integration/pseudobulked.per.study.sample.rds")




df= logFC.c%>% filter(grepl("Itg", gene))%>% select(-pct.1, -pct.2) %>%
  pivot_wider(names_from = study,values_from = avg_log2FC)%>% column_to_rownames("gene")

Heatmap(df,name = "avg_log2FC")


# pseudobulk_fibroblast_states ----------------------------------------------------------------

x= studies[1]

clusts= unique(meta.integrated$opt_clust_integrated)
pb.studies.states= lapply(studies, function(x){
  print(x)
  study.state= lapply(clusts,function(y){
    print(y)
    seu2= subset(seu, study==x & opt_clust_integrated == y)
    gex.count= GetAssayData(seu2, assay = "RNA", slot = "count")
    groups= seu2@meta.data[, c("orig.ident")]
    pb <- aggregate.Matrix(t(gex.count),
                           groupings = groups,
                           fun= "sum")
    return(t(pb))

  })
  names(study.state)= clusts
  return(study.state)
})


names(pb.studies.states) = studies

saveRDS(pb.studies.states, file ="output/fib_integration/pseudobulked_study_cellstates.rds")
