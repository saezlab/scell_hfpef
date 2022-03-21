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
## 1. Calculate Pseudobulk profiles per sample and per cell type.
## 2. Normalize profiles with TMM (edgeR), perform PCA and limma
## 3. Calculate cosine distances between and within groups and plot

library(Seurat)
library(tidyverse)
library(Matrix.utils)
library(WriteXLS)
library(ggrepel)
library(edgeR)
library(limma)

seu= readRDS( file = "output/seu.objs/integrated_cellstate_nameupdated.rds")

seu@meta.data= seu@meta.data%>% mutate(celltype= ifelse(grepl("Fib", celltype2), "Fibroblasts", celltype2))
groups= seu@meta.data[, c("orig.ident", "celltype")]

gex.count= GetAssayData(seu, slot= "counts")
head(gex.count)

seu<-FindVariableFeatures(seu, assay = "RNA")
# 1. Pseudobulking -------------------------------------------------------------------------------
## Pseudobulking:
pb <- aggregate.Matrix(t(gex.count),
                       groupings = groups,
                       fun= "sum")
dim(pb)

colnames(pb)
rownames(pb)
splitf <- sapply(stringr::str_split(rownames(pb),
                                    pattern = "_",
                                    n = 2),
                 `[`, 1)

pb <- split.data.frame(pb,
                       factor(splitf)) %>%
  lapply(function(u)
    magrittr::set_colnames(t(u),
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))


saveRDS(pb, "output/cell_specific_obj/cell_type_based_DE/samplewise_cellwise_PB.rds")
pb = readRDS("output/cell_specific_obj/cell_type_based_DE/samplewise_cellwise_PB.rds")


# 2. TMM, LIMMA and PCA on pseudobulks ------------------------------------------------------------------------

#celltype = "fibroblasts"
#pb= pb$ct1
run_TMM_limma= function(pb, celltype){


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
  DE_results %>% ggplot(., aes(x= adj.P.Val))+geom_histogram()
  p.volcano= DE_results%>% ggplot(., aes(x= logFC, y= -log10(P.Value)))+
    geom_point()

  #PCA

  PCA <- prcomp(t(voom_count) ,center = TRUE, scale. = T)

  plot.pca = PCA$x %>%
    as.data.frame %>%
    rownames_to_column("sample") %>%
    as_tibble()%>%
    mutate(group = ifelse((grepl("ct", sample)), "Ctrl", "HFpEf"))

  p.pca = ggplot(plot.pca,aes(x= PC1, y= PC2,color = group))+
    geom_point(size= 3)+
    theme_minimal()+
    labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
         y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
    ggtitle(paste0("pseudobulked ", celltype))+
    geom_text_repel(aes(label= sample),show.legend = FALSE)


  return(list("p" = p.volcano,
              "dea"= DE_results,
              "p2"= p.pca,
              "pca.df"= PCA,
              "pca.pca"= plot.pca,
              "dge"= dge,
              "voom"= v$E
              )
         )

}

#
# PCA= pb.de$Smooth$pca.df
# cells= colnames(pb$ct1)
#cells= c("Endothelial", "Fibroblasts", "T", "B","NK","SMC", "Granulocytes", "Macrophages")
cells= colnames(pb$ct1)
pb.de = lapply(cells, function(x){
  print(x)
  pb.dea= run_TMM_limma(pb,celltype =  x)
})

names(pb.de) = c("B.cells", "Cd4.T.cell", "Cd8.T.cell", "Endothelial","Fibroblasts", "Granulocytes",  "Macrophages", "NK.cells", "SMC/Pericytes")

saveRDS(object = pb.de,
        file= "output/cell_specific_obj/cell_type_based_DE/pseudobulk_celltype_DE.rds")

pb.de= readRDS(file= "output/cell_specific_obj/cell_type_based_DE/pseudobulk_celltype_DE.rds")

pls.pca= map(pb.de, function(x) (x$p2))
p.PCA = cowplot::plot_grid(plotlist = pls.pca)

pls.volc= map(pb.de, function(x) (x$p))
p.volc= cowplot::plot_grid(plotlist = pls.volc)

pls. = map(names(pb.de), function(x){
  pb.de[[x]]$dea%>%
    ggplot(aes(x= P.Value))+
    geom_histogram(bins= 80)+
    ggtitle(x)
})

p.hist= cowplot::plot_grid(plotlist= pls.)
pdf("output/figures/DEA/pseudobulk.overview.pdf",
    height= 12,
    width= 12)
p.PCA
p.volc
p.hist
dev.off()


