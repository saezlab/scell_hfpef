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
## do a pseudobulk DEG analysis of the 2x2 samples


library(Seurat)
library(tidyverse)
library(Matrix.utils)

seu= readRDS( file = "output/seu.objs/integrated_cellstate_noECs.rds")

groups= seu@meta.data[, c("orig.ident")]

gex.count= GetAssayData(seu, slot= "counts")
head(gex.count)


# Pseudobulking -------------------------------------------------------------------------------
## Pseudobulking:
pb <- aggregate.Matrix(t(gex.count),
                       groupings = groups,
                       fun= "sum")
dim(pb)

colnames(pb)
rownames(pb)

# LIMMA on pseudobulks ------------------------------------------------------------------------

library(edgeR)
library(limma)

#subset pb object to desired celltype

df= t(pb)


#### Filter & Normalize
#create DGE class object
group <- c("ct","ct","hf","hf")
print("make sure this aligns:")
print(group)
print(colnames(df))
dge <- DGEList(counts=df, group=group)

#detect and remove low expressed gene
keep <- filterByExpr(dge,group = group ,min.total.count= 30, min.count= 10, min.prop= 1)
#?filterByExpr
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)

#
voom_count= v$E

## Run PCA:
PCA <- prcomp(t(voom_count) ,center = TRUE, scale. = T)

plot.pca = PCA$x %>%
  as.data.frame %>%
  rownames_to_column("sample") %>%
  as_tibble()%>%
  mutate(group = ifelse((grepl("ct", sample)), "Ctrl", "HFpEF"))

p.pca = ggplot(plot.pca,aes(x= PC1, y= PC2,color = group))+
  geom_point(size= 3)+
  geom_label(aes(label= sample))+
  theme_minimal()+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
  ggtitle(paste0("pseudobulked ", " sample"))

p.pca

pdf("output/figures/DEA/full_pseudobulk_PCA.pdf")
p.pca
dev.off()

#Adjust a linear model to the expression data
f = factor(group, levels= c("ct","hf"))
design = model.matrix(~0+f)
colnames(design) = c("ct","hf")
fit = lmFit(voom_count, design)

#Define contrasts
cont.matrix = makeContrasts(hf_ct = hf-ct,
                            levels=design)

#Empirical Bayes
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

#obtain differentially expressed genes
DE_results = as.data.frame(topTable(fit2,adjust.method = "BH",number = Inf)) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(desc(abs(t))) %>%
  as_tibble()

DE_results  %>% print(n=100)

p.pval.dist= DE_results %>% ggplot(., aes(x= P.Value))+geom_histogram(bins = 100)
DE_results %>% ggplot(., aes(x= adj.P.Val))+geom_histogram()
p.volcano= DE_results%>% ggplot(., aes(x= logFC, y= -log10(P.Value)))+
  geom_point()

