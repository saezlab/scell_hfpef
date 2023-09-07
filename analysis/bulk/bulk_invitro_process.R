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

library(limma)
library(edgeR)
library(tidyverse)
library(ggrepel)
library(stringr)


target = data.frame("ID"= c(14,15,17,18,21,23), "group"= c("Ct","Ct","Aav", "Aav","Ct","Aav" ))
target = data.frame("ID"= c(14,15,17,18,21,23), "group"= c("Aav","Ct","Aav", "Aav","Ct","Ct" ))
counts= read.csv("data/BulkRNAseq/invitro/invitro_angptl4-rawdata.txt", sep = '\t')
colnames(counts)

df= column_to_rownames(counts, "X")
colnames(df)= target$ID

boxplot(df)
df["Angptl4", ]
target
# optional filter out
x= target%>% filter(factor1== "")%>% pull(sample)
target = target %>% filter(sample %in% x)
df= df[,x]

group= target$group
dge <- DGEList(counts=df, group = group)

#detect and remove low expressed gene
keep <- filterByExpr(dge, min.prop =0.4,min.count = 2 )
keep["Angptl4"]= T
#?filterByExpr
dge <- dge[keep,,keep.lib.sizes=FALSE]
#dge$counts= dge$counts + 5

dge <- calcNormFactors(dge)
# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)

### save R-objects
voom_count= v$E

#Adjust a linear model to the expression data
f = factor(group, levels= c("Ct", "Aav"))
#f2 = factor(target$factor1)#levels= c("", "1/2 hearts"))
design = model.matrix(~0+f)
#colnames(design) = c("Ct","HFpEF")
ExpMat= voom_count
fit = lmFit(ExpMat, design)

#Define contrasts
cont.matrix = makeContrasts(HF = fAav-fCt,
                            levels=design)

#Empirical Bayes
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

#obtain differentially expressed genes
DE_results = as.data.frame(topTable(fit2,adjust.method = "BH",number = Inf)) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(desc(abs(t))) %>%
  as_tibble()
DE_results %>% filter(grepl("Angp", gene))

counts["Angtpl4", ]

DE_results%>%
  mutate(sig= ifelse(P.Value <0.05, "y", "n"),
         hfpef= ifelse(gene %in% gene_sigs$total$HFpEF, gene, ""),
         angii= ifelse(gene %in% gene_sigs$total$AngII, gene, ""))%>%
  ggplot(.,aes(x= logFC, y= -log10(P.Value), col =sig, label = angii))+
  geom_point()+
    geom_label_repel(aes(label = hfpef), max.overlaps = 1000, col = "black")

## qq
hist(DE_results$P.Value)
null_model <- pnorm(rnorm(length(DE_results$gene)))
qq.plot= cbind("null"= sort(null_model), "p.val"= sort(DE_results$P.Value))%>%
  as_tibble()%>%
  ggplot(aes(x= null, y= p.val))+
  geom_point()+
  xlim(c(1,0))+
  ylim(c(1,0))+
  geom_abline(slope = 1, intercept= 0)+
  theme_bw()+
  theme(axis.text = element_text(size = 11, color ="black"),
        rect = element_rect(fill = NA, size = 1),
        legend.position = "none")

qq.plot
gene_sigs= readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")




PCA <- prcomp(t(v$E[,]) ,center = TRUE, scale. = F)

plot.pca = PCA$x %>%
  as.data.frame %>%
  rownames_to_column("ID") %>%
  as_tibble()%>%
  left_join(target%>% mutate(ID = as.character(ID)), by= "ID")

p.pca = ggplot(plot.pca,aes(x= PC1, y= PC2, col= group))+
  geom_point(size= 3)+
  theme_minimal()+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
  ggtitle(paste0(""))+
  geom_text_repel(aes(label= ID),show.legend = FALSE)
p.pca
# p.pca2 = ggplot(plot.pca,aes(x= PC3, y= PC4, col= Group, shape= factor1))+
#   geom_point(size= 3)+
#   theme_minimal()+
#   labs(x= paste0("PC3 (",as.character(round(PCA$sdev[3]^2/sum(PCA$sdev^2)*100)),"%)"),
#        y= paste("PC4 (",as.character(round(PCA$sdev[4]^2/sum(PCA$sdev^2)*100)),"%)"))+
#   ggtitle(paste0(""))+geom_text_repel(aes(label= sample),show.legend = FALSE)
# p.pca2
#
# cowplot::plot_grid(p.pca, p.pca2)
# gex.obj= list("meta"= target,
#               "exp"= v$E,
#               "de"= DE_results
#                 )
# saveRDS(gex.obj, "output/bulk_mouse_obj.rds")
