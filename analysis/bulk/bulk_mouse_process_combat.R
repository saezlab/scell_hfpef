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

# optional filter out
#x= target%>% filter(factor1== "")%>% pull(sample)
#target = target %>% filter(sample %in% x)
#df= df[,x]

group= target$Group
dge <- DGEList(counts=df, group = group)

#detect and remove low expressed gene
keep <- filterByExpr(dge, min.prop =0.5,min.count = 5 )
#?filterByExpr
dge <- dge[keep,,keep.lib.sizes=FALSE]
#dge$counts= dge$counts + 5

dge <- calcNormFactors(dge)
# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)

### save R-objects
voom_count= v$E

boxplot(voom_count)
boxplot(cpm(dge, log= T))

## combat seq
f1 = factor(group, levels= c("Control", "HFpEF"))
f2 = factor(target$Group, levels= c("Control",
                                   "5-weeks HFpEF",
                                   "10-weeks HFpEF",
                                   "15-weeks HFpEF"
))

design = model.matrix(~f)

counts_adj=ComBat_seq(dge$counts,
                      batch= target$factor1,
                      group= f2)



dge <- DGEList(counts=counts_adj, group = group)

#detect and remove low expressed gene
keep <- filterByExpr(dge, min.prop =0.5,min.count = 5 )
#?filterByExpr
dge <- dge[keep,,keep.lib.sizes=FALSE]
#dge$counts= dge$counts + 5

dge <- calcNormFactors(dge)

## Combat

#Adjust a linear model to the expression data
group= target$Group2
f = factor(group, levels= c("Control", "HFpEF"))
f2 = factor(target$factor1)#levels= c("", "1/2 hearts"))
f = factor(target$Group, levels= c("Control",
                                   "5_weeks",
                                   "10_weeks",
                                   "15_weeks"
                                   ))
design = model.matrix(~f)
design = model.matrix(~Group2, data = target)
#colnames(design) = c("Ct","HFpEF")

ExpMat= ComBat(cpm(dge, log = T),
       batch= target$factor1,
       mod= design)

##PCA
PCA <- prcomp(t(expm) ,center = TRUE, scale. = F)

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

gex.obj= list("meta"= target,
              "exp"= v$E,
              "de"= DE_results
                )
saveRDS(gex.obj, "output/bulk_mouse_obj.rds")



## limma + control for covariate


group= target$Group2
f = factor(group, levels= c("Control", "HFpEF"))
f2 = factor(target$factor1)
f3 = factor(target$Group, levels= c("Control",
                                   "5_weeks",
                                   "10_weeks",
                                   "15_weeks"
))
design = model.matrix(~0+f2+f3)

fit = lmFit(voom_count, design)

#Define contrasts
cont.matrix = makeContrasts(wk5 = "f35_weeks",
                            wk10= "f310_weeks",
                            wk15="f315_weeks",
                            levels=design)

#Empirical Bayes
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

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

expm= limma::removeBatchEffect(cpm(dge, log = T),batch = target$factor1)#, design = design2)
PCA <- prcomp(t(cpm(dge, log = T)) ,center = TRUE, scale. = F)
PCA <- prcomp(t(expm[, colnames(expm) != c("H56")]) ,center = TRUE, scale. = F)
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


## qq


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


gene_sigs= readRDS("output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")


##dea


DE_results %>%
  filter(gene %in% gene_sigs$total$AngII)%>%
  print(n=100)

x= target%>% filter(factor1== "")%>% pull(sample)
