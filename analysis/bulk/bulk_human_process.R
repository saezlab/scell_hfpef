## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-03-10
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## process DAS et al EBI

library(limma)
library(edgeR)
library(limma)
library(ggrepel)
library(biomaRt)


if(!exists("genemap")){
  genemap= readRDS("~/R-projects/SFB_arrythmia/data/gene_id_symbol_dictionary_July2021.rds")
}

x= read_tsv("data/_Hfpef_bulk_study/HF_pEF_Normal_normalized_EBI_submission_091118.txt")

boxplot(log10(x[,-1]))
## some genes are pseudogenes and will be exclueded

x %>% left_join(genemap %>% rename(GeneID= ensembl_gene_id), by= "GeneID") %>%
  filter(hgnc_symbol== "")%>% dplyr::select(hgnc_symbol,GeneID)

#map gene IDs
df= x %>% left_join(genemap %>% rename(GeneID= ensembl_gene_id), by= "GeneID") %>%
  dplyr::select(hgnc_symbol, everything(), -GeneID) %>% # drop_na() %>%
  group_by(hgnc_symbol)%>%
  summarise(across(everything(), sum)) %>%
  filter(hgnc_symbol!= "")

target= data.frame("sample"= colnames(df), "group"= grepl(x = colnames(df), "Normal"))
target=target %>% mutate(group = ifelse(group, "HFpEF", "Ct"))

df= column_to_rownames(df, "hgnc_symbol")

colnames(df)== target$sample
group= target$group
dge <- DGEList(counts=df, group=group)

#detect and remove low expressed gene
keep <- filterByExpr(dge,group = group, min.prop =0.3,min.count = 2 )
#?filterByExpr
dge <- dge[keep,,keep.lib.sizes=FALSE]
#dge$counts= dge$counts + 5

dge <- calcNormFactors(dge)
# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)

### save R-objects
voom_count= v$E

#Adjust a linear model to the expression data

f = factor(group, levels= c("Ct", "HFpEF"))
design = model.matrix(~0+f)
colnames(design) = c("Ct","HFpEF")
ExpMat= voom_count
fit = lmFit(ExpMat, design)

#Define contrasts
cont.matrix = makeContrasts(HF = HFpEF-Ct,
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
  rownames_to_column("sample") %>%
  as_tibble()%>%
  left_join(target)

p.pca = ggplot(plot.pca,aes(x= PC1, y= PC2,color = group))+
  geom_point(size= 3)+
  theme_minimal()+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
  ggtitle(paste0(""))+
  geom_text_repel(aes(label= sample),show.legend = FALSE)

pdf("output/figures/supp/bulk_dasetal_report.pdf",
    width= 5,
    height= 5)
p.pval.dist
p.volcano
p.pca
dev.off()

cowplot::plot_grid(p.pval.dist, p.volcano, p.pca)

saveRDS(list("DEA"= DE_results,
             "counts"= df,
             "vooom"= v,
             "target"= target
             ), "data/_Hfpef_bulk_study/HF_pEF_processed.rds")

