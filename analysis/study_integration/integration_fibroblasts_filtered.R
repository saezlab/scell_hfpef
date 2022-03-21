## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-10-19
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Integration of fibroblasts between different studies. HFpEF and AngII model
## Rerun the with filtering based on QC
## Plot cell state composition per group per study

library(tidyverse)
library(Seurat)
library(harmony)
library(cluster)
library(clustree)
library(cowplot)
library(purrr)

# # Load SC objects and unifydate meta data: --------------------------------------------------
# Ang II model from circulation:

circ= readRDS( file = "../sc-exploration/output/circ_obj/fibroblasts_seu.rds")
circ@meta.data= circ@meta.data %>%
  mutate(group = ifelse(group== "Ct", "ct", "hf"),
         study= "circ")

spam_cols <- grepl("snn_res",
                   colnames(circ@meta.data))
spam_cols2 <- grepl("seurat",
                    colnames(circ@meta.data))
circ@meta.data <- circ@meta.data[,!spam_cols]
circ@meta.data <- circ@meta.data[,!spam_cols2]
DefaultAssay(circ)= "RNA"

## MI model from forte et al:
forte= readRDS( file = "../sc-exploration/output/cell_mi_files/fibro_subset_integrated.rds")
forte@meta.data= forte@meta.data %>%
  mutate(group = ifelse(OP== "myocardial infarction", "hf", "ct"),
         study= "forte")
DefaultAssay(forte)= "RNA"

## hfpef model ::
hfpef= readRDS(file = "output/cell_specific_obj/fibroblasts.rds")
hfpef@meta.data = hfpef@meta.data %>% mutate(study= "hfpef")
# is hf .ct labeled?
unique(hfpef@meta.data$group)
# hfpef@meta.data = hfpef@meta.data %>% mutate(group= ifelse(group=="hfpef", "hf", "ct"),
#                                              study= "hfpef")
DefaultAssay(hfpef)= "RNA"
#elife= readRDs(file = "../sc-exploration/data/")




# Integration via Harmony ---------------------------------------------------------------------

## 1. Rename cells to facilitate back mapping later:
rename_cells= function(seu, string_add_to_cellid){
  seu = RenameCells(seu,  string_add_to_cellid)
}

hfpef= rename_cells(hfpef, "a")
forte= rename_cells(forte, "b")
circ= rename_cells(circ, "c")

cells_to_exclude = readRDS("output/seu.objs/study_integrations/cellid_to_exclude.rds")

hfpef= subset(hfpef, cells = cells_to_exclude, invert= T)
forte= subset(forte, cells = cells_to_exclude, invert= T)
circ= subset(circ, cells = cells_to_exclude, invert= T)

DefaultAssay(hfpef)= "RNA"
DefaultAssay(forte) = "RNA"
DefaultAssay(circ)= "RNA"

DimPlot(circ)
DimPlot(hfpef)
DimPlot(forte)


# Calculate HVGs  ------------------------------------------------
## Calculate HVGs for each sample in each study, and select those that appear most often

seu.obj= list(circ,
              hfpef,
              forte)
rm(circ)
rm(hfpef)
rm(forte)

calculate_hvgs= function(seu, sample){

  integrated_data <- SplitObject(seu, split.by = sample)

  # First get all the most variable genes per patient per batch

  hvg_list<- map(integrated_data, function(x) {

    DefaultAssay(x) <- "RNA"

    x <- FindVariableFeatures(x, selection.method = "vst",
                              nfeatures = 3000)

    x@assays[["RNA"]]@var.features

  }) %>% unlist()

  hvg_list <- table(hvg_list) %>%
    sort(decreasing = TRUE)

  hist(hvg_list)

  gene_selection= hvg_list[which(hvg_list==max(hvg_list))]
  print(paste0("hvg genes with max number=", length(gene_selection)))
  i= 1

  while(length(gene_selection)<2000){
    gene_selection= c(gene_selection,hvg_list[which(hvg_list==(max(hvg_list)-i))])
    i=i+1
    print(paste0("hvg genes incresed to ", length(gene_selection)))
  }

  return(names(gene_selection))

}

x= map(seu.obj, function(x){calculate_hvgs(x, "orig.ident")})

# These are now the hvg genes per study:
hvg_list <- table(unlist(x)) %>%
  sort(decreasing = TRUE)

hist(hvg_list)
# we will select only those hvg genes that ware hvgs in each study: for 3 studies-> 3
gene_selection= hvg_list[which(hvg_list>1)]
gene_selection= names(gene_selection)

# Read object to subset -----------------------------------------------------

#   # Genes that are stable among batches
# gene_selection <- hvg_list[hvg_list == 2] %>% names()
#
# Re-integrate data

integrated_data <- reduce(seu.obj,
                          merge,
                          merge.data = TRUE)

rm(seu.obj)

#integrated_data

DefaultAssay(integrated_data) = "RNA"

integrated_data@meta.data <- integrated_data@meta.data %>% select(orig.ident, nCount_RNA, nFeature_RNA, percent.mt, group, study)
colnames(integrated_data@meta.data)

## Scale -> PCA -> Harmony
integrated_data <- integrated_data %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(features = gene_selection,
         npcs = 30,
         verbose = FALSE)

# Original PCs
original_pca_plt <- DimPlot(object = integrated_data,
                            reduction = "pca",
                            pt.size = 1,
                            group.by = "orig.ident")

original_study <- DimPlot(object = integrated_data,
                          reduction = "pca",
                          pt.size = 0.4,
                          group.by = "study")

saveRDS(integrated_data, "output/seu.objs/study_integrations/harmony_fib_pre.rds")
integrated_data= readRDS( "output/seu.objs/study_integrations/harmony_fib_pre.rds")
integrated_data@meta.data %>% group_by(study) %>% count
integrated_data@meta.data %>% group_by(group, study) %>% count

# Integrate the data -----------------------

integrated_data <- RunHarmony(integrated_data,
                              group.by.vars= c("study", "orig.ident"),
                              plot_convergence = TRUE,
                              max.iter.harmony = 20,
                              project.dim = F,
                              assay.use= "RNA")

saveRDS(integrated_data, "output/seu.objs/study_integrations/harmony_fib_filt.rds")
integrated_data = readRDS( "output/seu.objs/study_integrations/harmony_fib_filt.rds")

corrected_pca_plt <- DimPlot(object = integrated_data,
                             reduction = "harmony",
                             pt.size =0.4,
                             group.by = "study")

#DefaultAssay(integrated_data)

# RUN UMAP and clustering:

integrated_data <- integrated_data %>%
  RunUMAP(reduction = "harmony", dims = 1:30,
          reduction.name= "umap_harmony_study")%>%
  RunUMAP(reduction = "pca", dims = 1:30,
          reduction.name = "umap_original")


integrated_data <- FindNeighbors(integrated_data,
                                 reduction = "harmony",
                                 dims = 1:30)

saveRDS(integrated_data, "output/seu.objs/study_integrations/harmony_fib_filt.rds")
integrated_data = readRDS( "output/seu.objs/study_integrations/harmony_fib_filt.rds")

# Create new clustering
print("Optimizing clustering")

seq_res <- seq(0.1, 1, 0.1)

integrated_data <- FindClusters(object = integrated_data,
                                resolution = seq_res,
                                verbose = F)
colnames(integrated_data@meta.data)colnames(integrated_data@meta.data)

pls= map(paste0("RNA_snn_res.", seq_res), function(x){
  print(x
        )
  DimPlot(integrated_data, group.by = x)

})

pdf("output/fib_integration/fib_cluster_depth.pdf")
pls
dev.off()
#DimPlot(integrated_data, group.by = "RNA_snn_res.0.2")integrated_data= readRDS(

clustree_plt <- clustree(integrated_data,
                         prefix = paste0(DefaultAssay(integrated_data), "_snn_res."))

# Optimize clustering ------------------------------------------------------
cell_dists <- dist(integrated_data@reductions$harmony@cell.embeddings,
                   method = "euclidean")


cluster_info <- integrated_data@meta.data[,grepl(paste0("RNA_snn_res"),
                                                 colnames(integrated_data@meta.data))] %>%
  dplyr::mutate_all(as.character) %>%
  dplyr::mutate_all(as.numeric)

silhouette_res <- apply(cluster_info, 2, function(x){
  si <- silhouette(x, cell_dists)
  if(!is.na(si)) {
    mean(si[, 'sil_width'])
  } else {
    NA
  }
})

rm(cell_dists)

integrated_data[["opt_clust_integrated"]] <- integrated_data[[names(which(silhouette_res == max(silhouette_res)))]]
#integrated_data[["opt_clust_integrated"]] = integrated_data[["RNA_snn_res.0.5"]]
Idents(integrated_data) = "opt_clust_integrated"

# Reduce meta-data -------------------------------------------------------------------------
spam_cols <- grepl(paste0(DefaultAssay(integrated_data), "_snn_res"),
                   colnames(integrated_data@meta.data))
integrated_data@meta.data <- integrated_data@meta.data[,!spam_cols]


p.harmony.pca= DimPlot(integrated_data, reduction = "harmony")
p.harmony.umap= DimPlot(integrated_data, reduction = "umap_harmony_study", group.by = "opt_clust_integrated",cols= col_vector,
                        label = T,
                        label.size = 4,
                        label.col= "black",
)+
  NoLegend()+
  ggtitle("")+
  labs(x= "UMAP1",
       y= "UMAP2")+
  theme(#axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )

p.harmony.study= DimPlot(integrated_data, reduction = "umap_harmony_study", group.by = "study")
p.harmony.sample= DimPlot(integrated_data, reduction = "umap_harmony_study", group.by = "orig.ident")


# plot proportional features ------------------------------------------------------------------

seu.meta= integrated_data[[]]
source("analysis/utils.R")
#seu.meta= int.fibs[[]]

seu.meta= seu.meta %>% mutate(study= ifelse(study== "forte", "MI", study),
                    study= ifelse(study== "circ", "AngII", study),
                    study= ifelse(study== "hfpef", "HFpEF", study))
# count total cells per group:
fibs.counts= seu.meta%>% group_by(study) %>% count()

# get proportions of study per group
df.prop= seu.meta%>% group_by(opt_clust_integrated, study) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  left_join(fibs.counts %>%
              rename(n.all= n),
            by= "study") %>%
  mutate(freq2= n/n.all) #  this is proportion of cluster within a study

p.int.fib.prop= ggplot(df.prop, aes(x= study, y= freq2, fill =opt_clust_integrated))+
  geom_col()+
  scale_fill_manual(values = col_vector)+
  theme_minimal()+
  labs(y= "cell state proportion per study (%)",
       x= "",
       fill= "integrated cluster")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

# now scale the proportion of study per cluster towards total cells per study
df.prop= df.prop %>% mutate(all.freqs= sum(freq2)) %>% ungroup() %>% mutate(sacaled.prop= freq2/all.freqs)

p.proportion.sclaed= ggplot(df.prop, aes(x= opt_clust_integrated, y= sacaled.prop, fill =study))+
  geom_col()+
  geom_hline(yintercept = c(1/3, 2/3))+
  theme_minimal()

# now plot how many hf-ct samples per cluster
df.prop.group= seu.meta%>% group_by(opt_clust_integrated, group) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

p.group.per.clust= df.prop.group %>%
  ggplot(., aes(x= opt_clust_integrated, y= freq, fill =group))+
  geom_col()+
  #scale_fill_manual(values = col_vector)+
  theme_minimal()

# plot per study the proportions of each cell state between groups
pls= map(unique(seu.meta$study), function(x){
  seu.meta.x= seu.meta%>% filter(study==x)
  fib.counts= seu.meta.x%>% group_by(opt_clust_integrated) %>% count()
  total_prop_per_cluster= seu.meta.x%>%
    group_by(opt_clust_integrated, group) %>%
    summarise(n = n())%>%
    left_join(fib.counts%>% rename(n.all= n),by="opt_clust_integrated")%>%
    mutate(cluster.prop= n/n.all)
  p.total= total_prop_per_cluster %>%
    ggplot(., aes(x= opt_clust_integrated, y= cluster.prop, fill =group))+
    geom_col()+
    #scale_fill_manual(values = col_vecto+r)+
    theme_minimal()+
    ggtitle(x)

  df= calc_mean_proportions_per_group(seu.meta = seu.meta.x,
                                      cluster.col = "opt_clust_integrated",
                                      group.col = "group",
                                      sample.col ="orig.ident"

                              )
 p.col=  df$groupwise %>%
    ggplot(., aes(x= opt_clust_integrated, y= mean.percent, fill =group))+
    geom_col(position = "dodge")+
    scale_fill_manual(values = hfpef_cols)+
    theme_minimal()+
   theme(axis.title.x = element_blank())+
    ggtitle(x)

  p.col
  #p.total


})

p.props.group= cowplot::plot_grid(pls[[2]]+theme(legend.position = "none"),
                                  pls[[1]]+theme(legend.position = "none"),
                                  pls[[3]]+theme(legend.position = "none"),
                                  nrow= 1)

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  pls[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12))
)

p.props.group= plot_grid(p.props.group, legend, rel_widths = c(3, .4))

p.props.group

p.marker = VlnPlot(integrated_data, c("Postn", "Meox1", "Cilp", "Igfbp3", "Dcn", "Col1a1","Fap"))

p.study=  calc_props(seu.meta, cluster.col = "opt_clust_integrated", "study")
p.sample=  calc_props(seu.meta, cluster.col = "opt_clust_integrated", "orig.ident")

pdf("output/figures/integration_studies/fibroblast_integration_report.pdf")
original_pca_plt
original_study
p.harmony.pca
p.harmony.sample
p.harmony.study
#clustree_plt
p.harmony.umap
p.study
p.sample
p.int.fib.prop
p.proportion.sclaed
p.group.per.clust
p.props.group
p.marker
dev.off()

# plot single pdfs:

pdf("output/figures/integration_studies/integrated_umap.pdf",
    width= 5,
    heigh= 5)
p.harmony.umap
dev.off()


pdf("output/figures/integration_studies/integrated_cluster_comp.pdf",
    width= 4,
    heigh= 4)
p.int.fib.prop
dev.off()

pdf("output/figures/integration_studies/integrated_cluster_prop_ss.pdf",
    width= 6,
    height = 2)
p.props.group
dev.off()

saveRDS(integrated_data, "output/seu.objs/study_integrations/harmony_fib_filt.rds")
integrated_data= readRDS( "output/seu.objs/study_integrations/harmony_fib_filt.rds")
