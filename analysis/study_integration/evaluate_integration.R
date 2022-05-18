## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-03-23
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## evaluate batch effects after integration by calculating

integrated_data = readRDS( "output/seu.objs/study_integrations/harmony_fib_filt.rds")
meta= integrated_data@meta.data

cell_dists <- dist(integrated_data@reductions$harmony@cell.embeddings,
                   method = "euclidean")
rm(integrated_data)

## see methods of https://www.nature.com/articles/s41592-021-01336-8#Sec11
calc_batch_mixing = function(batch_col, meta){
  dic= unique(meta[[batch_col]])
  names(dic)= c(1:length(dic))
  meta$cluster_info= match(meta[[batch_col]], dic)
  #check
  meta%>% distinct(cluster_info, !!as.symbol(batch_col))
  si <- cluster::silhouette(meta$cluster_info, cell_dists)

  asw=1-abs( mean(si[, 'sil_width']))
  return(list("sil"= si,
              "batch.asw"= asw))
}

group.asw= calc_batch_mixing("group", meta)
study.asw= calc_batch_mixing("study", meta)
sample.asw= calc_batch_mixing("orig.ident", meta)

group.asw[[2]]
study.asw[[2]]
sample.asw[[2]]

# map do seu obj and plot asw per bell
integrated_data= AddMetaData(integrated_data,
                             metadata = 1-abs(group.asw[1][[1]][, "sil_width"]),
                             col.name = "group.asw")
integrated_data= AddMetaData(integrated_data,
                             metadata = 1-abs(study.asw[1][[1]][, "sil_width"]),
                             col.name = "study.asw")
integrated_data= AddMetaData(integrated_data,
                             metadata = 1-abs(sample.asw[1][[1]][, "sil_width"]),
                             col.name = "sample.asw")

pdf("output/figures/supp/asw.integration.pdf",
    width= 8,
    height=8)

FeaturePlot(integrated_data,
            features =  c("group.asw", "sample.asw", "study.asw"),keep.scale = "all")

VlnPlot(integrated_data,
            features =  c("group.asw", "sample.asw", "study.asw"))
dev.off()
