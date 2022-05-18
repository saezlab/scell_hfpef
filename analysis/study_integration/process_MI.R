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
## processing of MI study Forte et. al



directory= "/net/data.isilon/ag-saez/bq_jlanzer/sc-mi"

file= as_tibble(readRDS(file.path(directory, "data/cell_paper/full_metainfo.rds")))

#filepaths in cluster:
filepaths <- list.files(path=file.path(directory, "data/cell_paper"), pattern="*identity.Rdata", full.names=T, recursive=T)

#update the meta=
file$file.path = filepaths

## select the sample for integration:
# 17004 and 17008 are two runs. we will only select one , based on QC metrics.
# MF17004_GT17-10999 has more cells than MF17004_GT17-13609
# MF17008_GT17-11003 has more cells than MF17008_GT17-13613
# remove 129S1 strain

file2 = file %>% filter(Characteristics.strain. == "C57BL/6J",
                        !filename.long %in% c("MF17008_GT17-13613","MF17004_GT17-13609" ))


# create object list for seurat input
obj_list= lapply(c(1:dim(file2)[1]), function(i){
  load(file2$file.path[i], X <- new.env())
  assign(file2$filename.long[i], X$data)
})

# seurat integration pipeline -----------------------------------------------------------------

anchors= FindIntegrationAnchors(object.list = obj_list, dims= 1:20)

int_cell_MI<- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(int_cell_MI)  = "integrated"

int_cell_MI= ScaleData(int_cell_MI, verbose = F )

#add meta data to seurat integrated object
x= int_cell_MI[[]]  %>% left_join(file %>% rename(orig.ident= filename.long) )
int_cell_MI = AddMetaData(int_cell_MI,x$filename, col.name = "filename")
int_cell_MI = AddMetaData(int_cell_MI,x$Factor.Value.time., col.name = "time")
int_cell_MI = AddMetaData(int_cell_MI,x$Factor.Value.injury., col.name = "OP")

int_cell_MI <- RunPCA(int_cell_MI, npcs = 30, verbose = FALSE)

int_cell_MI <- RunUMAP(int_cell_MI, reduction = "pca", dims = 1:20)
int_cell_MI <- FindNeighbors(int_cell_MI, reduction = "pca", dims = 1:20)
int_cell_MI <- FindClusters(int_cell_MI, resolution = 0.5)

saveRDS(int_cell_MI, file =filepath(directory, "output/integrated_seurat.rds"))


# find fibros ---------------------------------------------------------------------------------

library(cowplot)
p1= DimPlot(int_cell_MI, cols = col.set)
p2= FeaturePlot(int_cell_MI, features= c("Gsn", "Col1a1", "Pdgfra"), ncol= 3)
p3= DotPlot(int_cell_MI, features= c("Gsn", "Col1a1", "Pdgfra"))
p.c= plot_grid(plot_grid(p1, p3), p2, ncol = 1)
pdf("output/figures/fibro_subset_MI.pdf",
    width= 10,
    height = 8)
p.c
dev.off()

directory= "/home/jan/R-projects/sc-exploration/"
int_cell_MI = readRDS(file= paste0(directory, "output/cell_mi_files/integrateint_cell_MI.rds"))

#quick look up of t-cells
# DimPlot(int_cell_MI)
# FeaturePlot(int_cell_MI, features = c("Cd4","Cd8a", "Ctla4", "Tcf7", "Cd28"))
# int_cell_MI[[]]



colnames(int_cell_MI[[]])

#subset to cluster 0,5,6
fibros= subset(int_cell_MI, subset= integrated_snn_res.0.1 %in% c(0,5,6))

saveRDS(fibros, file= "output/seu.objs/MI_fibro_seu.rds")
