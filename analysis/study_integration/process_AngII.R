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
## Processing of AngII model (E-MTAB-8810). FASTQ files were preprocessed via cell ranger


# libs and data -------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)

#remote =
#directory = "/net/data.isilon/ag-saez/bq_jlanzer/sc-mi/"

#local =
directory= "/home/jan/R-projects/sc-exploration/"

# data, create meta_file information with paths and names

# filepath
filepath <- list.files(path= paste0(directory, "data/circ_MI"), pattern="*identity.Rdata", full.names=T, recursive=T)

# filenames(same as in meta)
filenames= list.files(path= paste0(directory, "data/circ_MI/"), pattern="*identity.Rdata", full.names=F, recursive=T)

# extract the long and short version of the file name (sample name)
m <- regexpr(".+?(?=_)", filenames, perl=TRUE)
filenames.short= regmatches( filenames, m)

m2 <- regexpr("^[^\\/]*", filenames)
filename.long = regmatches( filenames, m2)

# store all filename data in df:
file_meta = data.frame(filename.short= filenames.short,
                  filename.long = filename.long,
                  file.path = filepath)

print(file_meta)
saveRDS(file_meta, paste0(directory, "data/circ_MI/file_meta.rds"))

####
raw_meta=as.tibble( read.csv("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8810/E-MTAB-8810.sdrf.txt", sep = '\t'))

colnames(raw_meta)

meta = raw_meta %>% dplyr::select(Extract.Name,
                                  Factor.Value.compound.,
                                  Factor.Value.dose.,
                                  Unit.dose.unit.,
                                  Characteristics.strain.,
                                  Characteristics.age.,
                                  Characteristics.genotype.,
                                  Characteristics.sex.)

write.csv(meta, file ="data/circ_MI/meta_circMI.csv")

#### load object list from meta.file paths
# # create object list for seurat input
obj_list= lapply( seq(1,dim(file_meta)[1]), function(i){
  print(i)
  print("loading")
  print(file_meta$file.path[i])
  load(file_meta$file.path[i], X <- new.env())
   assign(file_meta$filename.long[i], X$data)
 })



# Integration ---------------------------------------------------------------------------------

anchors= FindIntegrationAnchors(object.list = obj_list, dims= 1:20)
saveRDS(anchors, paste0(directory, "data/circ_MI/anchros.rds"))
int_seu<- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(int_seu)  = "integrated"

int_seu= ScaleData(int_seu, verbose = F )

# run cluster on integrated -------------------------------------------------------------------

int_seu <- RunUMAP(int_seu, reduction = "pca", dims = 1:20)
int_seu <- FindNeighbors(int_seu, reduction = "pca", dims = 1:20)
int_seu <- FindClusters(int_seu, resolution = 0.5)

saveRDS(inte_seu, paste0(directory, "data/circ_MI/integrated_seurat_circMI.rds"))


# identify fibroblasts ---------------------------------------------------------------------------------------------

meta = read_csv(paste0(directory, "data/circ_MI/meta_circMI.csv"))
files = readRDS(paste0(directory, "data/circ_MI/file_meta.rds"))

meta= meta %>% left_join(files %>% select(-file.path) %>% rename(Extract.Name = filename.short))

seu_meta= seu[[]]
seu_meta= seu_meta %>%
  left_join(meta %>% rename(orig.ident= filename.long ))#%>%
  mutate(group= ifelse(treatment == "angiotensin II", "angII", "Ct"))

# add info to seu, mainly treatment info
seu = AddMetaData(seu, seu_meta$Factor.Value.compound., col.name = "treatment")
seu = AddMetaData(seu,seu_meta$Factor.Value.dose., col.name = "treatment_dose")
#seu = AddMetaData(seu,seu_meta$group, col.name = "group")


##
DimPlot(seu)
DimPlot(seu, group.by = "treatment")

#rerun  a coarse clustering to identify major celltype
DefaultAssay(seu) = "integrated"

seu <-RunPCA(seu,  npcs = 23, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:23)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:23)
seu <- FindClusters(seu, resolution = seq(0.1,0.5,0.1))

Idents(seu) =  "integrated_snn_res.0.1"


#seu = readRDS("data/circ_MI/integrated_seurat_circMI_processed.rds")

library(cowplot)
p1= DimPlot(seu, cols = col.set)
p2= FeaturePlot(seu, features= c("Gsn", "Col1a1", "Pdgfra"), ncol= 3)
p3= DotPlot(seu, features= c("Gsn", "Col1a1", "Pdgfra"))

p.c= plot_grid(plot_grid(p1, p3), p2, ncol = 1)
pdf("output/figures/fibro_subset_AngII.pdf",
    width= 10,
    height = 8)
p.c
dev.off()
# subset to fibroblasts -----------------------------------------------------------------------

fibro= subset(seu,  subset= integrated_snn_res.0.1 %in% c(0,8))

DimPlot(fibro)

fibro <-RunPCA(fibro,  npcs = 23, verbose = FALSE)
fibro <- RunUMAP(fibro, reduction = "pca", dims = 1:23)
fibro <- FindNeighbors(fibro, reduction = "pca", dims = 1:23)
fibro <- FindClusters(fibro, resolution = seq(0.1,0.5,0.1))

levels= c("integrated_snn_res.0.1" ,"integrated_snn_res.0.2", "integrated_snn_res.0.3", "integrated_snn_res.0.4","integrated_snn_res.0.5")
DimPlot(fibro, group.by = levels)

Idents(fibro) =  "integrated_snn_res.0.5"


saveRDS(fibro, file = "output/circ_obj/fibroblasts_seu.rds")

