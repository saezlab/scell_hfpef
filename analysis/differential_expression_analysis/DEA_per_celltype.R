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
## do a pseudobulk DEG analysis per cell type

library(Seurat)
library(tidyverse)
library(Matrix.utils)
library(WriteXLS)
library(ggrepel)

seu= readRDS( file = "output/seu.objs/integrated_cellstate_noECs.rds")

groups= seu@meta.data[, c("orig.ident", "celltype")]

#add group.labels per cell type to meta;
seu@meta.data= seu@meta.data %>% mutate(group= ifelse(grepl("ct", orig.ident), "ct", "hf"),
                                        label = paste0(celltype,".", group))
Idents(seu)= "label"
cells= unique(seu@meta.data$celltype)

# run DEA without downsampling ----------------------------------------------------------------

de.list.no.ds= lapply(cells, function(x){

  print(x)

  DefaultAssay(seu)= "RNA"

  fibs.de= FindMarkers(seu,
                       ident.1 = paste0(x, ".hf"),
                       ident.2 = paste0(x, ".ct"),
                       #max.cells.per.ident = min(cell.count$n),
                       test.use = "MAST",
                       random.seed = 21 ,
                       pseudocount.use = 5,
                       min.pct = 0.1,
                       logfc.threshold = 0.1)
  fibs.de = fibs.de %>%rownames_to_column("gene") %>% as_tibble()
})
names(de.list.no.ds) = cells

map(de.list.no.ds, function(x){
  x= x%>% filter( p_val_adj< 0.05)%>% pull(gene)
  length(x)})


# run DE on contrast for each celltype, with downsampling of cells ----------------------------------------------------------
#1. for each cell type
#2. subsample the same number of cells per sample

de.list= lapply(cells, function(x){

  print(x)

  #calculate cell number per group:
  cell.count= seu@meta.data %>%rownames_to_column("cellid") %>%
  group_by(orig.ident) %>%
    filter(celltype== x) %>%count()

  #select minimum number of cell counts per group to downsample the other
  print(cell.count)
  print(paste0("downsample to ", min(cell.count$n), " cells"))

  DefaultAssay(seu)= "RNA"

  fibs.de= FindMarkers(seu,
              ident.1 = paste0(x, ".hf"),
              ident.2 = paste0(x, ".ct"),
              max.cells.per.ident = min(cell.count$n),
              test.use = "wilcox",
              random.seed = 21 ,
              pseudocount.use = 5,
              min.pct = 0.1,
              logfc.threshold = 0.1)
  fibs.de = fibs.de %>%rownames_to_column("gene") %>% as_tibble()
})

names(de.list) = str_replace_all(cells, "/","")
names(de.list) = str_replace_all(names(de.list), "\\?","")
map(de.list, function(x){hist(x$p_val)})

map(de.list, function(x){
  x= x%>% filter( p_val_adj< 0.05)%>% pull(gene)
  length(x)})


library(WriteXLS)

WriteXLS(de.list, "output/cell_specific_obj/cell_type_based_DE/cells_DEA.xlsx")

saveRDS(de.list, "output/cell_specific_obj/cell_type_based_DE/cells_DEA.rds")

#
# gex.sc= GetAssayData(seu, slot= "counts")
# x= seu[[]]
# p.gex= enframe(gex.sc["Postn", ]) %>% left_join(x%>% rownames_to_column("name")%>% dplyr::select(name, group, orig.ident, celltype))%>%
#   filter(celltype== "fibroblasts")
#
# ggplot(p.gex , aes(x= factor(group), y= value))+
#   geom_boxplot()



# Run DEA multiple times with subsampling  ----------------------------------------------------

get_subsampled_DEG= function(seu, seed, genes.to.test){

  set.seed(seed)

  orig.ids= unique( seu@meta.data$orig.ident)

  cell.count= seu@meta.data %>%rownames_to_column("cellid") %>%
    group_by(orig.ident) %>%
    count()


  random.cell.ids= map(orig.ids , function(x){
    df= seu@meta.data%>% filter(orig.ident== x)
    sample(rownames(df), min(cell.count$n), replace = F)
  })


  #print(cell.count)
  map((random.cell.ids), length)
  print(paste0("downsample to ", min(cell.count$n), " cells"))

  seu.sub= subset(seu, cells = unlist(random.cell.ids))
  DefaultAssay(seu.sub)= "RNA"
  Idents(seu.sub) = "group"
  # DEA= FindMarkers(seu.sub,
  #                      ident.1 = "hf",
  #                      ident.2 = "ct",
  #                     # pseudocount.use = 5,
  #                      #logfc.threshold = 0.1,
  #                      #features =intersect(rownames(seu.sub), genes.to.test),
  #                   #     min.pct = 0,
  #                   # min.cells.feature = 0,
  #                   # min.cells.group = 0,
  #                   # min.diff.pct = 0,
  #                  #min.pct	= 0,
  #                  test.use = "MAST"
  #                       )
  DEA= FindMarkers(seu.sub,
                   ident.1 = "hf",
                   ident.2 = "ct",
                      # max.cells.per.ident = min(cell.count$n),
                       test.use = "wilcox",
                       random.seed = 21 ,
                       pseudocount.use = 5,
                       min.pct = 0.1,
                       logfc.threshold = 0.1)
  print(dim(DEA))
  DEA = DEA %>%rownames_to_column("gene") %>% arrange(gene) %>% as_tibble()
}

intersection <- function(x, y, ...){
  if (missing(...)) intersect(x, y)
  else intersect(x, intersection(y, ...))
}

intersects= map(cells, function(x){

  seu2 = subset(seu, celltype==x)
  #seu2= FindVariableFeatures(seu2)
  #genes.to.test=rownames(HVFInfo(seu2)[VariableFeatures(seu2), ])[1:4000]
  number_of_replicates= 20
  print(x)

  deg.array=  sapply(seq(1:number_of_replicates), function(x){
    get_subsampled_DEG(seu2, seed = x, genes.to.test = rownames(seu2))
  }, simplify = F)

  #filter for significant up and dn regulated genes
  up.= lapply(deg.array, function(x){
    single.df = x
    single.df %>%
      filter(p_val_adj<0.05,
             avg_log2FC> (0),
             pct.1>0.1) %>%
      pull(gene)
    #dn.= single.df %>% filter(p_val_adj<0.01)%>% filter(avg_log2FC<0)%>% pull(gene)
    #list(up., dn.)
  })

  dn.= lapply(deg.array, function(x){
    single.df = x
    single.df %>%
      filter(p_val_adj<0.05,
             avg_log2FC< (0),
             pct.1>0.1) %>%
      pull(gene)
    #dn.= single.df %>% filter(p_val_adj<0.01)%>% filter(avg_log2FC<0)%>% pull(gene)
    #list(up., dn.)
  })

  #report genes in more than 90% of downsamples

  t.up = table(unlist(up.))
  genes.up= names(t.up[t.up > (0.5 *number_of_replicates)])

  t.dn = table(unlist(dn.))
  genes.dn= names(t.dn[t.dn > (0.5 *number_of_replicates)])


  result= list("up"= genes.up,
               "dn"=  genes.dn)


  return(result)
})

names(intersects)= cells
intersects$fibroblasts
length_int= map(intersects, function(x){map(x, length)})

p.DEG.number= enframe(unlist(length_int), value = "number.of.DEGs")%>%
  ggplot(.,aes(x= name, y= number.of.DEGs))+
  geom_point()+
  coord_flip()+
  geom_col(width = 0, colour = "black", lwd = 0.2)+
  theme_minimal()+
  labs(x= "")


saveRDS(intersects, "output/DEG_downsampled_intersect_hfpef.rds")

map(cells, function(x){

   cell.count= seu@meta.data %>%rownames_to_column("cellid") %>%
     filter(celltype ==x)%>%
    group_by(orig.ident) %>%
    count()


})


# compare number of DEGs for different sampling approaches ------------------------------------

df= map(de.list.no.ds, function(x){
  x= x%>% filter( p_val_adj< 0.05)%>% pull(gene)
  length(x)})
df= enframe(df, value = "no.ds")%>% mutate(no.ds= unlist(no.ds))

names(de.list)= cells
df2= map(de.list, function(x){
  x= x%>% filter( p_val_adj< 0.05)%>% pull(gene)
  length(x)})
df2= enframe(df2, value= "ds")%>% mutate(ds= unlist(ds))

df3= map(length_int, function(x) {
  enframe(x) %>% mutate(value = unlist(value))%>% summarise(value= sum(value))
                        })
df3= enframe(df3)%>% unnest(value)%>% rename(repeated.ds= value)
enframe(length_int)
full_join(df, df2, by= "name")%>% full_join(df3, by= "name") %>%
  pivot_longer(cols= c("no.ds", "ds", "repeated.ds"),
               names_to = "approach")%>%
  ggplot(aes(x= name,y= value, fill = approach ))+
  geom_col(position = "dodge")

  intersec

# OLD code ------------------------------------------------------------------------------------


  deg.array=  sapply(seq(1:number_of_replicates), function(x){
    get_subsampled_DEG(seu, seed = x, genes.to.test = genes.to.test)
  }, simplify = F)


  tested.genes= map(deg.array, function(x){
    x$gene
  })
    sort(unique(unlist(tested.genes)))

  ncol(deg.array)
  gene.order= deg.array[1,1]$gene
  ngenes= length(gene.order)
  library(survcomp)

## collect p values per gene and run fisher method for ranking later
fish_p_vec= c()
mean_log_fc= c()
for(i in(1:ngenes)){
  gene_p_vec= c()
  gene_fc= c()
  for(j in (1:number_of_replicates)){
    gene_p_vec= c(gene_p_vec, deg.array[2,j]$p_val[i])
    gene_fc=c(gene_fc, deg.array[3,j]$avg_log2FC[i])

  }
  print(gene_p_vec)
  fish.p = combine.test(p= gene_p_vec, na.rm = T)
  fish_p_vec= c(fish_p_vec, fish.p)
  mean_log_fc= c(mean_log_fc, median(gene_fc, na.rm = T))
}

names(fish_p_vec)= gene.order
fish.p= cbind(enframe(fish_p_vec), mean_log_fc) %>% as_tibble()%>% arrange(value)# %>% print(n=100)

x= 1

#get the interesect of all DE.genes=

up.= lapply(deg.array, function(x){
  single.df = x
  up.= single.df %>% filter(p_val_adj<0.01)%>% filter(avg_log2FC>0.1) %>% pull(gene)
  #dn.= single.df %>% filter(p_val_adj<0.01)%>% filter(avg_log2FC<0)%>% pull(gene)
  #list(up., dn.)
})

up.= lapply(1:ncol(deg.array), function(x){
  single.df = as_tibble(as.data.frame(deg.array[,x]))
  up.= single.df %>% filter(p_val_adj<0.01)%>% filter(avg_log2FC>0) %>% pull(gene)
  #dn.= single.df %>% filter(p_val_adj<0.01)%>% filter(avg_log2FC<0)%>% pull(gene)
  #list(up., dn.)
})

dn.= lapply(1:ncol(deg.array), function(x){
  single.df = as_tibble(as.data.frame(deg.array[,x]))
  #up.= single.df %>% filter(p_val_adj<0.01)%>% filter(avg_log2FC>0) %>% pull(gene)
  dn.= single.df %>% filter(p_val_adj<0.01)%>% filter(avg_log2FC<0)%>% pull(gene)
  #list(up., dn.)
})


do.call(intersection, up.)
do.call(intersection, dn.)
fish.p

de.list= lapply(cells, function(x){

  print(x)

  #calculate cell number per group:
  cell.count= seu@meta.data %>%rownames_to_column("cellid") %>%
    group_by(orig.ident) %>%
    filter(celltype== x) %>%count()

  #select minimum number of cell counts per group to downsample the other
  print(cell.count)
  print(paste0("downsample to ", min(cell.count$n), " cells"))

  DefaultAssay(seu)= "RNA"

  fibs.de= FindMarkers(seu,
                       ident.1 = paste0(x, ".hf"),
                       ident.2 = paste0(x, ".ct"),
                       max.cells.per.ident = min(cell.count$n),
                       random.seed = 21 ,
                       pseudocount.use = 5,
                       logfc.threshold = 0.1)
  fibs.de = fibs.de %>%rownames_to_column("gene") %>% as_tibble()
})

