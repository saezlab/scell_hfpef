## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-03-14
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Calculate distance between samples

library(tidyverse)
library(Seurat)
library(edgeR)
library(ComplexHeatmap)
library(lsa)

seu = readRDS("output/seu.objs/integrated_cellstate_nameupdated.rds")
seu@meta.data= seu@meta.data %>% mutate(celltype= ifelse(grepl("Fib", celltype2), "Fibroblasts", celltype2))
seu.int= readRDS("output/seu.objs/study_integrations/harmony_fib_filt.rds")

pb.de= readRDS(file= "output/cell_specific_obj/cell_type_based_DE/pseudobulk_celltype_DE.rds")

pb.de.study= readRDS(
        "output/fib_integration/pseudobulked.per.study.sample.rds")

pb.studies.states=readRDS( "output/fib_integration/pseudobulked_study_cellstates.rds")

# run for cell type ---------------------------------------------------------------------------------------
#takes gex profiles and returns ratio of within/between sample distances

calc_distances= function(cells, pb.de, ratios= F ){

  require(lsa)

  # map pseudobulks per celltype to calculate cosine distances between samples
  # and return ratio of within and between groups
  res= map(cells, function(x){
    print(x)
    #TMM normalization (cpm with edgeR)
    pb.x= pb.de[[x]]$dge
    pb.x.cpm= cpm(pb.x, log = F)


    #get hvg:
    seu2= subset(seu, celltype== x)
    seu2= FindVariableFeatures(seu2 )

    #reduce pseudobulk to hvg
    pb.x.cpm= pb.x.cpm[rownames(pb.x.cpm) %in% rownames(HVFInfo(seu2))[1:3000], ]

    #eucledian:
    #cpm.dist= as.matrix(dist(t(pb.x.cpm)))

    #cosine similarity:
    cpm.dist= cosine((pb.x.cpm))
    cpm.dist= 1-cpm.dist


    # calculate mean distance within groups:
    within.= mean(c(cpm.dist["hf2","hf1"],
                    cpm.dist["ct1","ct2"]))

    if(!ratios){

      #calculate median distance between groups:
      across.= median(c(cpm.dist["hf1", "ct1"],
                        cpm.dist["hf1", "ct2"],
                        cpm.dist["hf2", "ct1"],
                        cpm.dist["hf2", "ct2"]))
      return(across./within.)

    }else{
      ratios= map(c(cpm.dist["hf1", "ct1"],
                       cpm.dist["hf1", "ct2"],
                       cpm.dist["hf2", "ct1"],
                       cpm.dist["hf2", "ct2"]),
                   function(x){within./x})

      return(ratios)
    }



  })

  names(res)= cells
  return(res)
}

names(pb.de)
names(pb.de) = str_replace_all(names(pb.de), "Cd", "CD")
names(pb.de) = str_replace_all(names(pb.de), "T.cell", "T.cells")
cells=names(pb.de)
cos.p= calc_distances(cells, pb.de)


cos.p2= calc_distances(cells, pb.de, ratios= T)
#cos.p= calc_distances(cells, pb.de, ratios = T)
# if ratios were calculated, use this function to perform one sample wilcox
ps= map(cos.p2, function(x){

  print(unlist(x))
  print(mean(unlist(x)))
  print(median(unlist(x)))
  wilcox.test(unlist(x), mu = 1)$p.value
})
ps

# plot median ratio:
p.distance_ratio= enframe(cos.p) %>%
  rowwise()%>%
  mutate(value= mean(unlist(value)))%>%
  mutate(name  = factor(name,levels = rev(c("Fibroblasts",
                                            "Endothelial",
                                            "NK.cells",
                                            "Macrophages",
                                            "CD4.T.cells",
                                            "CD8.T.cells",
                                            "B.cells",
                                            "Granulocytes",
                                            "SMC/Pericytes")))
         )%>%
  ggplot(., aes(x= name,y= value ))+
  geom_col(width = 0, colour = "black", lwd = 1)+
  geom_point()+
  #geom_col(width = 0, colour = "black", lwd = 2)+
  labs(y= "distance ratio (between/within)", x="")+
  theme_minimal()+
  coord_flip()+
  geom_hline(yintercept = 1, col= "grey")+
  theme(axis.text = element_text(size= 11, colour= "black"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

saveRDS(p.distance_ratio,"output/figures/main/p.distance_ratio_obj.rds")
p.distance_ratio
  coord_flip()

# heatmap?
enframe(cos.p) %>%
  rowwise()%>%
  mutate(value= mean(unlist(value)),
         dummy= "1")%>%
  ggplot(.,aes(x= dummy, y= name, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low= "green", high ="red", mid ="white")+
  theme()

pdf("output/figures/sample_distances_celltype.pdf",
    width= 5,
    height= 2.4)
p.distance_ratio
dev.off()


# get.distane.matrix --------------------------------------------------------------------------

calc_distances2= function(cells, pb.de, ratios= F ){

  require(lsa)

  # map pseudobulks per celltype to calculate cosine distances between samples
  # and return ratio of within and between groups
  res= map(cells, function(x){
    print(x)
    #TMM normalization (cpm with edgeR)
    pb.x= pb.de[[x]]$dge
    pb.x.cpm= cpm(pb.x, log = F)

    #get hvg:
    seu2= subset(seu, celltype== x)
    seu2= FindVariableFeatures(seu2 )

    #reduce pseudobulk to hvg
    pb.x.cpm= pb.x.cpm[rownames(pb.x.cpm) %in% rownames(HVFInfo(seu2))[1:3000], ]

    #eucledian:
    #cpm.dist= as.matrix(dist(t(pb.x.cpm)))

    #cosine similarity:
    cpm.dist= cosine((pb.x.cpm))
    cpm.dist= 1-cpm.dist

  })

  names(res)= cells
  return(res)
}
cos.p2= calc_distances2(cells, pb.de)

map(names(cos.p2), function(x){Heatmap(cos.p2[[x]], name = x)})
# study ---------------------------------------------------------------------------------------


calc_distances_study= function(studies,seu.int,  pb.de.study, ratios= F ){

    require(lsa)

    # map pseudobulks per celltype to calculate cosine distances between samples
    # and return ratio of within and between groups
    res= map(studies, function(x){
      print(x)
      #TMM normalization (cpm with edgeR)
      pb.x= pb.de.study$pseudobulkedGEX[[x]]
      dge <- DGEList(counts=pb.x)

      #detect and remove low expressed gene
      keep <- filterByExpr(dge, min.prop =0.5,min.count = 5 )
      #?filterByExpr
      dge <- dge[keep,,keep.lib.sizes=FALSE]
      dge$counts= dge$counts + 5

      dge <- calcNormFactors(dge)

      pb.x.cpm= cpm(dge, log = F)

      #get hvg:
      seu2= subset(seu.int, study== x)
      seu2= FindVariableFeatures(seu2 )

      #reduce pseudobulk to hvg
      pb.x.cpm= pb.x.cpm[rownames(pb.x.cpm) %in% rownames(HVFInfo(seu2))[1:3000], ]
      rm(seu2)

      #eucledian:
      #cpm.dist= as.matrix(dist(t(pb.x.cpm)))

      #cosine similarity:
      cpm.dist= cosine((pb.x.cpm))
      cpm.dist= 1-cpm.dist

      #get group.data
      group.df= seu2@meta.data%>% distinct(orig.ident, group)
      ct.= group.df%>% filter(group=="ct")%>% pull(orig.ident)
      hf. = group.df%>% filter(group=="hf")%>% pull(orig.ident)

      #within group.
      cpm.dist2[upper.tri(cpm.dist, diag= T)]<-NA

      within.u= mean(mean(cpm.dist[hf., hf.], na.rm = T), mean(cpm.dist[ct., ct.], na.rm = T))

      across.u= mean( mean(cpm.dist[hf., ct.], na.rm = T)
                      # mean(cpm.dist[ct., hf.], na.rm = T)
                       )
      boxplot(c(cpm.dist[hf., hf.],cpm.dist[ct., ct.]),
              c(cpm.dist[hf., ct.], cpm.dist[ct., hf.])
      )
      return(across.u/within.u)
      # calculate mean distance within groups:


      # if(!ratios){
      #
      #   #calculate median distance between groups:
      #   across.= median(c(cpm.dist["hf1", "ct1"],
      #                     cpm.dist["hf1", "ct2"],
      #                     cpm.dist["hf2", "ct1"],
      #                     cpm.dist["hf2", "ct2"]))
      #   return(across./within.)
      #
      # }else{
      #   ratios= map(c(cpm.dist["hf1", "ct1"],
      #                 cpm.dist["hf1", "ct2"],
      #                 cpm.dist["hf2", "ct1"],
      #                 cpm.dist["hf2", "ct2"]),
      #               function(x){within./x})
      #
      #   return(ratios)
      # }



    })

    names(res)= studies
    return(res)
}

calc_distances_study(studies= unique(seu.int$study),
                      seu.int= seu.int,
                     pb.de.study)



#study state


require(lsa)

# map pseudobulks per celltype to calculate cosine distances between samples

studies= names(pb.studies.states)
cellstates= names(pb.studies.states$circ)
 x="circ"
 y= "0"

res= lapply(studies, function(x){
  print(x)
  lapply(cellstates, function(y){
    print(x)
    #TMM normalization (cpm with edgeR)
    pb.x= pb.studies.states[[x]][[y]]
    dge <- DGEList(counts=pb.x)

    #detect and remove low expressed gene
    keep <- filterByExpr(dge, min.prop =0.5,min.count = 5 )
    #?filterByExpr
    dge <- dge[keep,,keep.lib.sizes=FALSE]
    dge$counts= dge$counts + 5

    dge <- calcNormFactors(dge)

    pb.x.cpm= cpm(dge, log = F)
    #pb.x.cpm= cpm(pb.studies.states[[x]][[y]], log = F)

    #get hvg:
    seu2= subset(seu.int, study==x & opt_clust_integrated == y)
    seu2= FindVariableFeatures(seu2 )

    pb.x.cpm= pb.x.cpm[rownames(pb.x.cpm) %in% rownames(HVFInfo(seu2))[1:3000], ]

    cpm.dist= cosine((pb.x.cpm))
    cpm.dist= 1-cpm.dist

    #get group.data
    group.df= seu2@meta.data%>% distinct(orig.ident, group)
    ct.= group.df%>% filter(group=="ct")%>% pull(orig.ident)
    hf. = group.df%>% filter(group=="hf")%>% pull(orig.ident)

    #within group.
    cpm.dist2= cpm.dist
    cpm.dist2[upper.tri(cpm.dist, diag= T)]<-NA

    #gather as df.
    #within

    df.1= enframe(c(cpm.dist2[hf., hf.], cpm.dist2[ct., ct.]))%>% mutate(name= c(rep("hf-hf", length(c(cpm.dist2[hf., hf.]))),
                                                                           rep("ct-ct", length(c(cpm.dist2[ct., ct.]))))
                                                                   )%>%
      drop_na

    #across
    df.2=   enframe(c(cpm.dist[hf., ct.])) %>% mutate(name="hf-ct")

    return(rbind(df.1, df.2)%>% mutate(cluster= y, study= x))

  }) %>% do.call(rbind, .)
})%>% do.call(rbind, res) %>%
  mutate(cross= ifelse(name == "hf-ct", "across", "within"))

#plot results:

#ratio of medians=
p.state.distances= res%>%
  group_by(study,cluster,  cross ) %>%
  summarise("m.dist"= median(value))%>%
  pivot_wider(names_from = cross, values_from= m.dist)%>%
  mutate(dist.ratio= across/within)%>%
  ggplot(aes(x=study, y= cluster, fill = log2(dist.ratio)))+
  geom_tile()+
  scale_fill_gradient(low= "white", high= "red")+
  geom_text(aes(label= round(dist.ratio,1)))


p.dists= res%>% ggplot(., aes(x= cluster, y= value, fill = cross))+
  facet_grid(rows= vars(study), scales = "free_y")+
  geom_boxplot()

pdf(file= "output/figures/integration_studies/sample.distances.cellstates.pdf",
    width= 5, height= 5)
p.state.distances
p.dists
dev.off()


# as suggested by Reviewer#3 perform correlation of biol. replicats ---------------------------
library(ComplexHeatmap)
pb.df= pb.de.study$pseudobulkedGEX$forte
target= seu.int[[]]%>%
  as_tibble()%>%
  distinct(study, orig.ident, group)

plot_corrs= function(pb.df, study.name){

  #pb.df = pb.de$Endothelial$voom
  l1= length(colnames(pb.df))

  df= data.frame(matrix(nrow = l1, ncol= l1),row.names = colnames(pb.df))
  colnames(df)= colnames(pb.df)
  for( x in rownames(df)){
    for( y in rownames(df)){

      df[x,y]= cor.test(as.matrix(pb.df)[,x], pb.df[,y], method = "pearson")$estimate
    }

  }

  target2=
    target %>% filter(study== study.name)
  x= target2$group
  names(x) = target2$orig.ident

  ha = HeatmapAnnotation(group= x, col = list(group= c("ct"= "grey", "hf"= "darkgreen")))
  #ha2 = rowAnnotation(group= x, row = list(group= c("ct"= "grey", "hf"= "darkgreen")))
  Heatmap(df, top_annotation = ha)
}
p1= plot_corrs(pb.de.study$pseudobulkedGEX$circ, "circ")
p2= plot_corrs(pb.de.study$pseudobulkedGEX$hfpef, "hfpef")
p3= plot_corrs(pb.de.study$pseudobulkedGEX$forte, "forte")

pdf(file = "output/figures/supp/correlation.pb.pdf")
p1
p2
p3
dev.off()


## add HVF per study

studies= unique(seu.int@meta.data$study)
hvf= map(studies, function(x){
  #get hvg:
  seu2= subset(seu.int, study== x)
  seu2= FindVariableFeatures(seu2 )

  #reduce pseudobulk to hvg
  rownames(HVFInfo(seu2))[1:2000]
})

names(hvf)= studies


plot_corrs2= function(pb.df, study.name){

  #pb.df = pb.de$Endothelial$voom
  l1= length(colnames(pb.df))

  df= data.frame(matrix(nrow = l1, ncol= l1),row.names = colnames(pb.df))
  colnames(df)= colnames(pb.df)
  for( x in rownames(df)){
    for( y in rownames(df)){
      print(length(pb.df[hvf[[study.name]],x]))
      df[x,y]= cor.test((pb.df)[hvf[[study.name]],x], pb.df[hvf[[study.name]],y],
                        method = "pearson")$estimate
    }

  }

  target2=
    target %>% filter(study== study.name)
  x= target2$group
  names(x) = target2$orig.ident

  ha = HeatmapAnnotation(group= x, col = list(group= c("ct"= "grey", "hf"= "darkgreen")))
  #ha2 = rowAnnotation(group= x, row = list(group= c("ct"= "grey", "hf"= "darkgreen")))
  Heatmap(df, top_annotation = ha)
}

p1= plot_corrs2(pb.de.study$pseudobulkedGEX$circ, study.name = "circ")
p2= plot_corrs2(pb.de.study$pseudobulkedGEX$hfpef, "hfpef")
p3= plot_corrs2(pb.de.study$pseudobulkedGEX$forte, "forte")


pdf(file = "output/figures/supp/correlation.pb.hvf.pdf")
p1
p2
p3
dev.off()


# run augur -----------------------------------------------------------------------------------
 # devtools::install_github("Bioconductor/MatrixGenerics")
 # devtools::install_github("const-ae/sparseMatrixStats")
 # devtools::install_github("neurorestore/Augur")

 library(Augur)

 augur = calculate_auc(seu, seu@meta.data,
                       cell_type_col = "celltype", label_col = "group")

 p1= plot_lollipop(augur) +
   geom_point(aes(color = cell_type), size = 3) +
    theme(axis.text.y =  element_text(size= 11),
         axis.text.x =  element_text(size= 11),
         legend.position = "none")

 saveRDS(augur, "output/augur.res.rds")
 augur= readRDS("output/augur.res.rds")



 DimPlot(seu, reduction = "umap_original_filt")

 aucs = augur$AUC
 size_sm = 11
 size_lg = 11
 range = range(aucs$auc)
 expand = abs(diff(range)) * 0.1
 p = aucs %>% ggplot(aes(x = reorder(cell_type, auc), y = auc)) +
   geom_hline(aes(yintercept = 0.5), linetype = "dotted",
              size = 0.3) +
   geom_point(size = 2) +
   geom_text(aes(label = format(auc,
                                digits = 3),
                 y = ifelse(auc < 0.5, 0.5, auc)), size = 3,
             nudge_y = expand, hjust = 0.5) +
   geom_segment(aes(xend = cell_type,
                    yend = 0.5)) +
   scale_y_continuous("AUC", limits = c(min(range[1] - expand, 0.5), range[2] + expand * 1.5)) +
   coord_flip() +
   theme_bw() + theme(axis.text.x = element_text(size = size_sm),
                      axis.text.y = element_text(size = size_sm), axis.title.x = element_text(size = size_lg),
                      axis.title.y = element_blank(), panel.grid = element_blank(),
                      strip.text = element_text(size = size_lg), strip.background = element_blank(),
                      axis.line.y = element_blank(), axis.line.x = element_blank(),
                      legend.position = "top", legend.text = element_text(size = size_sm),
                      legend.title = element_text(size = size_sm),
                      legend.key.size = unit(0.6,
                                                                                          "lines"),
                      legend.margin = margin(rep(0, 4)), legend.background = element_blank(),
                      plot.margin = margin(1,1,1.5,1.2, "cm"),
                      plot.title = element_text(size = size_lg, hjust = 0.5))
pdf("output/figures/supp/augur.pdf",
    height= 3,
    width = 5)
p
dev.off()
# run coda ------------------------------------------------------------------------------------

 #install pckgs
 #BiocManager::install(c("clusterProfiler", "DESeq2", "DOSE", "EnhancedVolcano", "enrichplot", "fabia", "GOfuncR", "Rgraphviz"))
 #devtools::install_github('kharchenkolab/cacoa')
 #devtools::install_github("kharchenkolab/sccore", ref="dev")
 #install.packages('coda.base')
 #BiocManager::install("DESeq2", force= T)
 library(cacoa)
#prepare input


 prep_cacoa_input= function(seu,
                            cell.label="celltype",
                            group.label= "group",
                            orig.ident = "orig.ident"){

   x= seu@meta.data%>% as_tibble()%>% distinct(!!as.name(orig.ident), !!as.name(group.label) )
   sample.groups= x[[group.label]]
   names(sample.groups)= x[[orig.ident]]
   sample.groups= as.factor(sample.groups)

   cell.groups= seu@meta.data[[cell.label]]
   names(cell.groups)= rownames(seu@meta.data)
   cell.groups= as.factor(cell.groups)

   sample.per.cell= seu@meta.data[[orig.ident]]
   names(sample.per.cell)= rownames(seu@meta.data)
   sample.per.cell= as.factor(sample.per.cell)

    return(list("sample.groups"= sample.groups,
             "cell.groups"=cell.groups,
             "sample.per.cell"=sample.per.cell)
          )
 }


input= prep_cacoa_input(seu)
input$sample.per.cell
Idents(seu)= seu$orig.ident
cao <- Cacoa$new(seu, ref.level="ct",
                 target.level="hfpef",
                 sample.groups = input$sample.groups,
                 cell.groups = input$cell.groups,
                 sample.per.cell = input$sample.per.cell
                  )

cao$estimateExpressionShiftMagnitudes()
cao$estimateCellLoadings()
cao$plotCellLoadings(show.pvals=T)
cao$estimateDEPerCellType()
cao$estimateDEStabilityPerCellType(top.n.genes = 50)
cao$plotDEStabilityPerCellType()

saveRDS(cao, "output/cao.obj.rds")
cao= readRDS("output/cao.obj.rds")

cao$plotExpressionShiftMagnitudes(notch= F)
cao$plotNumberOfDEGenes()
cao$plotCellGroupSizes()
cao$plotVolcano()

#save plots

pdf("output/figures/coa.res.pdf")
p1
cao$plotCellLoadings(show.pvals=T)
cao$plotExpressionShiftMagnitudes(notch= F)
cao$plotNumberOfDEGenes()
cao$plotCellGroupSizes()
cao$plotVolcano()
dev.off()
