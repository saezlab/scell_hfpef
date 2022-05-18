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
    pb.x.cpm= pb.x.cpm[rownames(pb.x.cpm) %in% rownames(HVFInfo(seu2))[1:5000], ]

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

#cos.p= calc_distances(cells, pb.de, ratios = T)
# if ratios were calculated, use this function to perform one sample wilcox
ps= map(cos.p, function(x){

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
