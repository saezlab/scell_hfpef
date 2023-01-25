## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-02-09
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Assess overlap between cell state marker and DEGs

library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(ggvenn)
library(circlize)

gene_signatures= readRDS( "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")

logFCs=readRDS(file = "output/fib_integration/marker_list/DEG_per_study_LOGFC.rds")

names(logFCs)=c("HFpEF", "MI", "AngII")

int_state_marker= readRDS("output/fib_integration/marker_list/integrated_marker.rds")

# load sc data
seu.objs=  readRDS("output/seu.objs/cell.state.obj.list.rds")
seu= seu.objs$fibroblasts
rm(seu.objs)
cut.off.state= 100
Idents(seu)= "cellstate"
marker.fibs= FindAllMarkers(seu, assay= "RNA")

loc_state_marker= marker.fibs%>%
  group_by(cluster) %>%
  filter(p_val_adj <0.05)%>%
  slice_max(order_by = (avg_log2FC), n = cut.off.state)

#function to calculate overlap between state marker and gene signatures.
#' @return percentage of state marker that are DEGs
calculate_per_cluster= function(cell_state_marker,gene_signatures, top_marks = 100){
  #x= gene_signatures$total[1]
  df=sapply(names(gene_signatures), function(x){
    #print(x)
    intersects= lapply(unique(cell_state_marker$cluster), function(y){
      #print(y)
      marks= cell_state_marker %>% filter(cluster== y,
                                   p_val_adj<0.05,
                                   avg_log2FC>0)%>%
        slice_max(., order_by = (avg_log2FC),n= top_marks)%>%pull(gene)
      #length(intersect(gene_signatures[[x]], marks))/top_marks
      length(intersect(gene_signatures[[x]], marks))/length(gene_signatures[[x]])
    })
    names(intersects)=unique(cell_state_marker$cluster)
   unlist(intersects)
  })
  return(df*100)
}

gene_signatures$total= gene_signatures$total[names(gene_signatures$total) != "MI"]
df= calculate_per_cluster(cell_state_marker = int_state_marker,
                      gene_signatures = gene_signatures$total,
                      top_marks = 100)

df2= calculate_per_cluster(cell_state_marker = loc_state_marker,
                          gene_signatures = gene_signatures$total,
                          top_marks = 100)

col_fun = colorRamp2(c( 0, max(df)), c( "white", "red"))

p.int= Heatmap(df, col = col_fun, name = "% overlap",
               cluster_columns = F,
               cluster_rows = F,
               cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.1f", df[i, j]), x, y, gp = gpar(fontsize = 10))
})


col_fun = colorRamp2(c( 0, max(df2)), c( "white", "red"))
#col_fun = colorRamp2(c( 0,2), c( "white", "red"))
df2= round(df2, 0)
p.loc= Heatmap(round(df2, 0), col = col_fun, name = "% overlap", cluster_columns = F,
               cluster_rows = F,
               cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.1f", df2[i, j]), x, y, gp = gpar(fontsize = 10))
})

p.loc
p.int

pdf("output/figures/supp/deg.hmaps/deg_overlap_state_marker.pdf",
    width= 2.9,
    height= 2.5)
p.loc
dev.off()

pdf("output/figures/supp/deg.hmaps/deg_overlap_int_state_marker.pdf",
    width= 2.8,
    height= 2.5)
p.int

dev.off()


# props for DEG -------------------------------------------------------------------------------
top_marks= 50
  calculate_per_sig= function(cell_state_marker,gene_signatures, top_marks = 100){

  df=sapply(names(gene_signatures), function(x){
    #print(x)
    marks= cell_state_marker %>%
      group_by(cluster)%>%
      filter(p_val_adj<0.05,
             avg_log2FC>0)%>%
      slice_max(., order_by = (avg_log2FC),n= top_marks)%>%pull(gene)

    length(intersect(gene_signatures[[x]], marks))/length(gene_signatures[[x]])

  })
  return(df*100)

}

df= calculate_per_sig(cell_state_marker= int_state_marker,gene_signatures =  gene_signatures$total, top_marks = 100)





# calculate ratio of FCs per gene -------------------------------------------------------------


studieS= map(names(logFCs), function(x){
  df1= logFCs[[x]]
    df1= df1[gene_signatures$total[[x]],]%>% rownames_to_column("gene")%>%
    rename(logFC_DEG =avg_log2FC)


  comb= marker.fibs%>%
    as_tibble() %>%
    group_by(gene)%>%
    summarise(max_fc= max(avg_log2FC))%>%
    #mutate(avg_log2FC = ifelse(avg_log2FC<0 , min(abs(avg_log2FC)), avg_log2FC))%>%
    #rownames_to_column("gene")%>%
    left_join(df1, by="gene")%>%
    mutate(r.fc=  logFC_DEG/ max_fc)%>% arrange(desc(r.fc))%>%
    drop_na()  %>% mutate(study= x)

})  %>% do.call(rbind, .)


studieS %>% ggplot(., aes(x= study, y= r.fc))+
  geom_boxplot()+
  geom_jitter()



df1= logFCs$HFpEF[ gene_signatures$total$HFpEF,]%>% rownames_to_column("gene")%>%
  rename(logFC_DEG =avg_log2FC)

studieS %>% filter(study== "HFpEF") %>% arrange(desc(logFC_DEG))
comb= marker.fibs%>% as_tibble() %>%
  group_by(gene)%>% summarise(max_fc= max(avg_log2FC))%>%
  #rownames_to_column("gene")%>%
  left_join(df1, by="gene")%>%
  mutate(r.fc=  logFC_DEG/ max_fc)%>% arrange(desc(r.fc))%>%
  drop_na() %>%
  print(n=100)

