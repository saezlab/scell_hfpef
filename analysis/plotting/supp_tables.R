## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-05-09
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## save suppp tables

# supp table 1- HFS marker
seu.objs=  readRDS("output/seu.objs/cell.state.obj.list.rds")

seu= seu.objs$fibroblasts
rm(seu.objs)
Idents(seu)= "cellstate"
df= FindAllMarkers(seu)
df = df %>% filter(p_val_adj<0.01, avg_log2FC>0.01)
df2= split(df, f= df$cluster)
WriteXLS::WriteXLS(df2, ExcelFileName = "output/supplement_table1_HFS_marker.xlsx",
                   SheetNames = names(df2))

# supp table 2- IFS marker --------------------------------------------------------------------------------

IFS_marker =readRDS("output/fib_integration/marker_list/integrated_marker.rds")
IFS_marker= IFS_marker %>% filter(p_val_adj<0.01, avg_log2FC>0.01)
df= split(IFS_marker, f= IFS_marker$cluster)
WriteXLS::WriteXLS(df, ExcelFileName = "output/supplement_table2_IFS_marker.xlsx",
                   SheetNames = names(df))

# supp table3 - disease signatures ---------------------------------------------------------------------

gene_signatures= readRDS( "output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")
df= readRDS( file = "output/fib_integration/marker_list/DEG_per_study_LOGFC.rds")
names(df)= names(gene_signatures$total)

sigs_= map(names(df), function(x){
  y= df[[x]][gene_signatures$total[[x]],]
  rownames_to_column(y, "gene")%>% arrange(desc(abs(avg_log2FC)))
})

WriteXLS::WriteXLS(sigs_, ExcelFileName = "output/supplement_table3_fibrosignatures.xlsx",
                   SheetNames = names(df))



