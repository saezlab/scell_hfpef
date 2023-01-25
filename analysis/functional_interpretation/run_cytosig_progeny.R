## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-05
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## run cytosig on fib data:

library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(decoupleR)
library(progeny)
library(ComplexHeatmap)
library(circlize)
library(stringr)

seu.objs= readRDS("output/seu.objs/cell.state.obj.list.rds")

# hfpef= readRDS(file = "output/cell_specific_obj/fibroblasts.rds")
#
# endos= readRDS(file = "output/cell_specific_obj/Endothelial.rds")




# utils:  -------------------------------------------------------------------------------------
source("code/utils_funcomics.R")


# load cytosig --------------------------------------------------------------------------------

cytosig_net= load_cytosig("data/prior_knowledge/cytosig_signature_centroid.csv",n_genes = 100)

ggplot(cytosig_net, aes(x= cytokine, y= weight))+
  geom_boxplot()
#capitalize network for quick check in mouse :

cytosig_net$target = str_to_title(cytosig_net$target)
cytosig_net %>% filter(cytokine== "IL3") %>% print(n=100)

# run on cell type/states----------------------------------------------------------------------------
seu= seu.objs$fibroblasts

wrap_funcs= function(seu, name, grouping2= "cellstate"){

  group_wise= wrap_cytosig_progeny_mlm(seu, "group")
  p.group= plot_cytosig_progeny(func_res = group_wise)

  sample_wise= wrap_cytosig_progeny_mlm(seu, "orig.ident")
  p.group= plot_cytosig_progeny(func_res = group_wise)

  state_wise=  wrap_cytosig_progeny_mlm(seu, "cellstate")
  p.state= plot_cytosig_progeny(state_wise)

  seu@meta.data= seu@meta.data %>% mutate(state_id= paste0(cellstate,".",group))
  state_group_wise= wrap_cytosig_progeny_mlm(seu, "state_id")
  p.state.group.wise= plot_cytosig_progeny(state_group_wise)
  p.state.group.wise+ theme(axis.text.x = element_text(angle= 40))

  #   pdf(paste0("output/figures/funcomics/cyto",name, ".pdf"),
  #     width= 6,
  #     height= 12)
  #
  # plot(p.group[[1]])
  # plot(p.state[[1]])
  # plot(p.state.group.wise[[1]])
  # dev.off()
  # #dev.off()
  #
  # pdf(paste0("output/figures/funcomics/progeny",name, ".pdf"),
  #     width= 6,
  #     height= 6)
  #
  # print(p.group[[2]])
  # print(p.state[[2]])
  # print(p.state.group.wise[[2]])
  # dev.off()
  # #dev.off()
 return(list(p.group, p.state,p.state.group.wise))
}


# run on all celltypes -----------------------------------------------------------------------------------
state_wise=  wrap_cytosig_progeny_mlm(hfpef = seu, ident. = "cellstate")
p.state= plot_cytosig_progeny(state_wise)
saveRDS(p.state[[2]], file= "output/figures/main/Fig2/progeny.rds")

x= wrap_funcs(seu =seu.objs$fibroblasts, name= "fibroblasts")
pdf("output/figures/funcomics/progeny_fibroblast_group.pdf",
        width= 3,
        height= 4)
  x[[1]][[2]]
    dev.off()

pdf("output/figures/funcomics/progeny_fibroblast_state.pdf",
        width= 4,
        height= 4)
    x[[2]][[2]]
dev.off()

pdf("output/figures/funcomics/progeny_fibroblast_groupstate.pdf",
        width= 5,
        height= 4)
    x[[3]][[2]]
dev.off()

y= wrap_funcs(seu.objs$Endothelial, "Endothelial")
m= wrap_funcs(seu.objs$macrophages, "macs")


group_wise$cytosig$celltype
# run on whole atlas: -------------------------------------------------------------------------

seu = readRDS("output/seu.objs/integrated_cellstate_noECs.rds")
whole.cell=  wrap_cytosig_progeny_mlm(seu, "celltype")
p.state= plot_cytosig_progeny(whole.cell)

seu@meta.data= seu@meta.data %>% mutate(type_id= paste0(celltype,".",group ))
state_group_wise= wrap_cytosig_progeny_mlm(seu, "type_id")
p.state.group.wise= plot_cytosig_progeny(state_group_wise)


pdf(paste0("output/figures/funcomics/cyto_celltypes.pdf"),
    width= 6,
    height= 15)
plot(p.state[[1]])
plot(p.state.group.wise[[1]])
dev.off()

pdf(paste0("output/figures/funcomics/progeny_celltypes.pdf"),
    width= 6,
    height= 6)

print(p.state[[2]])
print(p.state.group.wise[[2]])

dev.off()



# integrated_fib ------------------------------------------------------------------------------

seu_int= readRDS("output/seu.objs/study_integrations/harmony_fib_filt_2.rds")
unique(seu_int$opt_clust_integrated)
Seurat::DimPlot(seu_int)

whole.cell=  wrap_cytosig_progeny_mlm(seu_int, "opt_clust_integrated")
p.state= plot_cytosig_progeny(whole.cell)



#func_res$progeny= func_res$progeny[c("Col15a1+", "Igfbp3+", "Pi16+", "Cxcl1+", "Cilp+", "Wif1+"),]
col_fun = colorRamp2(c(min(whole.cell$progeny ,na.rm = T), 0, max(whole.cell$progeny, na.rm = T)), c("blue", "white", "red"))
Heatmap_progeny= Heatmap(t(whole.cell$progeny),
                        col = col_fun,
                        cluster_columns = F,
                        name = "PROGENy\nscore\n",
                        column_names_rot = 40,
                        border = T,
                        rect_gp = gpar(col = "darkgrey", lty = 1, size= 0.1),
                        row_names_gp = gpar(fontsize = 10),
                        column_names_gp = gpar(fontsize = 10),
                        row_names_side = "left",
                        show_row_dend = F)
Heatmap_progeny

pdf(paste0("output/figures/main/Fig2/progeny_celltypes.pdf"),
    width= 3.7,
    height=3)
print(Heatmap_progeny)
dev.off()

Heatmap_cytosig= whole.cell$cytosig%>%
  mutate(star= ifelse(adj_pvalue < 0.01,
                      "***",
                      ifelse(adj_pvalue< 0.05,
                             "** ",
                             ifelse(adj_pvalue < 0.1,
                                    "* ",
                                    " ")
                      )
  )
  ) %>%
  filter(cytokine != "IL3") %>%
  ggplot(., aes(x= celltype, y= cytokine, fill = NES, label= star) )+
  geom_tile()+
  geom_text(aes(label= star))+
  scale_fill_gradient2(low= "blue",mid= "white", high= "red")+
  theme_minimal()+
  theme(axis.text.x = element_text(size= 12, angle = 40, hjust= 1))
Heatmap_cytosig




