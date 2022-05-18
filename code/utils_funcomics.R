## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-12
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## functions to run progeny and cytosig for sc analysis

load_cytosig <- function(cytosig_path = "data/input/cytosig/cytosig_signature_centroid.csv",
                         n_genes=500){
  ## read CytoSig
  # models/signatures were obtained from https://github.com/data2intelligence/CytoSig/tree/master/CytoSig
  cyto_signatures <- read.table(cytosig_path,
                                header = TRUE,
                                row.names = 1)

  # Format Cytosig
  cytosig_net <- cyto_signatures %>%
    as.data.frame() %>%
    rownames_to_column("target") %>%
    pivot_longer(-target,
                 names_to = "cytokine",
                 values_to = "weight") %>%
    # keep top 500 genes per cytokine
    group_by(cytokine) %>%
    slice_max(n=n_genes, order_by = abs(weight)) %>%
    select(cytokine, target, weight) %>%
    mutate(mor = if_else(weight>0, 1, -1)) %>%
    mutate(weight = abs(weight)) %>%
    mutate(cytokine = if_else(cytokine=="A",
                              "Activin_A",
                              cytokine)) %>%
    ungroup()

  return(cytosig_net)
}

get_pseudobulk <- function(seurat_object,
                           assay,
                           expr_prop,
                           sum_count_thresh){
  # show celltypes considered
  levels(Idents(seurat_object)) %>%
    map(function(lev) message(lev))

  # Convert Seurat Object to SCE
  sce <- Seurat::as.SingleCellExperiment(seurat_object, assay = assay)
  colLabels(sce) <- SeuratObject::Idents(seurat_object)
  gc()

  # Pseudobulk - sum all counts by celltype (+ gene expression proportions)
  pseudobulk <- scuttle::summarizeAssayByGroup(sce,
                                               ids = colLabels(sce),
                                               assay.type = "counts", # raw counts
                                               statistics = c("sum", "prop"))
  pseudobulk_expr <- pseudobulk@assays@data$sum %>%
    as_tibble(rownames = "gene") %>%
    pivot_longer(-gene, names_to = "celltype", values_to = "sum_count")
  pseudobulk_prop <- pseudobulk@assays@data$prop.detected %>%
    as_tibble(rownames = "gene") %>%
    pivot_longer(-gene, names_to = "celltype", values_to = "prop")

  # Filter according to expression
  pseudo <- pseudobulk_expr %>%
    left_join(pseudobulk_prop, by = c("gene", "celltype")) %>%
    # filter genes not expressed in at least 10% of cells per celltype
    # and only keep those with summed count of at least 10
    filter(prop >= expr_prop) %>%
    filter(sum_count >= sum_count_thresh) %>%
    select(-prop)

  # looks good
  pseudo %>%
    mutate(celltype = as.factor(celltype)) %>%
    ggplot(aes(x=celltype, y = log2(sum_count))) +
    geom_violin(trim=FALSE)

  # Nest by Celltype, format, and normalize
  pseudo %<>%
    group_by(celltype) %>%
    group_nest(.key = "counts") %>%
    # format and log transform
    mutate(logcounts = counts %>%
             map(function(c) c %>%
                   as.data.frame() %>%
                   column_to_rownames("gene") %>%
                   as.matrix() %>%
                   log2()))
}


run_progeny= function(gex.profile, .label, organism= "Mouse", ...){

  prog.res= progeny(expr = gex.profile,
                    organism = organism,
                    scale =T,
                    perm= 1000,
                    z_scores =T)
  rownames(prog.res)= .label
  return(prog.res)
}

wrap_cytosig_progeny_mlm= function(hfpef, ident. = "opt_state", sample_id= "orig.ident"){

  Idents(hfpef)= ident.
  pseudo <- get_pseudobulk(seurat_object =hfpef,
                           assay = "RNA",
                           expr_prop = 0.05,
                           sum_count_thresh = 5)

  message("Pseudobulk Cytokine Enrichment")

  pseudo_cytosig <- pseudo %>%
    mutate(cytosig_res = logcounts %>%
             map(function(logc){
               run_mlm(
                 as.matrix(logc),
                 cytosig_net,
                 .source = "cytokine",
                 .target = "target",
                 .mor = "mor",
                 .likelihood = "weight")%>%
                 # times = 1000,
                 # seed = 1234,
                 # sparse = TRUE,
                 # randomize_type = "cols_independently") %>%
                 # keep only norm weighted mean
                 #filter(statistic == "corr_wmean") %>%
                 # rename
                 select(cytokine=source,
                        NES=score,
                        p_value) %>%
                 # correct p
                 mutate(adj_pvalue = p.adjust(p_value))
             }))

  pseudo_cytosig= pseudo_cytosig  %>% unnest(cytosig_res)


  ### add progeny:
  message("Pseudobulk Pathway Enrichment")

  pseudo_progeny <- pseudo %>%
    mutate(prog_res = logcounts %>%
             map(function(logc){
               run_progeny(
                 logc,
                 " celltype")
             }))

  plot.df= do.call(rbind, pseudo_progeny$prog_res)

  rownames(plot.df)= pseudo$celltype

  funcres= list("cytosig"= pseudo_cytosig,
                "progeny"= plot.df)


}

plot_cytosig_progeny = function(func_res){

  Heatmap_cytosig= func_res$cytosig%>%
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
    ggplot(., aes(x= celltype, y= cytokine, fill = NES, label= star) )+
    geom_tile()+
    geom_text(aes(label= star))+
    scale_fill_gradient2(low= "blue",mid= "white", high= "red")+
    theme_minimal()+
    theme(axis.text.x = element_text(size= 12, angle = 40, hjust= 1))
  library(circlize)
  func_res$progeny= func_res$progeny[c("Col15a1+", "Igfbp3+", "Pi16+", "Cxcl1+", "Cilp+", "Wif1+"),]
  col_fun = colorRamp2(c(min(func_res$progeny ,na.rm = T), 0, max(func_res$progeny, na.rm = T)), c("blue", "white", "red"))

  Heatmap_progeny= Heatmap(t(func_res$progeny),col = col_fun,  cluster_columns = F,
                           name = "progeny_score",
                           column_names_rot = 40,border = T)

  return(list(Heatmap_cytosig, Heatmap_progeny))
}

## function to clean names of GO terms
clean_names= function(df, col= "Term"){
  df[[col]]=  gsub(pattern = "\\(.*",replacement = "",x = df[[col]])
  df
}
