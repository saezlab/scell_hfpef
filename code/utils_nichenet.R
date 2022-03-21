## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-15
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## NICHNET utils


# function to prepare nichenet input
#' @return, list of (background genes, geneset oi and potential ligands)

prep_nichenet_input= function(receiver, sender, geneset_receiver= NULL,geneset_sender= NULL, ...){

  ##1) define expressed genes in sender, reciever
  # receiver:
  expressed_genes_receiver = get_expressed_genes(receiver, seu, pct = 0.10, assay_oi = "RNA")
  background_expressed_genes = expressed_genes_receiver %>%
    .[. %in% rownames(ligand_target_matrix)]

  # sender:
  list_expressed_genes_sender = sender  %>% unique() %>% lapply(get_expressed_genes, seu, 0.1, assay_oi= "RNA") # lapply to get the expressed genes of every sender cell type separately here

  expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()


  ##2) define the genset_oi

  #if a target geneset is provided (here upregulted marker after MI), use only those:
  # otherwise perform DEA on condition for receiver cells
  if(!is.null(geneset_receiver)){
    geneset_oi= intersect(geneset_receiver, expressed_genes_receiver)
  }else{  # if no geneset is provided, use upregulated markers in receiver
    Idents(seu)= seu$celltype.stim
    geneset_oi = rownames(get_up_marker(receiver, seu)%>% filter(p_val_adj<0.05))
    Idents(seu)= seu$celltype
  }

  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

  ## 3  define a set of potential ligands (expressed by sender, putative receptor expressed by receiver)
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()

  expressed_ligands = intersect(ligands,expressed_genes_sender)

  # If a sender geneset is provided, we will reduce potential ligands to those genes
  if(!is.null(geneset_sender)){
    expressed_ligands= intersect(geneset_sender,expressed_ligands)
  }

  expressed_receptors = intersect(receptors,expressed_genes_receiver)

  potential_ligands = lr_network %>%
    filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
    pull(from) %>%
    unique()


  return(list(geneset_oi= geneset_oi,
              background_expressed_genes= background_expressed_genes,
              potential_ligands= potential_ligands
  )
  )

}



run_nichenet= function(sender, receiver, geneset_receiver, geneset_sender){

  #PREP

  input= prep_nichenet_input(receiver = receiver,
                             sender = sender,
                             geneset_receiver =geneset_receiver,
                             geneset_sender = geneset_sender)

  message("nichenet_input_prepared")
  # MAIN
  # predict ligand activity
  ligand_activities = predict_ligand_activities(geneset = input$geneset_oi,
                                                background_expressed_genes = input$background_expressed_genes,
                                                ligand_target_matrix = ligand_target_matrix,
                                                potential_ligands = input$potential_ligands)

  ligand_activities = ligand_activities %>%
    arrange(-pearson) %>%
    mutate(rank = rank(desc(pearson)))

  ligand_activities

  ## plot - histogram of pearson corr for all ligands
  p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) +
    geom_histogram(color="black", fill="darkorange")  +
    # geom_density(alpha=.1, fill="orange") +
    geom_vline(aes(xintercept=min(ligand_activities %>% top_n(10, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) +
    labs(x="ligand activity (PCC)", y = "# ligands") +
    theme_classic()

  p_hist_lig_activity

  # this plot can help to select the best top ligands. Here we only have a few ligands we tested,
  # because of the small overlap of upregulated genes in T-cells AND ligand sin the matrix
  # therefore I will use all lignds for downstream analysis, as all ligands display corr>0.1

  best_upstream_ligands = ligand_activities %>%
    #top_n(10, pearson) %>%
    arrange(-pearson) %>%
    pull(test_ligand) %>%
    unique()

  # Data transformation for heatmap plotting
  active_ligand_target_links_df = best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links, geneset = input$geneset_oi, ligand_target_matrix = ligand_target_matrix) %>% bind_rows() %>% drop_na()

  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.1)

  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

  #  heatmap plotting
  p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))

  ##
 # Idents(seu) = seu$celltype.stim
  #p.dot.ligands= DotPlot(subset(seu, idents = c(paste0(sender,"_none"), paste0(sender,"_myocardial.infarction"))), features = best_upstream_ligands %>% rev(), cols = "RdYlBu") +coord_flip()

  Idents(seu)= seu$celltype


  return(list(plots= list(p_ligand_target_network,
                          plot_grid(p_hist_lig_activity, #p.dot.ligands ,
                                    nrow= 2))))




}

#' Basic function to convert human to mouse genesymbols (temporary solution)
#' @param op_resource omnipath_resource as obtained via `liana::select_resource`
#'
#' @details adapted from https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
convert_to_murine <- function(op_resource){

  # query biomaRt databases
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  # obtain tibble with human and murine genesymbol
  symbols_tibble <- getLDS(attributes = c("hgnc_symbol"),
                           filters = "hgnc_symbol",
                           values = union(op_resource$source_genesymbol,
                                          op_resource$target_genesymbol),
                           mart = human,
                           martL = mouse,
                           attributesL = c("mgi_symbol")) %>%
    dplyr::rename(human_symbol = HGNC.symbol,
                  murine_symbol = MGI.symbol) %>%
    as_tibble()

  # intentionally we introduce duplicates, if needed
  # these should be resolved when LIANA is called
  # as inappropriately matched genes will not be assigned any values
  op_resource %>%
    left_join(symbols_tibble, by=c("target_genesymbol"="human_symbol")) %>%
    mutate(target_genesymbol = murine_symbol, .keep = "unused") %>%
    left_join(symbols_tibble, by=c("source_genesymbol"="human_symbol")) %>%
    mutate(source_genesymbol = murine_symbol, .keep = "unused") %>%
    filter(!is.na(target_genesymbol) | !is.na(source_genesymbol)) %>%
    filter(!is.na(target_genesymbol)) %>%
    filter(!is.na(source_genesymbol))
}

