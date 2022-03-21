## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-09-07
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## utility functions



library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
hfpef_cols= c("#EF946C", "#785474")
up_dn_cols= c("#0353A4", "#C41E3D")

# Translation of mouse and human genes --------------------------------------------------------

# Basic function to convert mouse to human gene names
convertMouseToHuman <- function(x, human= NULL, mouse= NULL){
  require("biomaRt")

  if(is.null(human)){
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  }

  if(is.null(mouse)){
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  }

  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

  #now we add both dfs to also report back the NAs
  results= enframe(x, value = "MGI.symbol") %>% left_join(genesV2 ) %>% select(-name)
  # Print the first 6 genes found to the screen
  print(head(results))
  return(list("full.genes"= results,
              "df"= genesV2)
  )
}
# Basic function to convert human to mouse gene names
convertHumanToMouse <- function(x, human= NULL, mouse= NULL){
  require("biomaRt")

  if(is.null(human)){
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  }

  if(is.null(mouse)){
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  }

  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)

  #now we add both dfs to also report back the NAs
  results= enframe(x, value = "HGNC.symbol") %>% left_join(genesV2 )%>% select(-name)
  # Print the first 6 genes found to the screen
  print(head(results))

  return(list("full.genes"= results,
              "df"= genesV2)
  )
}



# get proportions per cluste ------------------------------------------------------------------

calc_props= function(seu.meta, cluster.col="opt_clust_integrated" , group.col= "orig.ident"){

  x= seu.meta%>% group_by(.data[[cluster.col]], .data[[group.col]]) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))

  n= length(unique(x[[cluster.col]]))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  p.OP= ggplot(x, aes(x= .data[[cluster.col]], y= freq, fill = .data[[group.col]] ))+
    geom_col()+
    scale_fill_manual(values = col_vector)+
    theme_minimal()+
    theme(axis.text = element_text(size= 12),
          axis.text.x = element_text(angle=70, hjust= 1))

}

## calculate scaled proportions

calc_scaled_proportion= function(seu.meta,cluster.col= "opt_state", group.col= "group" ){

  #get total cell count to scale on:
  total.counts= seu.meta%>%
    group_by(.data[[group.col]]) %>%
    summarise(n = n())

  # get proportions of study per group
  df.prop= seu.meta%>%
    group_by(.data[[cluster.col]], .data[[group.col]]) %>%
    summarise(n = n()) %>%
    mutate(freq.total = n / sum(n)) %>%
    left_join(total.counts %>%
                dplyr::rename(n.all= n),
              by= group.col) %>%
    mutate(freq.group= n/n.all) #  this is proportion of cells from cluster within a group

  #now we scale freq.group towards total cell counts

  df.prop= df.prop %>%
    mutate(freqs.sum= sum(freq.group)) %>% #just a dummy variable
    ungroup() %>%
    mutate(freq.group.scaled= freq.group/freqs.sum) %>%
    select(-freqs.sum)

  df.FC= df.prop %>% select(.data[[cluster.col]], group, freq.group) %>%
     pivot_wider( names_from = "group", values_from= "freq.group") %>%
    group_by(.data[[cluster.col]])%>%
    mutate(fc= hf/ct)%>%
    select(.data[[cluster.col]], fc)

  df.prop=  df.prop %>% left_join(df.FC, by= .data[[cluster.col]])
  return(df.prop)
}


calc_mean_proportions_per_group= function(seu.meta,
                                          cluster.col= "opt_state",
                                          group.col= "group" ,
                                          sample.col= "orig.ident"){

  #get total cell count to scale on:
  total.counts= seu.meta%>%
    group_by(.data[[cluster.col]],
             .data[[sample.col]]) %>%
    summarise(n = n())

  total.counts.sample= seu.meta%>%
    group_by( .data[[sample.col]]) %>%
    summarise(n.sample = n())

  merged.counts= total.counts %>% left_join(total.counts.sample) %>%
    mutate(props= n/n.sample) %>%
    left_join(seu.meta %>% distinct(.data[[group.col]],
                                    .data[[sample.col]]))


  merged.counts.sum= merged.counts %>%
    group_by( .data[[group.col]],
              .data[[cluster.col]],) %>%
    summarise(mean.percent = mean(props))


  return(list("samplewise"= merged.counts,
              "groupwise"= merged.counts.sum))
}


plot_mean_proportions =function(df,cluster.col = "cellstate",
                                main= "X"){
  plot. = df %>%
    ggplot(.,aes( x= .data[[cluster.col]], y= mean.percent, fill= group))+
    geom_bar(stat="identity", width=.5, position = "dodge") +
    scale_fill_manual(values= hfpef_cols)+
    theme_minimal()+
    theme(axis.text.x = element_text(size= 13, angle = 40, hjust= 1),
          axis.text.y= element_text(size=14))

  if(main != "X") {plot. = plot. +ggtitle(main)}
  return(plot.)
}
# use cluster label to plot basic features of cells:
get_QC_plots= function(seu.meta, cluster.col = "celltype"){

  require(cowplot)
  x= seu.meta

  precent.mt.p= ggplot(x, aes(x= get(cluster.col), y= percent.mt))+
    geom_point()+
    geom_boxplot()
  percent.rb.p= ggplot(x, aes(x= get(cluster.col), y= percent.rb))+
    geom_point()+
    geom_boxplot()
  precent.dis.p= ggplot(x, aes(x= get(cluster.col), y= dissociation_s1))+
    geom_point()+
    geom_boxplot()
  # s.score.p= ggplot(x, aes(x= get(cluster.col), y= S.Score))+
  #   geom_point()+
  #   geom_boxplot()
  # g2m.score.p= ggplot(x, aes(x= get(cluster.col), y= G2M.Score))+
  #   geom_point()+
  #   geom_boxplot()
  count.p= ggplot(x, aes(x= get(cluster.col), y= nCount_RNA))+
    geom_point()+
    geom_boxplot()
  nfeature.p= ggplot(x, aes(x= get(cluster.col), y= nFeature_RNA))+
    geom_point()+
    geom_boxplot()
  doublet.p= ggplot(x, aes(x= get(cluster.col), y= doublet_score))+
    geom_point()+
    geom_boxplot()

  sanity_p1= plot_grid(precent.mt.p,percent.rb.p, precent.dis.p,
                       ncol = 1, align = "v")

  sanity_p2= plot_grid(doublet.p, #g2m.score.p ,
                       #s.score.p,
                       nfeature.p,
                       count.p, ncol = 1, align = "v")

  return(list(sanity_p1, sanity_p2))
}


# module score addition -----------------------------------------------------------------------


#GET NABA genesets:
processNABA = function(filepath = "/home/jan/R-projects/sc-exploration/data/NABAgsets.xls") {
  con = file(filepath, "r")
  naba_gsets = list()
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    split_line = unlist(strsplit(line,split="\t"))
    naba_gsets[[split_line[1]]] = split_line[3:length(split_line)]
  }
  close(con)
  return(naba_gsets)
}

# add module score for listed gene sets to a seurat obj
add_geneset_scores = function(seurat_obj, genesets){

  #DefaultAssay(seurat_obj) = "RNA"
  set.names= names(genesets)

  for(gset in set.names){
    features = list(gset = genesets[[gset]])
    ctrl_genes = sum(genesets[[gset]] %in% rownames(seurat_obj))
    print(gset)
    print(ctrl_genes)

    if(ctrl_genes>500){ctrl_genes= 500}

    if(ctrl_genes>=2){
      seurat_obj = AddModuleScore(object = seurat_obj,
                                  features = features,
                                  name = gset,
                                  ctrl = ctrl_genes,
                                  seed = 77)
    }



  }

  return(seurat_obj)
}





# ORA analysis --------------------------------------------------------------------------------

GSE_analysis = function(geneList,Annotation_DB){
  library(dplyr)
  library(tidyr)
  library(tibble)

  geneList = geneList[geneList %in% unique(unlist(Annotation_DB))]

  ResultsDF = matrix(0,nrow = length(Annotation_DB),ncol = 5)
  rownames(ResultsDF) = names(Annotation_DB)
  colnames(ResultsDF) = c("GenesInPathway","GenesInList","GeneNames","p_value","corr_p_value")

  DB_genecontent = length(unique(unlist(Annotation_DB)))

  GenesDB = DB_genecontent
  SelectedGenes = length(geneList)

  for(gset in rownames(ResultsDF)){
    GP = length(((Annotation_DB[[gset]])))
    GL = length(intersect(Annotation_DB[[gset]],geneList))

    ResultsDF[gset,"GenesInList"] = GL
    ResultsDF[gset,"GenesInPathway"] = GP
    ResultsDF[gset,"GeneNames"] = paste(intersect(Annotation_DB[[gset]],geneList),collapse = ",")
    ResultsDF[gset,"p_value"] = phyper(q=GL - 1, m=GP, n=GenesDB-GP, k=SelectedGenes, lower.tail = FALSE, log.p = FALSE)
  }

  ResultsDF[,"corr_p_value"] = p.adjust(ResultsDF[,"p_value"],method = "BH")
  ResultsDF = data.frame(ResultsDF,stringsAsFactors = F)
  ResultsDF = ResultsDF[order(ResultsDF[,"p_value"]),]

  ResultsDF = ResultsDF %>%
    rownames_to_column("gset") %>%
    mutate_at(c("GenesInPathway","GenesInList",
                "p_value","corr_p_value"),
              as.numeric) %>%
    dplyr::arrange(corr_p_value,GenesInList)

  return(ResultsDF)

}

# function to plot ORA results:
plot_ORA= function(Ora_res, top.marker ){

  for (i in names(Ora_res)){
    Ora_res[[i]] = Ora_res[[i]] %>% mutate(cluster= i)
  }

  Ora_res= as_tibble(do.call(rbind, Ora_res))%>%
    #dplyr::rename(cluster= gset) %>%
    mutate(stars= ifelse(corr_p_value<0.05, "*", ""))



  p.p.ora= ggplot(Ora_res, aes(x=  cluster,
                               y= gset, fill= -log10(corr_p_value)))+
    geom_tile()+
    scale_fill_gradient2(low= "white" , high= "red")+
    geom_text(mapping = aes(label= stars))+
    theme_minimal()+
    theme(axis.text.x= element_text(angle=40, hjust= 1, size= 10),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )+coord_flip()
  p.p.ora
  return(p.p.ora)
}


