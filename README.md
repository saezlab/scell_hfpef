# sc_hfpef

**Collaborators**:
Laura Wienecke and Florian Leuschner. 

**Background**: Heart Failure with preserved ejection fraction (HFpEF) makes up about 50 % of the HF population and is a systemic syndrome characterized by accumulated risk factors and comorbidities and a noncompliant and stiff heart that is exposed to lower wall stress. HFpEF often is described as disease that is not "heart centric", i.e. the heart  falls victim to a pro inflammatory and/or dysbalanced metabolic state that arises from co-morbidities like diabetes or chronic kidney disease.
Multiple disease models have been proposed to study HFpEF (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5764178/#!po=25.6098), and only recently a promising model was established: https://www.nature.com/articles/s41586-019-1100-z. 
This model has been investigated for its symptomatic phenotype that resembles HFpEF characteristic but has not been described on the single cell level yet. 

**Methods**: whole single cell RNA sequencing, 10x genomics, done in sc facility Mannheim. 

**Design:** Overall: 2x Control mice vs 2x hfpef mice, time points 5-6 weeks after treatment start.

Project Overview: https://docs.google.com/document/d/1oSFvogHiYlhuFfhMWhfwn2UsAHQ2MFsN0yaTN8pQPYk/edit?usp=sharing


**Analysis Workflow**

Processing, QC and integration
1) Run sample wise processing.R [sample wise QC]
2) harmony_integration.R [first sample integration]
3) cluster_label.R [Cluster marker identification and manual labeling]
4) filter_cluster.R [based on marker genes, some clusters are removed and samples are reintegrated]

Cellstate definitions and Differntial gene expression
1) create_celltype_atlas.R 
2) marker_gene_cell_states.R [cell states are defined ]
3) pseudobulk_DE_celltypes [HfpEF DEA per celltype]
4) differntial_cell_proportions.R [calculate cell proportion changes with label permutation]

Funcomics
1) run_cytosig_progeny.R [pathway analysis and cytokine footprinting]
2) run_dorothea.R [TF activity estimation]
3) run_liana.R [Ligand Receptor analysis]
4) run_nichenet.R [Ligand - footprint analysis]
5) GO_enrich [Enrichment of Gene Ontology terms]

Fibroblast Atlas Integration
1) Process study1
2) Process study2
1) integration_fibroblasts.R + integration_fibroblasts_filtered.R [Harmony based integration of three studies with cells filtered]
3) composition_comparison.R [Compare cell HFpEF cell states mapping]
4) DEA_fibs_all_studies.R [Fibroblas Marker gene detection in each study]
5) map_study_DEG_to_integrated_cluster.R [fibroblast marker are mapped back to each study]
6) 





A [workflowr][] project.

[workflowr]: https://github.com/jdblischak/workflowr
