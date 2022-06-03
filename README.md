# Comparison of HFpEF specific fibrotic signatures

**Collaborators**:
Laura Wienecke and Florian Leuschner. 

**Background**: Heart Failure with preserved ejection fraction (HFpEF) makes up about 50 % of the HF population and is a systemic syndrome characterized by accumulated risk factors and comorbidities and a noncompliant and stiff heart that is exposed to lower wall stress. HFpEF often is described as disease that is not "heart centric", i.e. the heart  falls victim to a pro inflammatory and/or dysbalanced metabolic state that arises from co-morbidities like diabetes or chronic kidney disease.
We used a two hit hfpef mouse model (https://www.nature.com/articles/s41586-019-1100-z) to connect phenotypic changes to cellular changes on the single cell level. 

**Methods**: whole single cell RNA sequencing, 10x genomics, done in sc facility Mannheim. 

**Design:** Overall: 2x Control mice vs 2x hfpef mice, time points 5-6 weeks after treatment start.

**Analysis Workflow**

Processing, QC and integration
1) run_samplewise_processing.R [sample wise QC]
2) harmony_integration_across_sample.R [first sample integration]
3) annotate_celltypes.R [Cluster marker identification and manual labeling]
4) reintegrate_after_filter.R [based on marker genes, some clusters are removed and samples are reintegrated]
5) sample_distance.R [calculate distance between pseudobulk profiles]
6) differntial_cell_proportions.R [calculate cell proportion changes with label permutation]

Cellstate definitions and Differntial gene expression
1) harmony_integration_per_celltype.R 
2) annotate_cellstates.R [cell states are defined ]
3) fibroblast_state_interpretation.R
4) DEA_per_celltype.R [HfpEF DEA per celltype with downsampling]

Functional analysis
1) GO_enrich [Enrichment of Gene Ontology terms]
2) run_cytosig_progeny.R [pathway analysis and cytokine footprinting]
3) run_dorothea.R [TF activity estimation]
4) run_liana.R [Ligand Receptor analysis]
5) run_nichenet.R [Ligand-network analysis]
32  
Fibroblast Integration Meta 
1) Process_study_MI.R
2) Process_study_Ang2.R
3) integration_fibroblasts.R + integration_fibroblasts_filtered.R [Harmony based integration of three studies with cells filtered]
3) fibroblast_state_interpretation_integrated.R
4) composition_comparison.R [Compare cell HFpEF cell states mapping]

Fibroblast disease model signatures
1) DEA_fibs_all_studies.R [Fibroblas Marker gene detection in each study]
2) fib_signature_interpretation.R [GSEA, progeny, dorothea on signatures]
3) compare_study_degs.R [compare overlap and correlation between signatures]
4) map_study_DEG_to_integrated_cluster.R [fibroblast marker are mapped back to each study]
5) compare_statemarker_and_deg.R [relate statemarker and signatures]






A [workflowr][] project.

[workflowr]: https://github.com/jdblischak/workflowr
