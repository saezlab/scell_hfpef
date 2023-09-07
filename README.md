# Comparison of HFpEF fibrotic signatures

**Collaborators**:
Laura Wienecke and Florian Leuschner. 

**Background**:
Inflammation, fibrosis and metabolic stress critically promote heart failure with preserved ejection fraction (HFpEF). Exposure to high-fat diet and nitric oxide synthase inhibitor N[w]-nitro-l-arginine methyl ester (L-NAME) recapitulate features of HFpEF in mice. To identify disease specific traits during adverse remodeling, we performed single-cell RNAseq (scRNAseq) of interstitial cells in murine HFpEF. Diastolic dysfunction and fibrosis were accompanied by an activation of cardiac fibroblast and macrophage subsets. Comparison with scRNAseq datasets of murine heart failure with reduced ejection fraction (HFrEF) identified a specific HFpEF fibroblast disease signature, characterized by e.g. overexpression of basement membrane genes. While myo- and matrifibroblast activation is known to be crucial for cardiac fibrosis in HFrEF, we found this cell-state switch to be less important in HFpEF. Disease specific fibroblast signatures were corroborated in human myocardial bulk transcriptomes and main markers confirmed at protein level in murine samples. Lipoproteinlipase Angiopoietin-like 4 was identified to indicate HFpEF mediated fibroblast activation and may serve as a potential biomarker for disease progression. Taken together, our results provide a characterization of the interstitial cellular landscape of murine HFpEF, including specific characteristics of fibroblast activation.

**Methods**: Whole single cell RNA sequencing with 10x genomics and Illumina 

**Analysis Workflow**

Processing, QC and integration
1) sample wise QC [run_samplewise_processing.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/sample_integration/run_sample_wise_preprocessing.R) 
2) first sample integration [harmony_integration_across_sample.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/sample_integration/harmony_integration_across_sample.R)
3) Cluster marker identification and manual labeling [annotate_celltypes.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/sample_integration/annotate_celltypes.R) 
4) filtering of clusters and reintegration [round1](https://github.com/saezlab/scell_hfpef/blob/main/analysis/sample_integration/filter_clusters_round1.R) [round2](https://github.com/saezlab/scell_hfpef/blob/main/analysis/sample_integration/filter_clusters_round2.R)
5) calculate distance between pseudobulk profiles [sample_distance.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/sample_integration/sample_distance.R)
6) calculate cell proportion changes with label permutation [differntial_cell_proportions.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/functional_interpretation/differntial_cell_proportions.R)

Cellstate definitions and Differntial gene expression
1) [harmony_integration_per_celltype.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/sample_integration/harmony_integration_per_celltype.R)
2) [annotate_cellstates.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/sample_integration/annotate_cellstates.R)
3) [fibroblast_state_interpretation.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/functional_interpretation/fibroblast_state_interpretation.R)
4) HfpEF DEA per celltype with downsampling [DEA_per_celltype.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/differential_expression_analysis/DEA_per_celltype.R)

Functional analysis
1) Enrichment of Gene Ontology terms  [GO_enrich.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/functional_interpretation/GO_enrich.R)
2) pathway analysis and cytokine footprinting [run_cytosig_progeny.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/functional_interpretation/run_cytosig_progeny.R)
3) TF activity estimation [run_dorothea.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/functional_interpretation/run_dorothea.R)
4) Ligand Receptor analysis [run_liana.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/functional_interpretation/run_liana.R) 
5) Ligand-network analysis [run_nichenet.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/functional_interpretation/run_nichenet.R)

Fibroblast integration with different studies
1) [Process_study_MI.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/study_integration/process_MI.R)
2) [Process_study_AngII.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/study_integration/process_AngII.R)
3)  Harmony based integration of three studies with cells filtered [integration_fibroblasts.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/study_integration/integration_fibroblasts.R) + [integration_fibroblasts_filtered.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/study_integration/integration_fibroblasts_filtered.R)
4) [fibroblast_state_interpretation_integrated.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/functional_interpretation/fibroblast_state_interpretation_integrated.R)
5) [composition_comparison.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/functional_interpretation/differntial_cell_proportions.R) 

Fibroblast signature comparisons
1) Fibroblas Marker gene detection in each study [DEA_fibs_all_studies.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/differential_expression_analysis/DEA_fibs_all_studies.R)
2 GSEA, progeny, dorothea on signatures [fib_signature_interpretation.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/study_comparison/fib_signature_interpretation.R) 
3) compare overlap and correlation between signatures [compare_study_degs.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/study_comparison/compare_stuy_degs.R)
4) Fibroblast marker are mapped back to each study [map_study_DEG_to_integrated_cluster.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/study_comparison/map_study_DEG_to_integrated_cellstates.R)
5) relate statemarker and signatures [compare_statemarker_and_deg.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/study_comparison/compare_statemarker_and_deg.R)

Human bulk comparison
1) Process Human HFpEF [bulk_human_process.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/bulk/bulk_human_process.R)
2) [compare bulk signatures](https://github.com/saezlab/scell_hfpef/blob/main/analysis/study_comparison/bulk_validation_external.R)




A [workflowr][] project.

[workflowr]: https://github.com/jdblischak/workflowr
