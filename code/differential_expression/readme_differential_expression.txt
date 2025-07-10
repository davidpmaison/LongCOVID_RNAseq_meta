Differential Expression Analysis (DESeq2)
Author: David Maison
==========================================

This directory contains scripts for performing differential expression analysis using DESeq2, annotating results with gene symbols, and generating volcano plots.

Directory Contents:
--------------------
1. DESeq2_differential_expression_template.R
   - Main pipeline for differential expression analysis.
   - Reads a count matrix and metadata CSV.
   - Runs DESeq2 to identify differentially expressed genes.
   - Exports results as a CSV.

2. annotate_DEG_results_with_biomaRt.R
   - Annotates DESeq2 results by mapping ENSEMBL IDs to gene symbols and gene biotypes using biomaRt.
   - Outputs an annotated results file.

3. volcano_plot_template.R
   - Generates a volcano plot of DESeq2 results using EnhancedVolcano.
   - Saves the plot as a PDF.

Required Input Files:
---------------------
1. Counts matrix (generated from STAR quantifications) with ENSEMBL IDs as rows and sample IDs as columns.
2. Metadata file (`metadata.csv`) containing at minimum:
   - sample_id: sample name matching count matrix column names.
   - condition: experimental group (e.g., "CR", "LC").
   - study_id: optional covariate if combining datasets.

Workflow:
---------
1. Place all STAR ReadsPerGene.out.tab files into a directory.
2. Generate counts matrix and metadata.csv (see DESeq2_differential_expression_template.R for details).
3. Run DESeq2_differential_expression_template.R to obtain differential expression results.
4. Annotate results with annotate_DEG_results_with_biomaRt.R.
5. Generate volcano plots with volcano_plot_template.R.

Notes:
------
- Adjust file paths in the scripts as needed.
- This pipeline assumes paired group comparisons (e.g., LC vs CR).
- Scripts are designed for bulk RNA-seq datasets with gene-level counts.
