CASP9 RCC Analysis Pipelines

This repository contains scripts used for single-cell, spatial, SMR, and bulk transcriptome analyses in the CASP9 clear cell renal cell carcinoma (ccRCC) study.

Files

bulk_pipeline.R
Bulk RNA-seq pipeline: differential expression, survival analysis (univariate/multivariate Cox), LASSO modeling, risk score validation, immune infiltration, pathway, and checkpoint analysis.

geneset_pipeline.R
Gene set–based single-cell pipeline: scoring and grouping by gene sets, CellChat communication, pseudotime trajectory analysis, spatial deconvolution, and multi-scale modeling.

single_gene_pipeline.R
Single-gene–based single-cell pipeline (e.g., CASP9/8/10): expression grouping, cell–cell communication analysis, pseudotime, and spatial validation.

prep_data.R
Data preparation: download and process GEO/TCGA datasets, normalization, and generation of matrices or RDS objects for downstream analysis.

smr_prep_and_plots.R
SMR analysis: prepare GWAS/eQTL input, process SMR output, filter significant genes, and generate locus/regional plots.

commot.py
Python script for spatial transcriptomics communication and co-localization analysis (based on commot/scanpy).

deepsurv.py
Python script for deep survival modeling (DeepSurv) using PyTorch/PyCox.
