# Bioinformatics - RNAseq Analysis

This repository contains the main codes that I used for the analysis of RNAseq data from the laboratory of molecular biology of pathogens (LBMP) of the Federal University of São Paulo (UNIFESP).

## Description

* `up_or_down_genes` this code does a simple logical test using a .tabular file from DESeq2 as input. The parameters used were fdr and log2fc

* `sE_GO` this code uses the simplifyEnrichment library to perform a Gene Ontology Term Enrichment

* `volcano_plot` this code generates a .pdf with the volcano plot for the up, down or not significant. It is also possible to add the label for the most differentially expressed genes (by default it adds to the top 10)

* `deseq_analysis` this code uses the DESeq2 algorithm to perform a differential gene expression analysis between two conditions, in addition to making the necessary normalizations (PCA)
