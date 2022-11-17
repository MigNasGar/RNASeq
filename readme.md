# Bioinformatics - RNAseq Analysis

This repository contains the main codes that I used for the analysis of RNAseq data from the laboratory of molecular biology of pathogens (LBMP) of the Federal University of São Paulo (UNIFESP).

## Description

* `up_or_down_genes` this code does a simple logical test using a .tabular file from DESeq2 as input. The parameters used were fdr and log2fc

* `sE_GO` this code uses the simplifyEnrichment library to perform a Gene Ontology Term Enrichment

* `volcano_plot` this code generates a .pdf with the volcano plot for the up, down or not significant. It is also possible to add the label for the most differentially expressed genes (by default it adds to the top 10)

* `deseq_analysis` this code uses the DESeq2 algorithm to perform a differential gene expression analysis between two conditions, in addition to making the necessary normalizations (PCA)

* `RPKM_normalization` with this [spreadsheet](https://docs.google.com/spreadsheets/d/1Hcq98c6PZ5QeFtVHtRMj6VApoqYxXAnoNvw4ijTEke8/edit?usp=sharing) is possible to make a Log2RPKM normalization using reads counts from RNAseq runs. To use the worksheet, add to it the number of counts for each desired gene, the total sum of reads per experiment and the size of each gene (in kb). The worksheet is configured for 9 genes but can be expanded to more or less values. ⚠️ Make a copy to your google drive to use it

