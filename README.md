# LOGIC_Platyhelminthes  
scRNA-seq Analysis of Ligand-Gated Ion Channel Conservation (LOGIC) in Schistosoma mansoni

This repository contains the analysis of ligand-gated ion channel (LOGIC) conservation (n = 35) across Platyhelminthes using single-cell RNA sequencing (scRNA-seq) data from Schistosoma mansoni.

**Dataset**
The scRNA-seq data used in this study were generated by George Wendt et al. (DOI: 10.1126/science.abb7709).
Following the publisher’s recommendations, RNA assay data were used for the analysis.

**Repository Contents**
This repository includes the following files:

**celltypes_cluster_res5_map_PMID_32973030.txt** – Mapping of clusters to cell types.
**Gene lists** – Categorized across five groups.
**R script** – Differential gene expression analysis for each cluster against all others.
**Python script** – Code for generating a heatmap summarizing results in muscle and neuronal cell types.

**Usage**
**Differential expression analysis:** Run the provided R scripts to compare expression levels across clusters.
**Heatmap visualization:** Use the Python script to generate a heatmap summarizing key findings.

