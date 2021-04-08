# Mitogenome evolution and implications on the phylogeny of the Lacerta viridis complex

Welcome to the **Lacerta viridis mitogenome** repository!

This repository is a collection of several scripts and mini-tutorials guiding you through the methods and analyses which were performed in the paper by [Jauss et al., 2021](https://doi.org/10.1080/14772000.2021.1912205).

The raw data can be downloaded [here](), plots and figures were generated with the final trees and miscellaneous files accessible in the folder [00_Data](00_Data/). 

## Table of Content
*Relevant scripts and tutorials can be found in the corresponding markdown files, which are* **linked and highlighted** *for each category.*

### 00 Data
This folder contains the trees and miscellaneous data used for the subsequent analyses and visualisations. Intermediate files are not provided here, but you can generate them yourself by following the next steps.

### 01 Metabarcoding Pipeline
**[In this pipeline](01_Pipeline/Readme.md)**, you find the neccessary scripts to generate the phylogenetic trees from raw .fastq files.

### 02 Constrained Topology Analysis
In this folder you find all input and output files for the constrained topology analysis performed with `iqtree`. The constraints were performed with RAxML based on the topologies indicated in the two text files.

### 03 Visualisation
This is a collection of several scripts neccessary for the visualisation of the phylogeny, biogeography and sequence motif evolution:
- **[Plot the phylogenetic tree](03_Visualisation/Readme.md#plot-phylogram-with-ggtree)** and **[project it on a geographic map](03_Visualisation/Readme.md#plot-phylogeny-projected-on-map)**
- **[Project the pairwise distance in a heatmap](03_Visualisation/Readme.md#pairwise-distance-heatmap)**
- **[Visualise the mitochondrial control region](03_Visualisation/Readme.md#control-region-visualisation)**
- **[Plot the sequence motif evolution](03_Visualisation/Readme.md#motif-evolution)**
