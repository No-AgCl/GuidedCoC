#  Guided Co-clustering Transfer Across Unpaired and Paired Single-cell Multi-omics Data

## Introduction
GuidedCoC is an unsupervised framework for structure-guided transfer learning from unpaired scRNA-seq to paired scRNA-seq/scATAC-seq data, enabling joint clustering and cross-modal feature alignment.

![image](https://github.com/No-AgCl/GuidedCoC/blob/main/overview/GuidedCoC.png)

## Getting Started

This repository contains MATLAB and R scripts for **GuidedCoC**, a structure-guided co-clustering framework for integrating unpaired scRNA-seq and paired scRNA-seq/scATAC-seq data.  
The core algorithm is implemented in **MATLAB**, while **R** scripts are used for data visualization.

---

### Requirements

####  MATLAB(for the GuidedCoC model)

- Version: **R2019b** or later
- Required: **Parallel Computing Toolbox** (for multi-core acceleration)


> Note: The core algorithm supports **parallel computation** using `parfor` to accelerate clustering and matching steps.  
Please ensure the **Parallel Computing Toolbox** is installed and `parpool` can be initialized before running.



####  R (for the data visualization)

- R version ≥ 4.1.0  


###  Runtime Environment

The experiments were conducted on a high-performance workstation with the following specifications:

- **CPU**: Intel Core i7-14700K (water cooled)  
- **Memory**: 64 GB DDR5 5600 MHz  
- **GPU**: NVIDIA RTX 4090D 24 GB  
- **Power Supply**: 1250 W  
- **Operating System**: Linux 64-bit

> Note: The multi-core CPU and large memory significantly accelerate parallel operations in MATLAB.

The code has also been tested on Windows 10 with:

- MATLAB R2022b or later
- R 4.3.0 or later

##Main functions

- GuidedCoC.m: impletation of the GuidedCoC algorithm
- match.m: performs optimal label matching to predict target cell types based on source clusters by our criterion
- Eval.m: calculate and show the ARI and NMI values of clustering results
- visualization.R: display the heatmaps of data after clustering and matching by GuidedCoC

##  Datasets

We evaluated GuidedCoC on four real-world single-cell datasets, each consisting of an unpaired scRNA-seq source dataset and a paired scRNA-seq + scATAC-seq target dataset. All datasets are publicly available:

| Example | Source Data (scRNA-seq) | Target Data (paired RNA+ATAC) | Type | Link |
|---------|--------------------------|--------------------------------|------|------|
| **Ex1** | 5k PBMCs from 10x Genomics (Donor 1, RNA only) | 10x Genomics PBMC multiome (3k, RNA+ATAC) | Human PBMC | [Link](https://www.10xgenomics.com/datasets/5k_Human_Donor1_PBMC_3p_gem-x) / [Link](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-3-k-1-standard-1-0-0) |
| **Ex2** | E18 Mouse Brain scRNA-seq (10x) | E18 Mouse Brain Multiome (RNA+ATAC) | Mouse brain | [Link](https://www.10xgenomics.com/datasets/5k_Human_Donor1_PBMC_3p_gem-x) / [Link](https://www.10xgenomics.com/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-1-0-0) |
| **Ex3** | Mouse Lymph Node scRNA-seq (Tabula Muris) | Human Lymph Node Multiome (RNA+ATAC, 10x) | Cross-species | [Link](https://www.10xgenomics.com/datasets/Mixture-of-cells-from-mouse-lymph-nodes-and-spleen-stained-with-totalseqc-mouse-universal-cocktail) / [Link](https://www.10xgenomics.com/datasets/fresh-frozen-lymph-node-with-b-cell-lymphoma-14-k-sorted-nuclei-1-standard-1-0-0) |
| **Ex4** | Human Pancreatic Islet scRNA-seq (HPAP) | Human Islet Multiome (RNA+ATAC, 10x) | Human pancreas | [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133) / [Link]( https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200044) |

## Examples

Please check the <a href="vignette/vignette.md"><u>vignette</u></a> for a tutorial. A few examples are contained for a quick start of GuidedCoC. 

 

##Acknoledgement

We utilized and modified part of the code of CoupleCoC+ algorithm, therefore please refer to the CoupleCoC+ code at [https://github.com/cuhklinlab/coupleCoC_plus](https://github.com/cuhklinlab/coupleCoC_plus).


The baseline methods used in this study are referenced from the benchmarking study published in *Nature Biotechnology*, titled **"Benchmarking algorithms for joint integration of unpaired and paired single-cell RNA-seq and ATAC-seq data"**.  
[https://pubmed.ncbi.nlm.nih.gov/37875977/](https://pubmed.ncbi.nlm.nih.gov/37875977/)