# A quick guide to GuidedCoC

## 1. Data Preparatioon

The data folder contains all datasets used in the four examples of the paper. It is organized into two subfolders:

- data/processed data/: Contains the processed RNA and ATAC matrices, along with their corresponding cell labels. There are the input matrices (p, q, q0) to the GuidedCoC algorithm.

- data/graph data/: Contains the clustering results and the matrices used for visualization (e.g. heatmaps). These include Cx_truth, Cy, Cz, and the p, q, q0 matrices used for plotting.

Each example (ex1–ex4) has its own set of files following consistent naming:

- *_source_data_RNA.mat,  *_(labeled)_target_data_RNA.mat,  *_(labeled)_target_data_ATAC.mat: matrices used as input to the algorithm and used for plotting heatmaps.

- *_source_data_cell_label.mat,  *_target_data_cell_label.mat: ground truth labels.

- *_Cx_truth.mat, *_Cy.mat, *_Cz.mat: clustering results for plotting.


### Note on Example4 Data

Due to GitHub’s file size limitations, three large matrices in data/graph data/ for Example 4 are not included in this repository. Specifically:

data/graph data/ex4_source_data_RNA.mat (contains variable p)

data/graph data/ex4_labeled_target_data_RNA.mat (contains variable q)

data/graph data/ex4_labeled_target_data_ATAC.mat (contains variable q0)

These matrices are **only used for visualization and evaluation**, not for running the GuidedCoC algorithm. You can regenerate them with the following steps.

#### Step 1: Regenerate q and q0 for Graph Data
You can either run the provided script **R/ex4_labeled_cells_selected.R**, or copy and run the following R code:

```r
library(R.matlab)
library(data.table)


base_dir <- "data/processed data"
input_dir <- "data/processed data"
output_dir <- "data/graph data"


all_barcodes <- fread(file.path(base_dir, "ex4_target_data_barcodes.txt"), header = FALSE)[[1]]
valid_barcodes <- fread(file.path(base_dir, "ex4_labeled_cells_barcodes.txt"), header = FALSE)[[1]]


row_idx <- which(all_barcodes %in% valid_barcodes)
cat("Extracting", length(row_idx), "rows for labeled cells\n")


process_and_save <- function(mat_path, row_idx, var_name, out_filename) {
  cat(" Processing:", mat_path, "\n")
  

  mat_list <- readMat(mat_path)
  
  if (!(var_name %in% names(mat_list))) {
    stop(paste("Variable", var_name, "not found in", mat_path))
  }
  
  mat <- mat_list[[var_name]]
  

  mat_labeled <- mat[row_idx, , drop = FALSE]
  

  out_path <- file.path(output_dir, out_filename)
  

  args <- list()
  args[[var_name]] <- mat_labeled
  do.call(writeMat, c(list(con = out_path), args))
  
  cat("Saved to:", out_path, "\n\n")
}


process_and_save(file.path(input_dir, "ex4_target_data_RNA.mat"), row_idx, var_name = "q", out_filename = "ex4_target_data_RNA.mat")
process_and_save(file.path(input_dir, "ex4_target_data_ATAC.mat"), row_idx, var_name = "q0", out_filename = "ex4_target_data_ATAC.mat")


```

#### Step 2: Copy p from Processed Data

```r
file.copy("data/processed data/ex4_source_data_RNA.mat", 
          "data/graph data/ex4_source_data_RNA.mat", overwrite = TRUE)

```

## 2. Implementation of GuidedCoC

```matlab
clear 
clc 
close all

% Load the processed data in example 1 in the paper.
% Note: The "data" folder includes the preprocessed datasets (ex1–ex4) used in the paper.

%%%%%%%%%%%%%%%%%%%%% GuidedCoC algorithm %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Input %%%%%%%%%%%%%%%
% p, q and q0 represent data R⁽ˢ⁾, data R⁽ᵗ⁾, and data A⁽ᵗ⁾, respectively
% nrowcluster: number of cell clusters in  data R⁽ᵗ⁾ and data A⁽ᵗ⁾
% ncolcluster: number of feature clusters in data R⁽ˢ⁾, data R⁽ᵗ⁾, and data A⁽ᵗ⁾
% iter: number of iterations (default: 20)
% beta: weight of source data R⁽ˢ⁾
% alpha: weight of unlinked features in A⁽ᵗ⁾
% Cx_truth: known ground truth labels of data R⁽ˢ⁾

%%%%% Output %%%%%%%%%%%%%%%
% Cy: predicted cell clusters for target data
% Cz: feature clusters for both source and target data
% cluster_p, cluster_q, cluster_q0: joint distributions for cells/features in Rˢ, data Rᵗ, and Aᵗ
% obj: vector of objective values per iteration
% match_result: optimal mapping between source and target clusters
% matm: the list of matching results from all shuffling trials.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('data/processed_data/ex1_source_data_RNA.mat');
load('data/processed_data/ex1_target_data_RNA.mat');
load('data/processed_data/ex1_target_data_ATAC.mat');
load('data/processed_data/ex1_source_data_cell_label.mat'); 
load('data/processed_data/ex1_target_data_cell_label.mat'); 
 
% Hyperparameter settings
nrowcluster = 8; 
ncolcluster = 12; 
iter = 20; 
beta = 1; 
alpha = 0.8;
ntrials = 25; 
jsd_threshold = 0.45;
n_shuffles = 25;
numCols = size(p,2);

addpath('Matlab'); % The Matlab folder contains function files for GuidedCoC

epsilon = 1e-6;
p_eps  = p  + epsilon;  
q_eps  = q  + epsilon;  
q0_eps = q0 + epsilon; 

% Run the GuidedCoC algorithm
[Cx_truth_new, Cy, Cz, cluster_p, cluster_q, cluster_q0, matm, match_result, obj] = ...
    GuidedCoC(p_eps, q_eps, q0_eps, Cx_truth, nrowcluster, ncolcluster, iter, beta, alpha, ntrials, jsd_threshold, n_shuffles, numCols);

% Evaluation
[TAB_Y, Eval_tab2] = Eval(Cy_truth, Cy); % Generates contingency table and metrics (ARI, NMI)
cell_sums = sum(TAB_Y, 1);
disp(match_result);
```

## 3. Heatmaps of clustering results by GuidedCoC

```r
source("R/visualization.R") 
library('gplots') 
library('R.matlab')

scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(50)

Y_link <- as.matrix(readMat("data/graph_data/ex1_target_data_RNA.mat")[[1]])
CY_best <- as.vector(readMat("data/graph_data/ex1_Cy.mat")[[1]])
CZ_best <- as.vector(readMat("data/graph_data/ex1_Cz.mat")[[1]])

CY_best <- reorder_label(CY_best, 1:8)
CZ_best <- reorder_label(CZ_best, 1:12)

gene_order_Z <- clu_sep(CZ_best)
cell_order_Y <- clu_sep(CY_best)
colsep_ST <- line_sep(CZ_best)
rowsep_Y <- line_sep(CY_best)

raw_data_y <- strong_signal_iter(Y_link[cell_order_Y, gene_order_Z], 
                                 CY_best, 
                                 clu_num(CY_best), 
                                 n_pick = 20, 
                                 n_iter = 5)

heatmap_fun_rna(raw_data_y, scaleyellowred, colsep_ST, rowsep_Y)
```

![alt text](https://github.com/No-AgCl/GuidedCoC/blob/main/images/heatmap_ex1.png)







