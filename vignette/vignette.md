# A quick guide to GuidedCoC and baseline methods

## 1. Implementation of GuidedCoC

```matlab
clear 
clc 
close all

%% Load the processed data in example 1 in the paper.
% Note: The "data" folder includes all the preprocessed datasets (ex1–ex3) used in the paper.

%%%%%%%%%%%%%%%%%%%%% GuidedCoC algorithm %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Input %%%%%%%%%%%%%%%
% p, q and q0 represent data R<sup>(s)</sup>, data R<sup>(t)</sup>, and A<sup>(t)</sup>, respectively
% nrowcluster: number of cell clusters in data R<sup>(t)</sup> and A<sup>(t)</sup>
% ncolcluster: number of feature clusters in data R<sup>(s)</sup>, R<sup>(t)</sup>, and A<sup>(t)</sup>
% iter: number of iterations (default: 20)
% beta: weight of source data R<sup>(s)</sup>
% alpha: weight of unlinked features in A<sup>(t)</sup>
% Cx_truth: known ground truth labels of data R<sup>(s)</sup>

%%%%% Output %%%%%%%%%%%%%%%
% Cy: predicted cell clusters for target data
% Cz: feature clusters for both source and target data
% cluster_p, cluster_q, cluster_q0: joint distributions for cells/features in R<sup>(s)</sup>, data R<sup>(t)</sup>, and A<sup>(t)</sup>
% obj: vector of objective values per iteration
% match_result: optimal mapping between source and target clusters
% matm: all mapping trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('data/processed_data/ex1_source_data_RNA.mat');
load('data/processed_data/target_data/ex1_target_data_RNA.mat');
load('data/processed_data/ex1_target_data_ATAC.mat');
load('data/processed_data/ex1_source_data_cell_label.mat'); 
load('data/processed_data/ex1_target_data_cell_label.mat'); 
 
% Hyperparameter settings
nrowcluster = 8; 
ncolcluster = 12; 
iter = 10; 
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

## 2. Heatmaps of clustering results by GuidedCoC

```r
source("R/visualization.R") 
library('gplots') 
library('R.matlab')

scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(50)

Y_link <- as.matrix(readMat("data/graph_data/q_log.mat")[[1]])
U <- as.matrix(readMat("data/graph_data/q0.mat")[[1]])
CY_best <- as.vector(readMat("data/graph_data/Cy.mat")[[1]])
CZ_best <- as.vector(readMat("data/graph_data/Cz.mat")[[1]])

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

## 3. One example for using baseline methods in the paper

### MultiVI

```python
# TODO: Add MultiVI usage example
```

### Cobolt

```python
# TODO: Add Cobolt usage example
```

### scMoMaT

```python
# TODO: Add scMoMaT usage example
```

### Seurat v4

```python
# TODO: Add Seurat example
```





