# Multi-batch integration with scPSM

## Overview

This repository contains the code for the paper **Propensity Score Matching enables batch effect corrected imputation in single-cell RNA-seq analysis** by Xu _et al._

ScPSM (a **p**ropensity **s**core **m**atching method for **sc**RNA-seq data) is a statistical tool useful for simultaneously correcting batch effect, imputing dropout and denoising gene expression.

## Arguments

- **`batches`** A list of one or more log-normalized-expression matrices where genes correspond to rows and cells correspond to columns. Each matrix should contain the same number of rows, corresponding to the same genes in the same order. Each matrix represents a batch.
- **`markers`** A vector specifying which features used as marker genes to compute propensity scores.
- **`hvg`** A vector specifying which features used as HVGs for identifying MNN group.
- **`k.self`** An integer scalar specifying the number of nearest neighbors in searching KNNs.
- **`k.mnn`** An integer scalar specifying the number of nearest neighbors in matching MNN pairs.
- **`correct.all`** A logical scalar specifying whether correction should be applied to all genes, even if only a subset is used for the MNN group identification.
- **`merge.order`** An integer vector containing the linear merge order of batches.

## Value

Returns a dgCMatrix of integrated data with rows corresponding to the same genes and columns corresponding to the same cells as in the argument **`batches`**.

## Usage

To perform scPSM, first run the source file [`scPSM_utils.R`](./script/scPSM_utils.R), then run the script [`scPSM_main.R`](./script/scPSM_main.R).

Please refer to a toy example. To run the example require the software *R >= 4.0.0*, *batchelor >= 1.4.0*, *BiocNeighbors >= 1.6.0* and *BiocParallel >= 1.22.0*.

## Data
The original data for the toy example is available in the **`data`** folder
- Load [`pancreas_expression_matrix.rds`](./data/pancreas_expression_matrix.rds) for a dgCMatrix with rows corresponding to 34363 genes and columns corresponding to 6321 cells.
- Load [`pancreas_metadata.rds`](./data/pancreas_metadata.rds) to get the batch (the "tech" item) and celltype information for all cells .
- [`HVGs_1000.txt`](./data/HVGs_1000.txt) can be extracted from *adata.var.index[adata.var["highly_variable"] == True].values* by implementing the python function *sc.pp.highly_variable_genes(adata, n_top_genes=1000, batch_key='tech')* by *impoting `scanpy` as sc*, or from *obj[["RNA"]]@var.features* by implementing the R function *FindVariableFeatures(obj, nfeatures = 1000)* by *library(`Seurat`)*.

