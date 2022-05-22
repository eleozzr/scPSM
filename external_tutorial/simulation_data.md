Load necessary R packages.

``` r
library(splatter)
library(dplyr)
library(Seurat)
```

# 1. Data of Simulation 1

This simulation data is related to Fig 2.

Set parameters: 1000 genes, one cell type, two batches, 100 cells for
each batch, almost no batch effect, and additional dropout to the whole
dataset.

``` r
params <- newSplatParams()
params <- setParam(params, "nGenes", 1000)
params <- setParam(params, "batchCells", c(100,100))
params <- setParam(params, "batch.facLoc", 0.001)
params <- setParam(params, "batch.facScale", 0.001)
params <- setParam(params, "dropout.type", "experiment") 
params <- setParam(params, "out.prob", 0) # no outlier
params <- setParam(params, "dropout.mid",3) 
params <- setParam(params, "seed", 1) 
```

Create true counts without dropout.

``` r
sim <- splatSimulate(params, method="single", verbose=FALSE)
# set counts of the SCE object equal to true counts
sim0 <- sim
counts(sim0) <- assays(sim0)$TrueCounts
# this normalization is no use, which will be replaced by NormalizeData() in Seurat object
sim0 <- scater::logNormCounts(sim0) 
input0 <- Seurat::as.Seurat(sim0)
input0 <- NormalizeData(input0)
```

Create counts with dropout in Batch2

``` r
sim1 <- sim
# set counts of Batch 1 equal to true counts, while keep dropout in Batch2
counts(sim1)[,1:100] <- assays(sim1)$TrueCounts[,1:100]
sim1 <- scater::logNormCounts(sim1) # no use
input <- Seurat::as.Seurat(sim1)
input <- NormalizeData(input)
```

Save data as the input for integration methods.

``` r
save(input, file = "raw1.rda")
SeuratDisk::SaveH5Seurat(input, filename = "raw1.h5Seurat")
SeuratDisk::Convert("raw1.h5Seurat", dest = "h5ad")
```

# 2. Data of Simulation 2

This simulation data is related to Fig 2.

Set parameters: 1000 genes, two cell types, two batches,
batch.factor=0.8, de.probability=0.3.

``` r
params <- newSplatParams()
params <- setParam(params, "nGenes", 1000)
params <- setParam(params, "batchCells", c(50,150)) # batch assignment
params <- setParam(params, "batch.facLoc", 0.8)
params <- setParam(params, "batch.facScale", 0.8)
params <- setParam(params, "group.prob", c(2/3,1/3)) # cell type assignment
params <- setParam(params, "de.prob",c(0.3,0.3))
params <- setParam(params, "out.prob", 0) # no outlier
params <- setParam(params, "seed", 1) 
```

Create counts in Seurat object.

``` r
sim <- splatSimulate(params, method="groups", verbose=FALSE)

input <- sim %>%
  scater::logNormCounts() %>% # no use
  as.Seurat() %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 200)

table(input@meta.data$Batch, input@meta.data$Group)
```

    ##         
    ##          Group1 Group2
    ##   Batch1     31     19
    ##   Batch2    103     47

Save data as the input for integration methods.

``` r
save(input, file = "raw2.rda")
SeuratDisk::SaveH5Seurat(input, filename = "raw2.h5Seurat")
SeuratDisk::Convert("raw2.h5Seurat", dest = "h5ad")
```

# 3. Data of supplementary simulation

This simulation data is related to Fig S4.

Set parameters: batch.factor = 0.01, 0.05, 0.1, 0.5. Create counts in
Seurat object.

``` r
batch.loc <- batch.scale <- c(0.01, 0.05, 0.1, 0.5)
sim <- lapply(1:4, function(x){
  params <- newSplatParams()
  params <- setParam(params, "nGenes", 1000)
  params <- setParam(params, "batchCells", c(100,100))
  params <- setParam(params, "batch.facLoc", batch.loc[x])
  params <- setParam(params, "batch.facScale", batch.scale[x])
  params <- setParam(params, "seed", 1) 
  res <- splatSimulate(params, method="single", verbose=FALSE) %>%
    scater::logNormCounts() %>% 
    Seurat::as.Seurat() 
  return(res)
})
```
