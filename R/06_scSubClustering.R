##############################
## Cellcept: Subclustering
## R version 4.1.0
## 
##############################
## Example for Cluster 1

library("Seurat")
library("pheatmap")

##------------------------------------------------------------ Step 1
## Load RDS with clustered cells

DATA.i<-readRDS("Cluster1.rds")
dataset_metadata<-DATA.i@meta.data
