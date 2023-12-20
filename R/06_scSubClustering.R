##############################
## Cellcept: Subclustering
## R version 4.1.0
## 
##############################
## Example for Cluster 1

library("Seurat")
library("pheatmap")

DATA.i<-readRDS("Cluster_1.rds") ## Load RDS with clustered cells
dataset_metadata<-DATA.i@meta.data

##------------------------------------------------------------ Step 1
## Clustering

resolution<-0.11 ## Set cluster resolution
cat("\nPerforming Clustering")
DATA.i <- FindNeighbors(object = DATA.i,reduction = "pca",dims = 1:40)
DATA.i <- FindClusters(DATA.i, resolution = resolution, algorithm = 1) ## celltype
DATA.i <- RunUMAP(DATA.i,dims = 1:40, n.components = 2, verbose = FALSE, future.seed = NULL)

print(table(DATA.i@meta.data$seurat_clusters))
plot(DimPlot(DATA.i, reduction = "umap", label = TRUE, pt.size = 0.00001,label.size = 3.5))

##------------------------------------------------------------ Step 2
## Heatmap top10 markers
top<-10

clusters<-unique(as.character(DATA.i@meta.data$seurat_clusters))
DATA_markers <-getDEG(DATA.i,filterAdj = T,
                      selectClusters=clusters,
                      nameFile="/Cluster_1_comparisons.csv")

DATA_markers %>% group_by(cluster) %>% top_n(top,avg_log2FC) -> top10
genes<-unique(top10$gene)
m <- DATA.i@assays$RNA@counts
m<-m[genes,rownames(DATA.i@meta.data)[DATA.i@meta.data$seurat_clusters %in% clusters]]

x<-DATA_markers
sel<-ifelse(x$p_val_adj<=0.05,T,F)
x<-x[sel,]
write.table(x, row.names = F,file = "/Cluster_1_comparisons.csv",
            sep="\t",quote = FALSE)

clusterpos<-as.character(DATA.i@meta.data[colnames(m),"seurat_clusters"])
mplot<-matrix(data=0,ncol=length(clusters),nrow=length(genes))
colnames(mplot)<-clusters; rownames(mplot)<-genes
for(i in 1:length(genes)){
  for(j in 1:length(clusters)){
    tmp<-m[i,clusterpos==clusters[j]]
    mplot[i,j]<-mean(as.numeric(tmp),na.r=T)
}}
colors<-brewer.pal(n = 9, name = 'Blues')

p1<-pheatmap(t(mplot),scale="column",cluster_rows = T,cluster_cols = T,
             color=c("white",colors), breaks = seq(-1.6,1.6,length.out = 9),
             clustering_method = "complete",fontsize=4.5,border_color="black")



