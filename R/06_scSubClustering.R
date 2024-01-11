##############################
## Cellcept: Subclustering
## R version 4.1.0
## 
##############################
## Example for Cluster 1
## Aproximate computational time (for bigger cluster): 

library("Seurat")
library("pheatmap")

DATA.i<-readRDS("Cluster_1.rds") ## Load RDS with clustered cells. Change PATH to file
dataset_metadata<-DATA.i@meta.data

##------------------------------------------------------------ Step 1
## Clustering

## Define clustering resolution
All.PostSCT <- FindNeighbors(object = DATA.i,reduction = "pca",assay = "RNA",dims = 1:40)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.01)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.05)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.1)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.2)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.3)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.4)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.5)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.6)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.7)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.8)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.9)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 1)

unsup.clust.colors <- pal_igv("defaul", alpha = 1)(30)
clustree(All.PostSCT, prefix = "RNA_snn_res.")

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


##------------------------------------------------------------ Step 3
## Get Gene-signature expression of response/non-response across clusters

degs<-readRDS("/geneListDEGs.RData") ##Change PATH to file

## degs can be any list of individual genes or gene-signatures
## I.e. Signatures to functionally annotate the clusters:

# degs<-list("DN4"= c("HOPX","PDE4D","IGHE","SELL"),
#            "DN3" = c("RHOB","VIM","CTSH"),
#            "DN2"= c("EMP3","CIB1","PSAP","CD72","DAPP1","HCK","ZEB2","RHOB","TNFRSF1B","FCRL3","FCRL5","FGR","MPP6"),
#            "DN1" = c("TAGLN2","IGHA2","JCHAIN","IGHA1","S100A10"),
#            "Cmem2" = c("LTB","TAGLN2","AHNAK","ITGB1","CRIP1","S100A10"),
#            "Cmem1" = c("LTB","MT-ATP8"),
#            "Mmem1" = c("LTB","IGHM"),
#            "Mmem2" = c("LTB","VIM","IGHM","AHNAK",""),
#            "naive" = c("SELL","PLPP5","FCER2","ILR4","IGHM"),
#            "trans" = c("VPREB3","IGHD","IIGLL5","TCL1A"))


PLOTS<-list()
count<-1
maxs<-NULL
for(i in 1:length(degs)){ 
  x<-degs[[i]]
  for(j in 1:length(x)){
    features<-intersect(x[[j]],rownames(DATA.i@assays$SCT))
    
    pbmc <- AddModuleScore(DATA.i,features = list(features),
                           name= paste0(names(degs)[i],sep="_",names(x)[j]),
                           nbin = length(features),
                           assay = "RNA",ctrl = 10)
    
    p1<-FeaturePlot(pbmc,features = paste0(names(degs)[i],sep="_",names(x)[j],"1"),
                    label = TRUE, repel = TRUE,raster = FALSE) +
      scale_colour_gradient2(
        low = muted("#132B43"),
        mid = "#FAFAFA",
        high = muted("#67041f"),
        midpoint = 0,
        space = "Lab",
        na.value = "#FAFAFA",
        guide = "colourbar",
        aesthetics = "colour")
    
    metadata<-pbmc@meta.data; metadata<-metadata[,c("seurat_clusters",paste0(names(degs)[i],sep="_",names(x)[j],"1"))]
    colnames(metadata)<-c("seurat_clusters","value")
    
    metadata<- metadata %>%
      group_by(seurat_clusters) %>%
      dplyr::summarise(median_expression = median(value)) 
    
    m<-as.data.frame(metadata[order(metadata$median_expression,decreasing = T),])
    colnames(m)<-c("cluster","score")
    
    maxs<-c(maxs,max(m$score))
    names(maxs)[length(maxs)]<- paste0(names(degs)[i],sep="_",names(x[j])) 
    
    p2<-ggplot(data=m,mapping = aes(x=cluster, y=score,fill=cluster) ) + 
      geom_bar(stat="identity",colour="black",size=0.3) + theme_classic()+
      theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))+
      geom_hline(yintercept =  median(m$score),linetype="dashed",color="black",alpha=0.6,size = 0.5)+
      scale_fill_manual(values = ClustersColors) +ggtitle(paste0(names(degs)[i],sep="_",names(x)[j]))
    
    plots<-list(p1,p2); names(plots)<-c("umap","barr")
    
    PLOTS[[count]]<-plots
    names(PLOTS)[count]<-paste0(names(degs)[i],sep="_",names(x)[j])
    count<-count+1
  }
}



