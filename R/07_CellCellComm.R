##############################
## Cell-Cell comunication using CellChat
## R version 4.1.0
## 
##############################
#devtools::install_github("sqjin/CellChat")
## Aproximate computational time (testing with 25% of randomly selected cells):
## NOTE: include all cells has high computational cost

setwd("/rdata/")  ## Change PATH to single cell RData files of each cluster/ cell type 

library("CellChat")
library("stringi")
library("stringr")
library("NMF")
library("biomaRt")

##------------------------------------------------------------ Step 1
## Load data

files<-list.files()
files<-files[str_detect(files,"m_")] ## sc data from each cell type clustered 

## Load data from clusters
DATA<-NULL
MET<-NULL
for(i in 1:length(files)){
  
  load(files[i])
  ## Use these commented code to test in a subsets of randomly selected cells
  #randomcells<-sample(1:ncol(data),ncol(data)*0.4) 
  #data<-data[,randomcells]
  met<-met[colnames(data),]
  
  if(i==1){
    DATA<-data
    MET<-met
  }else{
    MET<-rbind(MET,met)
    sel<-intersect(rownames(DATA),rownames(data))
    DATA<-cbind(DATA[sel,],data[sel,])
  }
}
data<-DATA
met<-MET

##------------------------------------------------------------ Step 2
## Impute cell-cell communications
cellchat <- createCellChat(object = data, meta = met, group.by = "seurat_clusters")
future::plan("multisession", workers = 12) # do parallel

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB ## using all database
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

## Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

##------------------------------------------------------------ Step 3
## See all cell-cell communicationd and select manually networks (paths) that connect clusters relevant for non-response to one or another drug.
paths<-c("SEMA4", "APP", "BAG", "GAS", "FLT3", "BAFF", "IL1", "CCL") 

## CCL: AZA
## IL1: AZA
## BAG
## BAFF: MMF
## FLT3: MMF (...)
## GAS: MMF (...)
## APP: MMF ---
## SEMA4: MMF

## Using biomaRt to translate gene ids from entrez to gene symbol (gene targets form each network/path)
mart = useMart("ensembl", dataset = paste0(casefold("hsapiens"),"_gene_ensembl"),host="https://jul2019.archive.ensembl.org")
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="https://jul2019.archive.ensembl.org")
genome <- getBM(attributes=c('external_gene_name','entrezgene_id'),
                mart = ensembl)
genome <- na.omit(genome)

DRUG<-NULL
for(i in 1:length(paths)){ 
  x<-cellchat@DB$interaction[cellchat@DB$interaction$pathway_name==paths[i],]
  entrz<-unique(c(x$receptor,x$ligand))
  tmp<-genome[genome$external_gene_name %in% entrz,]
  tmp<-paste(tmp$entrezgene_id,collapse=",") 
  tmp<-as.data.frame(c(paths[i],tmp))
  DRUG<-rbind(DRUG,t(tmp))
}
rownames(DRUG)<-NULL
colnames(DRUG)<-c("target","entrez")

saveRDS(DRUG,"DRUGS.rds")

