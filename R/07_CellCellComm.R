
#devtools::install_github("sqjin/CellChat")
## Cell-Cell comunication using CellChat

setwd("/rdata/")

library("CellChat")
library("stringi")
library("stringr")
library("NMF")
library("biomaRt")

files<-list.files()
files<-files[str_detect(files,"m_")] ## sc data from each cell type clustered 

DATA<-NULL
MET<-NULL
for(i in 1:length(files)){
  
  load(files[i])
  #print(i)
  #random25<-sample(1:ncol(data),ncol(data)*0.4)
  #data<-data[,random25]
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

paths<-c("SEMA4", "APP", "BAG", "GAS", "FLT3", "BAFF", "IL1", "CCL") 

## CCL: AZA
## IL1: AZA
## BAG
## BAFF: MMF
## FLT3: MMF (...)
## GAS: MMF (...)
## APP: MMF ---
## SEMA4: MMF

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


############################################

pathways.show <- unique(df.net$pathway_name)

for(i in 1: length(pathways.show)){
  netVisual_chord_cell(cellchat, signaling = pathways.show[i],
                       #targets.use = targets,
                       #sources.use = targets,
                       #group = group.cellType, 
                       title.name =pathways.show[i])
}

paths<-c("CXCL","IL1","IL16","SEMA4")

targets<-c("Bcell_2","Bcell_4","Mono1_2","Mono1_4","Mono1_6","Mono6_1","Mono6_2",
           "NK_3","NK_4")

netVisual_chord_cell(cellchat, signaling = pathways.show,
                     targets.use = c("Bcell_2","Mono1_2","Mono1_6","NK_3"),
                     sources.use = 
                     #group = group.cellType, 
                     title.name ="All")
