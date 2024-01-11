##----------------------------------------------------------------------- STEP 1
## Load data

library("hipathia")

load("DEGS.RData")  ## Change PATH to file

H<-data[,rownames(metadata)[metadata$Diagnosis!="SLE"]]
MMF<-DATA$MMF$data
AZA<-DATA$AZA$data

drugs<-readRDS("/DRUGS.rds") ## Output 07_CellCellComm. ## Change PATH to file
drugs<-as.data.frame(drugs)
colnames(drugs)<-c("Target","Entrez")

##----------------------------------------------------------------------- STEP 2
## Hypathia - Inference of circuit activity

## Get circuits activity (Hipathia)
pathways <- load_pathways(species = "hsa") 

colnames(H)<-paste0(colnames(H),"_healthy")
colnames(MMF)<-paste0(colnames(MMF),"_mmf")
colnames(AZA)<-paste0(colnames(AZA),"_aza")

## Annotation
data<-as.matrix(cbind(H,MMF,AZA))
data<- translate_data(data,species = "hsa")

## Create matrix for each condition
data.MMF <- data[,str_detect(string = colnames(data),pattern = "_mmf")]
data.AZA <- data[,str_detect(string = colnames(data),pattern = "_aza")]
data.H <- data[,str_detect(string = colnames(data),pattern = "_healthy")]

data.SLE<-cbind(data.MMF,data.AZA)

## Add drug modifiers
drugs.matrix<-list()
for(i in 1:nrow(drugs)){
  
  tmp<-data.SLE
  targets<-as.character(str_split(drugs$Entrez[i],",")[[1]])
  targets<-intersect(rownames(tmp),targets)
  
  if(length(targets)>0){
    tmp[targets,]<-tmp[targets,]*0.1 ## Inhibition
    colnames(tmp)<-paste0(colnames(tmp),"_",drugs[i,"Target"])
    drugs.matrix[[drugs[i,"Target"]]]<-tmp
  }else{
    print(paste0("NO targets shared for: ",drugs$Target[i]))
  }
}

colnames(data.SLE)<-paste0(colnames(data.SLE),"_base")

all.data<-cbind(data.H,data.SLE,do.call("cbind",drugs.matrix))
all.data.norm<-normalize_data(all.data)

## Get circuits activity
tmp <- hipathia(all.data.norm, pathways) 
data.ALL.paths <- get_paths_data(tmp, matrix = TRUE) 
All.pathways<- normalize_paths(data.ALL.paths, metaginfo = pathways)



##----------------------------------------------------------------------- STEP 3
## Hypathia - Get response score by patient


## Calculate distance between two vectors
CalculateEuclideanDistance <- function(vect1, vect2){
  
  distances<-NULL
  for(i in 1:ncol(vect1)){
    distances<-c(distances,sqrt(sum((vect1[,i] - vect2[,i])^2)))
  }
  names(distances)<-colnames(vect1)
  return(distances)
} 

H<-All.pathways[,str_detect(string = colnames(All.pathways) ,pattern = "_healthy")]
BASE<-All.pathways[,str_detect(string = colnames(All.pathways) ,pattern = "_base")]
colnames(BASE)<-gsub("_base","",colnames(BASE))


#
base<-BASE


clinResp<-data.frame("id"= c(DATA$mmf.prcr$clin$sampleID,DATA$aza.sri4$clin$sampleID),
                     "drug"= c(rep("mmf",nrow(DATA$mmf.prcr$clin)),rep("aza",nrow(DATA$aza.sri4$clin))),
                     "response" = c(DATA$mmf.prcr$clin$Response,DATA$aza.sri4$clin$Response))

Results.table<-NULL
Results.plots<-NULL
Results.freq<-NULL
Results.means<-NULL

for(i in 1:length(drugs.matrix)){
  
  tmp.data<-NULL
  tmp.data<-All.pathways[,str_detect(string = colnames(All.pathways) ,pattern = paste0("_",sep="",names(drugs.matrix)[i]))]
  #print(tmp.data[1:10,1:4])
  colnames(tmp.data)<-gsub(paste0("_",sep="",names(drugs.matrix)[i]),"",colnames(tmp.data))
  
  x<-CalculateEuclideanDistance(vect1 = tmp.data, vect2 = base) ## Here, changes in transcriptome are calculated as a distance between after and before target inhibition (for each patient)
  
  names(x)<-gsub("_aza",replacement = "",names(x))
  names(x)<-gsub("_mmf",replacement = "",names(x))
  x<-x[clinResp$id]
  
  tmp<-data.frame(clinResp,"dist"=as.numeric(x),
                  "group"=paste0(clinResp$drug,sep="_",clinResp$response))
  
  # tmp<-data.frame("patients"=names(x),"dist"=as.numeric(x),
  #                 "cluster"=ids[names(x),"cluster"]) 
  #print(tmp)
  ## tmp$dist values of change in transcriptome for all patients
  
  pl1<-data.frame(table(tmp[tmp$dist>=mean(tmp$dist),"group"])[c("mmf_NO","mmf_YES","aza_NO","aza_YES")] / table(tmp$group)[c("mmf_NO","mmf_YES","aza_NO","aza_YES")])
  rownames(pl1)<-pl1$Var1
  pl1$Var1<-1:4
  pl1[is.na(pl1)]<-0
  
  Results.freq<-rbind(Results.freq,pl1$Freq)
  #pl1$Freq <- pl1$Freq*100
  
  pl1<-pl1[str_detect(string = rownames(pl1),"_NO"),]
  pl1$Var1<-rownames(pl1)
  
  p1<-ggplot(pl1,aes(x=Var1,y=Freq))+geom_bar(stat = "identity")+ theme_bw() + ylim(0,1)+ggtitle(names(drugs.matrix)[i])
  #p2<-ggplot(tmp,aes(x=group,y=dist))+geom_jitter()+geom_boxplot(alpha=0.2)+theme_classic()+ggtitle(names(drugs.matrix)[i])
  
  res<-wilcox.test(tmp[tmp$group=="mmf_NO","dist"],tmp[tmp$group=="aza_NO","dist"])$p.value
  
  print(paste0(names(drugs.matrix)[i],sep="_",res))
  plot(p1)
 
}

names(Results.plots)<-names(drugs.matrix)
rownames(Results.freq)<-names(drugs.matrix)
rownames(Results.means)<-names(drugs.matrix)
colnames(Results.freq)<-c("Cl1.freq","Cl2.freq","Cl3.freq")
colnames(Results.means)<-c("Cl1.mean(sd)","Cl2.mean(sd)","Cl3mean(sd)")

Results.table<-na.omit(Results.table)
Results.plots<-Results.plots[rownames(Results.table)]
Results.freq<-Results.freq[rownames(Results.table),]

colnames(Results.table)<-paste0("pval_",sep="",colnames(Results.table))

Results.freq<-round(Results.freq*100,digits = 2)

Results<-cbind(Results.freq,Results.table)

Results = data.frame("Target"=rownames(Results),Results)
Results.means<-Results.means[rownames(Results),]
Results.means = data.frame("Target"=rownames(Results.means),Results.means)

write.table(Results,"hypathiaResults.txt",sep="\t",quote = T,row.names = F)
write.table(Results.means,"hypathiaMeans.txt",sep="\t",quote = T,row.names = F)

for(i in 1:length(Results.plots)){
  plot(Results.plots[[i]]$boxplot)
}

for(i in 1:length(Results.plots)){
  plot(Results.plots[[i]]$barrplot)
}

Results.table

