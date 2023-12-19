##########################
## Get DEG genes beetwen responders/ non-responders
## R version 4.3.1 (2023-06-16)
##########################

## Load Gene expression data and metadata from:
## https://drive.google.com/drive/folders/15KVn3tckVPZKieCsSZizdi8bAQFZkOIC?usp=sharing

load("Dataset.RData")

DEGS<-list()
DATA<-list()

##--------------------------------------------------- Step1
## Get DEGs

## MMF ------
met.tmp<-metadata[metadata$Drug_group=="MMF" & metadata$Diagnosis!="Healthy",]
dat.tmp<-data[,rownames(met.tmp)]

res<-limma.DEG(data=dat.tmp,
                     metadata = met.tmp,
                     covars = ~Sex+Act+MMF_or_AZA_dosis+Prednisone_dosis,
                     padj = "bonferroni")
table(ifelse(res$adj.P.Val<=0.05,T,F))

data.list<-list(dat.tmp,met.tmp); names(data.list)<-c("data","clin")
DEGS[["MMF"]]<-res
DATA[["MMF"]]<-data.list

## AZA ------
met.tmp<-metadata[metadata$Drug_group=="AZA" & metadata$Diagnosis!="Healthy",]
dat.tmp<-data[,rownames(met.tmp)]

res<-limma.DEG(data=dat.tmp,
                     metadata = met.tmp,
                     covars = ~Sex+Act+MMF_or_AZA_dosis+Prednisone_dosis,
                     padj = "bonferroni")
table(ifelse(res$adj.P.Val<=0.05,T,F))

data.list<-list(dat.tmp,met.tmp); names(data.list)<-c("data","clin")
DEGS[["MMF"]]<-res
DATA[["MMF"]]<-data.list

## HC ------
met.tmp<-metadata[metadata$Drug_group=="HC" & metadata$Diagnosis!="Healthy",]
dat.tmp<-data[,rownames(met.tmp)]

res<-limma.DEG(data=dat.tmp,
                     metadata = met.tmp,
                     covars = ~Sex+Act,
                     padj = "bonferroni")
table(ifelse(res$adj.P.Val<=0.05,T,F))

data.list<-list(dat.tmp,met.tmp); names(data.list)<-c("data","clin")
DEGS[["HC"]]<-res
DATA[["HC"]]<-data.list

## SOC ------
met.tmp<-metadata[(metadata$Drug_group=="HC" | metadata$Drug_group=="PHC") & metadata$Diagnosis!="Healthy",]

## Remove patients with cytotoxic drugs
met.tmp<-met.tmp[!is.na(met.tmp$MMF_or_AZA_dosis),]
cyt<-unique(met.tmp[ifelse(met.tmp$MMF_or_AZA_dosis!=0,T,F),"Patient_ID"])
met.tmp<-met.tmp[!met.tmp$Patient_ID %in% cyt,]

dat.tmp<-data[,rownames(met.tmp)]

res<-limma.DEG(data=dat.tmp,
                     metadata = met.tmp,
                     covars = ~Sex+Act,
                     padj = "bonferroni")
table(ifelse(res$adj.P.Val<=0.05,T,F))

data.list<-list(dat.tmp,met.tmp); names(data.list)<-c("data","clin")
DEGS[["SOC"]]<-res
DATA[["SOC"]]<-data.list

rm(list=setdiff(ls(),c("DATA","DEGS","metadata","data")))

## Preparing Gene-signatures
degs<-list()

## MMF
tmp<-DEGS$MMF
up<-rownames(tmp)[ifelse(tmp$adj.P.Val<=0.05 & tmp$logFC>0,T,F)]
down<-rownames(tmp)[ifelse(tmp$adj.P.Val<=0.05 & tmp$logFC<0,T,F)]
mmf<-list(up,down); names(mmf)<-c("nR_down","nR_up") 
degs[["mmf"]]<-mmf

## AZA
tmp<-DEGS$AZA
up<-rownames(tmp)[ifelse(tmp$adj.P.Val<=0.05 & tmp$logFC>0,T,F)]
down<-rownames(tmp)[ifelse(tmp$adj.P.Val<=0.05 & tmp$logFC<0,T,F)]
aza<-list(up,down); names(aza)<-c("nR_down","nR_up")
degs[["aza"]]<-aza

## HC
tmp<-DEGS$HC
up<-rownames(tmp)[ifelse(tmp$adj.P.Val<=0.05 & tmp$logFC>0,T,F)]
down<-rownames(tmp)[ifelse(tmp$adj.P.Val<=0.05 & tmp$logFC<0,T,F)]
hc<-list(up,down); names(hc)<-c("nR_down","nR_up")
degs[["hc"]]<-hc

## SOC
tmp<-DEGS$SOC
up<-rownames(tmp)[ifelse(tmp$adj.P.Val<=0.05 & tmp$logFC>0,T,F)]
down<-rownames(tmp)[ifelse(tmp$adj.P.Val<=0.05 & tmp$logFC<0,T,F)]
soc<-list(up,down); names(soc)<-c("nR_down","nR_up")
degs[["soc"]]<-soc

saveRDS(degs,"geneListDEGs.RData")

rm(list=setdiff(ls(),c("DATA","DEGS","metadata","data","degs")))
save.image("DEGS.RData")

##--------------------------------------------------- Step2
## PLOTS (UpSet, Volcano, Heatmap, Ratios




##--------------------------------------------------- Step3
## GSEA (similarity between Gene-signatures)




##--------------------------------------------------- Step4
## Functional analysis (qusage)

library("qusage")
library("tmod")
library("stringi")

load("sysdata.rda") ## File contains tmod (and others) database
## Conection between genes and functional pathways

QS.results<-list()
for(i in 1:length(DATA)){
  
  tmp<-DATA[[i]]
  dat<-tmp$data
  met<-tmp$clin
  dat<-dat[,rownames(met)]
  label<-met$Response
  contrast = "YES-NO"
  
  Modules.list <- genesetsData$tmod
  Modules.ann <- as.data.frame(tmod$gs)
  Modules.ann<-Modules.ann[,c("ID","Title","SourceID")]
  colnames(Modules.ann)<-c("ID","Function","Database")
  rownames(Modules.ann)<-Modules.ann$ID
  Modules.list<-Modules.list[rownames(Modules.ann)]
  
  qs.results = qusage(dat, label, contrast, Modules.list)
  Modules.ann$pos = pdf.pVal(qs.results,alternative = "greater")
  Modules.ann$neg = pdf.pVal(qs.results,alternative = "less")
  
  paths<-Modules.ann[ifelse(Modules.ann$pos<=0.05 | Modules.ann$neg<=0.05,T,F),]
  #paths$pos<-(-log10(paths$pos)); paths$neg<-(-log10(paths$neg))
  paths<-paths[ifelse(paths$Function=="TBA" | paths$Function=="Undetermined",F,T),]
  #paths$ids<-paste0(paths$ID," -> ",paths$Function)
  paths<-paths[order(paths$pos,decreasing = F),]
  paths$drug<-rep(names(DATA)[i],nrow(paths))
  
  QS.results[[i]]<-paths
}
names(QS.results)<-names(DATA)

saveRDS(QS.results,"QusageResults.rds")
