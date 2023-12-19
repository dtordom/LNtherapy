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


##--------------------------------------------------- Step2
## PLOTS



##--------------------------------------------------- Step3
## GSEA (similarity between Gene-signatures)



##--------------------------------------------------- Step4
## Functional analysis (qusage)













