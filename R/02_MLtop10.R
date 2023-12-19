##########################
## Get ML-based models to predict Response to LN drugs
## R version 4.3.1 (2023-06-16)
##########################

set.seed(12345678)
library("pathMED")
library("stringr")
library("parallel")
library("doParallel")

## ······································································· Step 1 
## Load data
load("DEGS.RData") ## Output of 01_GetDEGs.R

degMMF<-DEGS$MMF[DEGS$MMF$adj.P.Val<=0.05,]
degAZA<-DEGS$AZA[DEGS$AZA$adj.P.Val<=0.05,]
degHC<-DEGS$HC[DEGS$HC$adj.P.Val<=0.05,]
degSOC<-DEGS$SOC[DEGS$SOC$adj.P.Val<=0.05,]

degMMF<-degMMF[!str_detect(rownames(degMMF), str_c(c("LOC","LINC","-AS1","orf"), collapse = "|")),]
degAZA<-degAZA[!str_detect(rownames(degAZA), str_c(c("LOC","LINC","-AS1","orf"), collapse = "|")),]
degHC<-degHC[!str_detect(rownames(degHC), str_c(c("LOC","LINC","-AS1","orf"), collapse = "|")),]
degSOC<-degSOC[!str_detect(rownames(degSOC), str_c(c("LOC","LINC","-AS1","orf"), collapse = "|")),]

degMMF<-degMMF[rownames(degMMF[order(abs(degMMF$logFC),decreasing = T),])[1:10],]
degAZA<-degAZA[rownames(degAZA[order(abs(degAZA$logFC),decreasing = T),])[1:10],]
degHC<-degHC[rownames(degHC[order(abs(degHC$logFC),decreasing = T),])[1:10],]
degSOC<-degSOC[rownames(degSOC[order(abs(degSOC$logFC),decreasing = T),])[1:10],]

dataMMF<-DATA$MMF$data[rownames(degMMF),]
dataAZA<-DATA$AZA$data[rownames(degAZA),]
dataHC<-DATA$HC$data[rownames(degHC),]
dataSOC<-DATA$SOC$data[rownames(degSOC),]

clinMMF<-DATA$MMF$clin
clinAZA<-DATA$AZA$clin
clinHC<-DATA$HC$clin
clinSOC<-DATA$SOC$clin

## ······································································· Step 2 
## Sub-sampling using same patients for train/test and ML using getML (pathMED)

ML<-list("MMF"=list(dataMMF,clinMMF),
         "AZA"=list(dataAZA,clinAZA),
         "HC"=list(dataHC,clinHC),
         "SOC"=list(dataSOC,clinSOC))

ML.results<-list()
for(i in 1:length(ML)){
  data.tmp<-ML[[i]][[1]]
  clin.tmp<-ML[[i]][[2]]
  resp<-unique(clin.tmp[clin.tmp$Response=="YES",]$Patient_ID)
  noresp<-unique(clin.tmp[clin.tmp$Response=="NO",]$Patient_ID)
  
  subpats<-unname((lapply(1:5,function(x){
    res<-c(resp[sample(1:length(resp),size = round(length(resp)*0.7,digits = 0))],
           noresp[sample(1:length(noresp),size = round(length(noresp)*0.7,digits = 0))])
    res<-res[order(res)]})))
  
  clin.tmp$pnt<-1:nrow(clin.tmp)
  subsamples<-unname(lapply(subpats,function(x){
    res<-clin.tmp[clin.tmp$Patient_ID %in% x,]$pnt}))
  
  cl <- makePSOCKcluster(8)
  registerDoParallel(cl)
  start_time <- Sys.time()
  ML.results[[i]] <- getML(expData=data.tmp,
                           metadata=clin.tmp,
                           var2predict="Response",
                           models = methodsML(outcomeClass = "character"),
                           Koutter = subsamples,
                           Kinner = 4,
                           repeatsCV = 10,
                           continue_on_fail = TRUE,
                           positiveClass = "YES",
                           filterFeatures = TRUE,
                           saveLogFile = "logFit.txt")
  stopCluster(cl)
  registerDoSEQ()
  print(Sys.time() - start_time) 
}
names(ML.results)<-names(DATA)

rm(list=setdiff(ls(),c("ML.results")))
save.image("MLresults.RData")

