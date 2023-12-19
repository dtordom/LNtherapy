##########################
## Association between drug response and clinical variables
## R version 4.3.1 (2023-06-16)
##########################

set.seed(123456)

## Load Data
## https://drive.google.com/drive/folders/15KVn3tckVPZKieCsSZizdi8bAQFZkOIC?usp=sharing
load("ClinicalData.RData")
load("DEGS.RData") ## Output of 01_GetDEGs.R

##--------------------------------------------------- Step1
## Get significance

## SLEDAI
res.sledai<-lapply(1:length(DATA),function(drug){
  res<-GetStats(clin.tmp = DATA[[drug]]$clin,
                       var.table = SLEDAI)
  return(res)
})
names(res.sledai)<-names(DATA)

## PGA
res.pga<-lapply(1:length(DATA),function(drug){
  res<-GetStats(clin.tmp = DATA[[drug]]$clin,
                var.table = PGA)
  return(res)
})
names(res.pga)<-names(DATA)

## TREATS
res.treats<-lapply(1:length(DATA),function(drug){
  res<-GetStats(clin.tmp = DATA[[drug]]$clin,
                var.table = TREATS)
  return(res)
})
names(res.treats)<-names(DATA)

## OTHERS
res.others<-lapply(1:length(DATA),function(drug){
  res<-GetStats(clin.tmp = DATA[[drug]]$clin,
                var.table = OTHERS,
                bySamples = T)
  return(res)
})
names(res.others)<-names(DATA)

## SERO
res.sero<-lapply(1:length(DATA),function(drug){
  res<-GetStats(clin.tmp = DATA[[drug]]$clin,
                var.table = SERO,
                bySamples = F)
  return(res)
})
names(res.sero)<-names(DATA)

## OTHER2
res.other2<-lapply(1:length(DATA),function(drug){
  res<-GetStats(clin.tmp = DATA[[drug]]$clin,
                var.table = OTHER2[,c("Sex.UC","Age.Onset","FAMHX")],
                bySamples = F)
  return(res)
})
names(res.other2)<-names(DATA)

## WEIGHT
res.we<-lapply(1:length(DATA),function(drug){
  res<-GetStats(clin.tmp = DATA[[drug]]$clin,
                var.table = WEIGHT,
                bySamples = F)
  return(res)
})
names(res.we)<-names(DATA)


## RACE 
race<-OTHER2[,c("Race.UC","Race.UC","Race.UC")]
colnames(race)<-c("AA","C","O")
race$AA<-ifelse(race$AA=="AA","YES","NO")
race$C<-ifelse(race$C=="C","YES","NO")
race$O<-ifelse(race$O=="O","YES","NO")

res.race<-lapply(1:length(DATA),function(drug){
  res<-GetStats(clin.tmp = DATA[[drug]]$clin,
                var.table = race,
                bySamples = F)
  return(res)
})
names(res.race)<-names(DATA)
