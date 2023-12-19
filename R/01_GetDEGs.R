##########################
## Get DEG genes beetwen responders/ non-responders
## R version 4.3.1 (2023-06-16)
##########################

## Load Gene expression data and metadata from:
## https://drive.google.com/drive/folders/15KVn3tckVPZKieCsSZizdi8bAQFZkOIC?usp=sharing

set.seed(123456)
library("pheatmap")
library("fgsea")
library("ggplot2")
library("stringr")
library("UpSetR")
library("ggrepel")
library("ggpubr")
library("cowplot")
library("ComplexHeatmap")

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

## UpSet plot
lt = list(MMF = unname(as.character(unlist(degs$mmf))),
          AZA = unname(as.character(unlist(degs$aza))),
          HC = unname(as.character(unlist(degs$hc))),
          SOC = unname(as.character(unlist(degs$soc))))
m1 = make_comb_mat(lt,mode = "distinc")
UpSet(m1, set_order = c("MMF", "AZA", "HC","SOC"), comb_order = order(comb_size(m1),decreasing = T))

## Volcano plot
PLOTS<-list()
for(j in 1:length(DEGS)){
  volcan<-DEGS[[j]]
  volcan$P.Value<-(-log10(volcan$P.Value))
  volcan$adj.P.Val<-(-log10(volcan$adj.P.Val))
  volcan$genes<-rownames(volcan)
  
  volcan$color<-ifelse(volcan$adj.P.Val>1.30103,"#619bc2",
                       ifelse(volcan$P.Value>1.30103,"#86cdeb","#a8a8a8"))
  volcan$label<-volcan$genes
  for(i in 1:nrow(volcan)){
    if(volcan$adj.P.Val[i]<1.30103){
      volcan$label[i]=""
    }
  }
  p1<-ggplot(data=volcan, aes(x=logFC, y=P.Value, color=color,label=label)) + 
    geom_point(size=0.7) + theme_classic()+ 
    xlim(min(volcan$logFC)+(min(volcan$logFC)*0.15),max(volcan$logFC)+(max(volcan$logFC)*0.15))+
    geom_text_repel(color="black",size=1.75,point.padding = 1)+
    scale_colour_manual(values = c("#619bc2","#86cdeb","#a8a8a8"))
  PLOTS[[j]]<-p1
}
names(PLOTS)<-names(DEGS)

## Heatmaps



## Ratios
PLOTS<-list()
for(i in 1:length(DEGS)){
  
  x<-DEGS[[i]]
  up<-x[x$adj.P.Val<=0.05 & x$logFC > 0,]
  down<-x[x$adj.P.Val<=0.05 & x$logFC < 0,]
  
  plots<-list()
  for(j in 1:length(DATA)){
   tmp <-data.frame("group"=DATA[[j]]$clin[colnames(DATA[[j]]$data),"Response"],
                 "ratio"=apply(DATA[[j]]$data[rownames(up),],2,mean) / 
                   apply(DATA[[j]]$data[rownames(down),],2,mean),
                 "drug"=rep(names(DATA)[j],ncol(DATA[[j]]$data)))
   
   p.tmp<-ggplot(tmp,aes(x=group,y=ratio)) + geom_jitter(color = "#333333",size=0.7) + 
     geom_boxplot(outlier.alpha = 0,width = 0.5,fill = "lightgrey",alpha=0.8,fatten=0.5,lwd=0.5) + 
     theme_classic() + 
     theme(axis.line = element_line(linewidth=0.5),
           axis.ticks = element_line(linewidth = 0.5))+
     ggtitle(wilcox.test(tmp$ratio[tmp$group=="YES"],tmp$ratio[tmp$group!="YES"])$p.value) +
     theme(plot.title = element_text(size = 4)) + ylim(0.7,1.2) 
   plots[[j]]<-p.tmp
  }
  names(plots)<-names(DATA)

  pbox<-ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],
                  ncol = 4, nrow = 1,legend = "none")
  PLOTS[[i]]<-pbox
}
names(PLOTS)<-c("mmf","aza","hc","soc")



##--------------------------------------------------- Step3
## GSEA (similarity between Gene-signatures)

## 1. Comparison between drugs
genes<-list(
  "MMF.up"= rownames(DEGS$MMF[DEGS$MMF$logFC>0 & DEGS$MMF$adj.P.Val<=0.05,]),
  "MMF.down"= rownames(DEGS$MMF[DEGS$MMF$logFC<0 & DEGS$MMF$adj.P.Val<=0.05,]),
  "AZA.up" = rownames(DEGS$AZA[DEGS$AZA$logFC>0 & DEGS$AZA$adj.P.Val<=0.05,]),
  "AZA.down" = rownames(DEGS$AZA[DEGS$AZA$logFC<0 & DEGS$AZA$adj.P.Val<=0.05,]),
  "HC.up" = rownames(DEGS$HC[DEGS$HC$logFC>0 & DEGS$HC$adj.P.Val<=0.05,]),
  "HC.down" = rownames(DEGS$HC[DEGS$HC$logFC<0 & DEGS$HC$adj.P.Val<=0.05,]),
  "SOC.up" = rownames(DEGS$SOC[DEGS$SOC$logFC>0 & DEGS$SOC$adj.P.Val<=0.05,]),
  "SOC.down" = rownames(DEGS$SOC[DEGS$SOC$logFC<0 & DEGS$SOC$adj.P.Val<=0.05,]))

resfgsea<-lapply(DEGS,function(x){
  tmp<- x[order(x$logFC,decreasing = T),]
  stats<-tmp$logFC; names(stats)<-rownames(tmp)
  res<-fgsea(pathways = genes,stats = stats)
  return(res)})

M<-do.call("rbind",lapply(names(resfgsea),function(x){
  x<-cbind(rep(x,nrow(resfgsea[x][[1]])),resfgsea[x][[1]])}))
colnames(M)[1]<-"Drug"
M<-as.data.frame(M[,c("Drug","pathway","ES","NES")])
M<-M[!str_detect(M$pathway,M$Drug),]

## PLOT
ggplot(M,aes(x=pathway,y=Drug))+theme_linedraw()+
  geom_point(aes(colour = ES,size = abs(NES))) + 
  theme(axis.text.x = element_text(size =8,angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_colour_gradient2(low = "#2786a0",mid="white",high = "#c7522b") + 
  geom_point(aes(size=abs(NES)),color="black",shape=1)


## 2. Comparison between results using different response metrics (UrCrPrRatio and SRI-4)
MMFresp2<-readRDS("MMF_resp2.rds")
res<-limma.DEG(data=MMFresp2$MMF_resp2$data,
               metadata = MMFresp2$MMF_resp2$clin,
               covars = ~sex+act+mmfDosis+prednisone,
               padj = "bonferroni")

res<-res[order(res$logFC,decreasing = T),]
stats<-res$logFC; names(stats)<-rownames(res)

genes<-degs$mmf
res<-fgsea(pathways = genes,stats = stats)

## PLOT
p1<-plotEnrichment(genes$nR_down,stats) 
p2<-plotEnrichment(genes$nR_up,stats)
p1$layers[[1]]$aes_params$colour <- '#c7522b'  
p1$layers[[3]]$aes_params$colour <- '#c7522b'
p2$layers[[1]]$aes_params$colour <- '#2786a0'
p2$layers[[3]]$aes_params$colour <- '#6dbc86'
p1
p2

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
