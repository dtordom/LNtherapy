##########################
## Comparison of cell proportions between responder and non-responders
## R version 4.3.1 (2023-06-16)
##########################
## Get associations between cell proportions and response/non respose to LN drugs

library("ggplot2")
linrary("ggpubr")
library("pheatmap")
set.seed(123456788)

## ······································································· Step 1 
## Impute cells from Cibersort
load("DEGS.RData") ## Output of 01_GetDEGs.R

## Download from google drive, see README file
source("/cibersort/cibersort.r") 
lmd22 = read.table('/cibersort/LM22.txt', sep = '\t', header = T, row.names = 1)

lmd22<-lmd22[,c("B.cells.naive","B.cells.memory","Plasma.cells","T.cells.CD8","T.cells.CD4.naive",
                "T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper",
                "T.cells.regulatory..Tregs.","T.cells.gamma.delta","NK.cells.resting","NK.cells.activated",
                "Monocytes","Dendritic.cells.resting","Dendritic.cells.activated","Eosinophils","Neutrophils",
                "Macrophages.M0","Macrophages.M1","Macrophages.M2")]

results.cbsrt <- CIBERSORT(sig_matrix = lmd22,
                     mixture_file = data,
                     perm=1000,QN=TRUE,
                     absolute = F,
                     abs_method = 'sig.score')

## Remove cells with poor quality
results<-results.cbsrt[ifelse(results.cbsrt[,22]<0.65,F,T),1:20]

allpatients<-NULL
for(i in 1:length(DATA)){
  tmp<-data.frame("ids" = rownames(DATA[[i]]$clin),
                  "drug" = rep(names(DATA)[i],nrow(DATA[[i]]$clin)),
                  "response" = DATA[[i]]$clin$Response)
  allpatients<-rbind(allpatients,tmp)
}
allpatients<-allpatients[!duplicated(allpatients$ids),]
rownames(allpatients)<-allpatients$ids

sel<-intersect(rownames(allpatients),rownames(results))
allpatients<-allpatients[sel,]
results.all<-as.data.frame(t(results[sel,]))

results.all<-results.all[apply(results.all,1,sd)!=0,]

## ······································································· Step 2
## Comparing cell proportion between responder and non-responders

PLOTS<-list()
for(i in 1:length(DATA)){
  
  clin.tmp<-DATA[[i]]$clin[intersect(rownames(DATA[[i]]$clin),rownames(results)),]
  cells<-data.frame("Response"=clin.tmp$Response,
                    results[rownames(clin.tmp),])
  
  pvals<-NULL
  count<-1
  plotList<-list()
  for(tc in 2:ncol(cells)){
    p<-tc
    nr<-as.numeric(cells[cells$Response=="NO",tc])*100
    r<-as.numeric(cells[cells$Response!="NO",tc])*100
    nr<-nr[ifelse(nr==0,F,T)]
    r<-r[ifelse(r==0,F,T)]
    
    if(length(nr)>2 & length(r)>2){
      p1<-wilcox.test(nr,r,alternative = "greater")
      p2<-wilcox.test(nr,r,alternative = "less")
      p<-min(p1$p.value,p2$p.value)
      pvals<-c(pvals,p)
    }else{
      pvals<-c(pvals,NA)
    }
    ## Plot
    if(p<=0.05){
      
      m<-data.frame("Response"=c(rep("NO",length(nr)),rep("YES",length(r))),
                    "value"=c(nr,r))
      m$value<-round(m$value,digits = 1)
      m$color<-ifelse(m$Response=="NO","#8f9aa2","#f0c35e")
      
      plotList[[count]]<-ggplot(m,aes(x=Response,y=value,fill=color))+theme_classic()+
        geom_jitter(size=0.5,color=m$color,alpha=0.8)+theme(axis.title.x = element_blank())+
        theme(axis.title.y = element_blank(),axis.text.y = element_text(size=6),
              axis.ticks = element_line(size=0.25),
              axis.line=element_line(size=0.25),
              plot.title = element_text(size = 4))+
        ggtitle(paste0(colnames(cells)[tc],sep=" ",p))+
        geom_boxplot(lwd=0.2,alpha=1,outlier.colour = NA,width = 0.5,alpha=0.8)+ scale_fill_manual(values=c("#8f9aa2","#f0c35e"))+
        theme(axis.text.x = element_blank())
      names(plotList)[count]<-colnames(cells)[tc]
      count<-count+1
    }
  }
  
  names(pvals)<-colnames(cells)[2:ncol(cells)]
  pvals<-pvals[ifelse(is.na(pvals)==T,F,T)]
  pvals<-pvals[pvals<=0.05]
  
  PLOTS[[i]]<-plotList
}
names(PLOTS)<-names(DATA)


## ······································································· Step 3
## Rich/poor analysis

rpPLOTS<-list()
for(dr in 1:length(DATA)){
  
  clin.tmp<-DATA[[dr]]$clin[intersect(rownames(DATA[[dr]]$clin),rownames(results)),]
  cells<-data.frame("Response"=clin.tmp$Response,
                    results[rownames(clin.tmp),])
  
  sel<-ifelse(apply(cells[2:ncol(cells)],2,sum)==0,F,T)
  cells<-cells[,c(TRUE,sel)]
  meanCells<-apply(cells[2:ncol(cells)],2,mean)
  
  rpcells<-do.call("cbind",lapply(names(meanCells),function(x){
    return(ifelse(cells[,x]>=meanCells[x],"Rich","Poor"))
  }))
  colnames(rpcells)<-names(meanCells); rownames(rpcells)<-rownames(cells)
  rpcells<-data.frame("Response" = cells$Response, rpcells)
  
  ## Get significance
  plotList<-list()
  count<-1
  pvals<-NULL
  for(i in 2:ncol(rpcells)){
    
    tmp<-data.frame("NO" = c(as.numeric(table(rpcells[ifelse(rpcells[,i]=="Poor",T,F),"Response"])["NO"]),
                             as.numeric(table(rpcells[ifelse(rpcells[,i]=="Rich",T,F),"Response"])["NO"])),
                    "YES" = c(as.numeric(table(rpcells[ifelse(rpcells[,i]=="Poor",T,F),"Response"])["YES"]),
                              as.numeric(table(rpcells[ifelse(rpcells[,i]=="Rich",T,F),"Response"])["YES"])))
    rownames(tmp)<-c("Poor","Rich")
    tmp[is.na(tmp)]<-0
    
    pval<-min(fisher.test(tmp,alternative = "less")$p.value,
              fisher.test(tmp,alternative = "greater")$p.value)
    pvals<-c(pvals,pval)
    if(pval<=0.05 & table(tmp>3)["TRUE"]==4){
      
      m<-data.frame("group"=rownames(tmp),
                    "value"=c(as.numeric((tmp["Poor","YES"]/apply(tmp,1,sum)[1])*100),
                              as.numeric((tmp["Rich","YES"]/apply(tmp,1,sum)[2])*100)))
      m$color<-ifelse(m$group=="Poor","#b38fc2","#85c37c")
      m$value<-round(m$value,digits = 2)
      
      m$texto<-c(paste0(tmp["Poor","YES"],sep="","/",apply(tmp,1,sum)[1]),
                 paste0(tmp["Rich","YES"],sep="","/",apply(tmp,1,sum)[2]))
      
      plotList[[count]]<-ggplot(m,aes(x=group,y=value,fill=color))+theme_classic()+
        theme(axis.title.x = element_blank())+
        theme(axis.title.y = element_blank(),axis.text.y = element_text(size=6),
              axis.line=element_line(size=0.4),
              plot.title = element_text(size = 4))+
        ggtitle(paste0(colnames(rpcells)[i],sep= " ",pval))+ylim(0,120)+
        geom_bar(stat = "identity")+ scale_fill_manual(values=c("#b38fc2","#85c37c"))+
        geom_text(aes(label = paste0(value,sep="","%")),size=1.4, vjust = -0.2,)+
        geom_text(aes(label = texto),size=1.4, vjust = 4,color="white",)+
        theme(axis.text.x = element_blank())
      
      names(plotList)[count]<-colnames(rpcells)[i]
      count<-count+1
      
    }
  }
  rpPLOTS[[dr]]<-plotList
  names(pvals)<-colnames(rpcells)[-1]  
}
names(rpPLOTS)<-names(DATA)




