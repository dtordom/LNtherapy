######################################################################
## Apply lm models using limma
#' @param data Gene-expression matrix. Samples in columns and genes in rows
#' @param metadata Dataframe with metadata. Samples in rows and variables in 
#' columns (rownames in metadata and colnames in data must be the same)
#' @param covars Formula to adjust the model (i.e. ~sex+age+act+dosis+prednisone)
#' @param padj Method to adjust p-value (i.e. "fdr", "bonferroni","BH",..)

limma.DEG<-function(data,
                    metadata,
                    covars,  
                    padj="bonferoni"){
  
  require("limma")
  
  data<-data[,rownames(metadata)]
  design.disease <- model.matrix(~Response,data = metadata)
  design.covariates <- model.matrix(covars,data = metadata) 
  
  ebatch <- removeBatchEffect(data[,rownames(design.covariates)],
                              covariates = design.covariates,
                              design = design.disease[rownames(design.covariates),])
  
  DEG <- eBayes(lmFit(ebatch, model.matrix(~Response,data=metadata[colnames(ebatch),])))
  DEG<-topTable(DEG,number = nrow(ebatch),adjust.method = padj)
  
  return(DEG)
}

######################################################################
## Function to extract sample order in a heatmap
#' @param exp Gene-expression matrix. Samples in columns and genes in rows
#' @param met Dataframe with metadata. Samples in rows and variables in 
#' columns (rownames in metadata and colnames in data must be the same)
#' @param var variable taken into account to sort samples

orderHeatmap<-function(exp,met,var){
  vars<-names(table(met[,var]))
  ordenPats<-NULL
  for(i in 1:length(vars)){    
    exp.i<-exp[,rownames(met[met[,var]==vars[i],])]
    pl<-pheatmap::pheatmap(as.matrix(exp.i),scale="row",show_colnames = F,cluster_cols = T,
                 cluster_rows = T, show_rownames = T, border_color = NA,
                 breaks=seq(-1.5,1.5,length.out = 100), fontsize = 5,
                 color = colorRampPalette(c("deepskyblue4","white","coral2"))(100))  
    ordenPats<-c(ordenPats,colnames(exp.i[pl$tree_col$order]))
  }
  return(ordenPats)
}


######################################################################
## Function to perform clinical association analysis
GetStats<-function(clin.tmp,var.table,bySamples=TRUE){
  RESULTS<-list()
  res.num<-do.call("rbind",lapply(1:ncol(var.table),function(x){
    
    if(bySamples){
      samples<-intersect(rownames(clin.tmp),rownames(var.table))
      clin.tmp.x<-clin.tmp[samples,]
      var.x<-var.table[samples,x]
      resp<-na.omit(var.x[clin.tmp.x$Response=="YES"])
      noresp<-na.omit(var.x[!clin.tmp.x$Response=="YES"])
      
      #nYES <- nrow(clin.tmp[clin.tmp$Response=="YES"])
      #nNO <- nrow(clin.tmp[clin.tmp$Response!="YES"])
      
    }else{
      clin.tmp.x<-clin.tmp[,c("patient","Response")]
      clin.tmp.x<-clin.tmp.x[!duplicated(clin.tmp.x$patient),]
      
      #nYES<-nrow(clin.tmp[clin.tmp$Response=="YES",])
      #nNO<-nrow(clin.tmp[clin.tmp$Response=="NO",])
      
      samples<-intersect(rownames(var.table),clin.tmp.x$patient)
      var.x<-var.table[samples,];
      rownames(clin.tmp.x)<-clin.tmp.x$patient
      clin.tmp.x<-clin.tmp.x[samples,]
      
      resp<-na.omit(var.x[clin.tmp.x$Response=="YES",x])
      noresp<-na.omit(var.x[clin.tmp.x$Response!="YES",x])
    }
    
    
    if(length(resp)>3 & length(noresp)>3){
      
      if(class(resp)!="character"){
        
        results<-data.frame("pvalue" = wilcox.test(resp,noresp)$p.value,
                            "mean_resp" = mean(resp,na.rm=T),
                            "sd_resp" = sd(resp,na.rm=T),
                            "mean_noresp" = mean(noresp,na.rm=T),
                            "sd_noresp" = sd(noresp,na.rm=T))
        rownames(results)<-colnames(var.table)[x]
        return(results)
        print(results)
      }
    }
  }))
  
  ###..................
  
  res.cat<-do.call("rbind",lapply(1:ncol(var.table),function(x){
    
    if(bySamples){
      samples<-intersect(rownames(clin.tmp),rownames(var.table))
      clin.tmp.x<-clin.tmp[samples,]
      var.x<-var.table[samples,x]
      resp<-na.omit(var.x[clin.tmp.x$Response=="YES"])
      noresp<-na.omit(var.x[!clin.tmp.x$Response=="YES"])
      
      nYES <- nrow(clin.tmp[clin.tmp$Response=="YES",])
      nNO <- nrow(clin.tmp[clin.tmp$Response!="YES",])
      
    }else{
      clin.tmp.x<-clin.tmp[,c("patient","Response")]
      clin.tmp.x<-clin.tmp.x[!duplicated(clin.tmp.x$patient),]
      
      nYES<-nrow(clin.tmp.x[clin.tmp.x$Response=="YES",])
      nNO<-nrow(clin.tmp.x[clin.tmp.x$Response=="NO",])
      
      samples<-intersect(rownames(var.table),clin.tmp.x$patient)
      var.x<-var.table[samples,];
      rownames(clin.tmp.x)<-clin.tmp.x$patient
      clin.tmp.x<-clin.tmp.x[samples,]
      
      resp<-var.x[clin.tmp.x$Response=="YES",x]
      noresp<-var.x[clin.tmp.x$Response!="YES",x]
    }
    
    if(length(resp)>3 & length(noresp)>3){
      
      if(class(resp)=="character"){
        
        vars<-na.omit(unique(c(resp,noresp)))
        
        if(length(vars)>1){
          
          tabl<-matrix(data=0,ncol=2,nrow=length(vars))
          colnames(tabl)<-c("resp","noresp"); rownames(tabl)<-vars
          
          for(v in vars){
            tabl[v,1]<-length(resp[resp==v])
            tabl[v,2]<-length(noresp[noresp==v])
          }
          
          pvalue = min(fisher.test(tabl,alternative = "less")$p.value,
                       fisher.test(tabl,alternative = "greater")$p.value)
          
          results<-data.frame("pvalue" = pvalue,
                              "n_resp" = tabl["YES","resp"],
                              "perc_resp" = (tabl["YES","resp"] / nYES)*100,
                              "n_noresp" = tabl["YES","noresp"],
                              "perc_noresp" = (tabl["YES","noresp"] / nNO)*100)
          rownames(results)<-colnames(var.table)[x]
          return(results)
          #print(results)
          
        } else {
          ## valor NA, guardar props
          results<-data.frame("pvalue" = 1,
                              "n_resp" = length(resp[resp=="YES"]),
                              "perc_resp" = (length(resp[resp=="YES"]) / nYES)*100,
                              "n_noresp" = length(noresp[noresp=="YES"]),
                              "perc_noresp" = (length(noresp[noresp=="YES"]) / nNO)*100)
          rownames(results)<-colnames(var.table)[x]
          return(results)
          p#rint(results)
        }
        
      }
    }
  }))
  
  RESULTS[["cat"]]<-res.cat
  RESULTS[["num"]]<-res.num
  return(RESULTS)
}



