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




