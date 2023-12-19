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

