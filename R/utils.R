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
      clin.tmp.x<-clin.tmp[,c("Patient_ID","Response")]
      clin.tmp.x<-clin.tmp.x[!duplicated(clin.tmp.x$Patient_ID),]
      clin.tmp.x$Patient_ID<-paste0("P",clin.tmp.x$Patient_ID)
      
      #nYES<-nrow(clin.tmp[clin.tmp$Response=="YES",])
      #nNO<-nrow(clin.tmp[clin.tmp$Response=="NO",])
      
      samples<-intersect(rownames(var.table),clin.tmp.x$Patient_ID)
      var.x<-var.table[samples,];
      rownames(clin.tmp.x)<-clin.tmp.x$Patient_ID
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
      clin.tmp.x<-clin.tmp[,c("Patient_ID","Response")]
      clin.tmp.x<-clin.tmp.x[!duplicated(clin.tmp.x$Patient_ID),]
      clin.tmp.x$Patient_ID<-paste0("P",clin.tmp.x$Patient_ID)
      
      nYES<-nrow(clin.tmp.x[clin.tmp.x$Response=="YES",])
      nNO<-nrow(clin.tmp.x[clin.tmp.x$Response=="NO",])
      
      samples<-intersect(rownames(var.table),clin.tmp.x$Patient_ID)
      var.x<-var.table[samples,];
      rownames(clin.tmp.x)<-clin.tmp.x$Patient_ID
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

#################################################################
## Function to plot expression of markers in single cell clusters

markerDots<-function(seurat,
                     genes,
                     wd=1000,
                     ht=680,
                     stroke=0.3,
                     path){
  
  require('pheatmap')
  metadata<-seurat@meta.data
  metadata$cell_id<-rownames(metadata)
  
  #genes<-str_to_title(genes)
  genes<-intersect(genes,rownames(seurat@assays$RNA@data))
  expression_matrix <- seurat@assays$RNA@data[genes,,drop=F]
  
  cell_markers <- expression_matrix %>%
    Matrix::t() %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    pivot_longer(cols = genes, names_to = "gene", values_to = "expression") %>%
    left_join(metadata[,c("cell_id","seurat_clusters")])
  
  count_per_clusters <- unique(cell_markers[,c("cell_id","seurat_clusters")]) %>%
    group_by(seurat_clusters) %>%
    dplyr::summarise(cells = n())
  
  percentage_expression <- cell_markers %>%
    mutate(bool = expression > 0) %>%
    filter(bool) %>%
    group_by(seurat_clusters,gene) %>%
    dplyr::summarise(count = n()) %>%
    left_join(count_per_clusters) %>%
    mutate(percentage = (count / cells)*100)
  
  # Escalar la expresión por gen
  expression_vals <- rbindlist(mclapply(unique(cell_markers$gene), function(gene_to_use){
    expression_vals <- cell_markers[cell_markers$gene == gene_to_use,]
    expression_vals$scale_vals <- scale(expression_vals$expression)
    expression_vals
  })) %>%
    # Calcular la mediana de los valores escalados
    group_by(seurat_clusters,gene) %>%
    dplyr::summarise(median_expression = median(scale_vals)) %>%
    left_join(percentage_expression)
  
  expression_vals$gene <- factor(expression_vals$gene, levels = unique(cell_markers$gene))
  m<-expression_vals[,c("seurat_clusters","percentage","gene")]
  m[is.na(m)]<-0; m<-as.data.frame(m)
  
  clusters<-as.numeric(names(table(m$seurat_clusters)))
  tmp<-matrix(ncol=length(clusters),nrow=length(genes))
  colnames(tmp)<-clusters; rownames(tmp)<-genes
  
  for(i in 1:length(genes)){
    for(j in clusters){
      punt<-j+1
      tmp[i,punt]<-m[m$seurat_clusters==clusters[punt] & m$gene==genes[i],]$percentage
    }
  }
  cl<-pheatmap(tmp,cluster_rows = F)
  cl<-cl$tree_col$order
  
  expression_vals$seurat_clusters<-factor(x = expression_vals$seurat_clusters,
                                          levels=colnames(tmp)[cl])
  
  tiff(filename = path,width = wd,height = ht,res = 300)
  gg4 <- ggplot(expression_vals, aes(x = seurat_clusters, y = gene))+
    geom_point(aes(fill = median_expression, size = percentage),shape = 21,stroke=stroke)+
    theme_classic()+coord_flip()+
    theme(axis.title = element_blank(),
          legend.key.size = unit(0.2,"cm"),
          legend.title = element_text(size = 5),
          legend.text = element_text(size = 5),
          axis.text.x = element_text(size = 5.8,angle = 90,hjust = 1),
          axis.text.y = element_text(size = 5.8,angle = 0,hjust = 1),
          axis.ticks.x = element_line(colour = 'black', size = 0.3),
          axis.ticks.y = element_line(colour = 'black', size = 0.3),
          axis.line = element_line(colour = 'black', size = 0.3))+
    scale_fill_gradient(low = "white", high = "darkred")+
    scale_size(range=c(0.5,3.3))+
    labs(fill = "Expression", size = "Percentage of\npositive cells")
  plot(gg4)
  invisible(dev.off())
  return(gg4)
}



