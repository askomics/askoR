NormCountsMean<-function(glmfit, ASKOlist, context){
  
  lib_size_norm<-glmfit$samples$lib.size*glmfit$samples$norm.factors                          # normalization computation of all library sizes 
  set_condi<-ASKOlist$context$condition[ASKOlist$context$context==context]                    # recovery of condition names associated to context
  
  for (condition in set_condi){
    sample_name<-rownames(glmfit$samples[glmfit$samples$condition==condition,])               # recovery of the replicate names associated to conditions
    subset_counts<-data.frame(row.names = row.names(glmfit$genes))                            # initialization of data frame as subset of counts table
    for(name in sample_name){
      lib_sample_norm<-glmfit$samples[name,"lib.size"]*glmfit$samples[name,"norm.factors"]    # normalization computation of sample library size 
      subset_counts$c<-glmfit$counts[,name]                                                   # addition in subset of sample counts column
      subset_counts$c<-subset_counts$c*mean(lib_size_norm)/lib_sample_norm                    # normalization computation of sample counts
      colnames(subset_counts)[colnames(subset_counts)=="c"]<-name                             # to rename the column with the condition name
    }
    mean_counts<-rowSums(subset_counts)/ncol(subset_counts)                                   # computation of the mean
    ASKOlist$stat.table$mean<-mean_counts                                                     # subset integration in the glm_result table 
    colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="mean"]<-paste(context,condition,sep = "/")             
  }                                                                                           # to rename the column with the context name
  return(ASKOlist$stat.table)                                                                 # return the glm object
}
