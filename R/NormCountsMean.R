#' @title NormCountsMean
#'
#' @description Calculation mean counts for two contrast or all matrix.
#'
#' @param glmfit fitted linear model object.
#' @param ASKOlist list of data.frame contain condition, contrast and context information made by asko3c.
#' @param context coefficient/contrast tested.
#' @return one of this:
#' \itemize{
#'    \item matrix with mean counts,
#'    \item meanValue for one context/Condition.
#' }
#'
#' @examples
#' \dontrun{
#'     # calculate mean counts in contrast contx1_vs_contx2
#'     mean1<-NormCountsMean(glmfit, ASKOlist, contx1)   # in the 1st context
#'     mean2<-NormCountsMean(glmfit, ASKOlist, contx2)   # in the 2nd context
#'
#'     # for all conditions
#'     n_count<-NormCountsMean(glmfit, ASKOlist, context=NULL)
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
NormCountsMean <- function(glmfit, ASKOlist, context=NULL){
  lib_size_norm<-glmfit$samples$lib.size*glmfit$samples$norm.factors                          # normalization computation of all library sizes
  if(is.null(context)==TRUE){
    set_condi<-row.names(ASKOlist$condition)
  }else{
    set_condi<-ASKOlist$context$condition[ASKOlist$context$context==context]                  # retrieval of condition names associated to context
  }
  table_c_norm <- data.frame(row.names = row.names(glmfit$counts))

  for (condition in set_condi){
    sample_name<-rownames(glmfit$samples[glmfit$samples$condition==condition,])               # retrieval of the replicate names associated to conditions
    subset_counts<-data.frame(row.names = row.names(glmfit$counts))                           # initialization of data frame as subset of counts table
    for(name in sample_name){
      lib_sample_norm<-glmfit$samples[name,"lib.size"]*glmfit$samples[name,"norm.factors"]    # normalization computation of sample library size
      subset_counts$c<-glmfit$counts[,name]                                                   # addition in subset of sample counts column
      subset_counts$c<-subset_counts$c*mean(lib_size_norm)/lib_sample_norm                    # normalization computation of sample counts
      colnames(subset_counts)[colnames(subset_counts)=="c"]<-name                             # to rename the column with the condition name
    }

    mean_counts<-rowSums(subset_counts)/ncol(subset_counts)                                   # computation of the mean
    table_c_norm$m <- mean_counts

    if(is.null(context)==TRUE){
      colnames(table_c_norm)[colnames(table_c_norm)=="m"]<-condition
    }else{
      colnames(table_c_norm)[colnames(table_c_norm)=="m"]<-paste(context,condition,sep = "/")
    }
  }

  return(table_c_norm)
}
