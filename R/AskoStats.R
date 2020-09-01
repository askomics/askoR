#' @title AskoStats
#'
#' @description Based on result contained in "glm_test":
#' \itemize{
#'    \item print summary result of differential expression analysis
#'    \item format all results in tabulate out file followed parameters given
#'    \item plot heatmap of top differential expressed genes
#' }
#' Create one file by contrast, each file contains for each genes: fold-change and
#' log fold-change values, PValue, Expression, Significance, logCPM, LR, FDR and
#' significance value for each condition/context.
#' By default, LR and logCPM were not displayed, you can switch this parametres
#' to TRUE for display.
#'
#' @param glm_test, tests for one or more coefficients in the linear model (likelihood ratio tests or empirical Bayes quasi-likelihood F-tests).
#' @param fit, fitted linear model object.
#' @param contrast, coefficient/contrast names tested.
#' @param ASKOlist, list of data.frame contain condition, contrast and context informations made by asko3c.
#' @param dge, large DGEList with normalized counts by GEnorm function.
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @return none
#'
#' @examples
#' \dontrun{
#'     AskoStats(glm_test, fit, contrast, ASKOlist, dge, parameters)
#' }
#'
#' @export
AskoStats <- function (glm_test, fit, contrast, ASKOlist, dge, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  asko_dir = paste0(study_dir, "Askomics/")
  image_dir = paste0(study_dir, "images/")

  contrasko<-ASKOlist$contrast$Contrast[row.names(ASKOlist$contrast)==contrast]   # to retrieve the name of contrast from Asko object
  contx1<-ASKOlist$contrast$context1[row.names(ASKOlist$contrast)==contrast]      # to retrieve the name of 1st context from Asko object
  contx2<-ASKOlist$contrast$context2[row.names(ASKOlist$contrast)==contrast]      # to retrieve the name of 2nd context from Asko object

  ASKO_stat<-glm_test$table
  ASKO_stat$Test_id<-paste(contrasko, rownames(ASKO_stat), sep = "_")             # addition of Test_id column = unique ID
  ASKO_stat$contrast<-contrasko                                                   # addition of the contrast of the test
  ASKO_stat$gene <- row.names(ASKO_stat)                                          # addition of gene column = gene ID
  ASKO_stat$FDR<-stats::p.adjust(ASKO_stat$PValue, method=parameters$p_adj_method)       # computation of False Discovery Rate

  # Between context1 and context2 :
  ASKO_stat$Significance=0
  ASKO_stat$Significance[ASKO_stat$logFC < -parameters$threshold_logFC & ASKO_stat$FDR <= parameters$threshold_FDR] = -1  # Significance values = -1 for down regulated genes
  ASKO_stat$Significance[ASKO_stat$logFC > parameters$threshold_logFC  & ASKO_stat$FDR <= parameters$threshold_FDR] = 1   # Significance values =  1 for up regulated genes

  # addition of column "expression"
  ASKO_stat$Expression=NA
  ASKO_stat$Expression[ASKO_stat$Significance==-1]<-paste(contx1, contx2, sep="<")  # the value of attribute "Expression" is a string
  ASKO_stat$Expression[ASKO_stat$Significance==1]<-paste(contx1, contx2, sep=">")   # this attribute is easier to read the Significance
  ASKO_stat$Expression[ASKO_stat$Significance==0]<-paste(contx1, contx2, sep="=")   # of expression between two contexts

  if(parameters$Expression==TRUE){colg="Expression"}else{colg=NULL}
  if(parameters$logFC==TRUE){cola="logFC"}else{cola=NULL}
  # computation of Fold Change from log2FC
  if(parameters$FC==TRUE){colb="FC";ASKO_stat$FC <- 2^abs(ASKO_stat$logFC)}else{colb=NULL}
  if(parameters$Sign==TRUE){colc="Significance"}
  if(parameters$logCPM==TRUE){cold="logCPM"}else{cold=NULL}
  if(parameters$LR==TRUE){cole="LR"}else{cole=NULL}
  if(parameters$FDR==TRUE){colf="FDR"}else{colf=NULL}

  # adding table "stat.table" to the ASKOlist
  ASKOlist$stat.table<-ASKO_stat[,c("Test_id","contrast","gene",cola,colb,"PValue",colg,colc,cold,cole,colf)]

  if(parameters$mean_counts==TRUE){                            # computation of the mean of normalized counts for conditions
    mean1<-NormCountsMean(fit, ASKOlist, contx1)            # in the 1st context
    mean2<-NormCountsMean(fit, ASKOlist, contx2)            # in the 2nd context
    ASKOlist$stat.table<-cbind(ASKOlist$stat.table, mean1)
    ASKOlist$stat.table<-cbind(ASKOlist$stat.table, mean2)
  }
  print(table(ASKO_stat$Expression))
  colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="gene"] <- paste("is", "gene", sep="@")                  # header formatting for askomics
  colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="contrast"] <- paste("measured_in", "Contrast", sep="@") # header formatting for askomics
  o <- order(ASKOlist$stat.table$FDR)                                                                                   # ordering genes by FDR value
  ASKOlist$stat.table<-ASKOlist$stat.table[o,]

  utils::write.table(ASKOlist$stat.table,paste0(asko_dir, parameters$organism, "_", contrasko, ".txt"), sep=parameters$sep, col.names = TRUE, row.names = FALSE, quote=FALSE)

  # for image size
  nsamples<-ncol(dge$counts)
  sizeImg=15*nsamples
  if(sizeImg < 480) {sizeImg=480}

  # heatmap of Most Differential Genes Expression
  if(parameters$heatmap==TRUE){
    cpm_gstats<-edgeR::cpm(dge, log=TRUE)[o,][seq(parameters$numhigh),]
    grDevices::png(paste0(image_dir, contrast, "_topDGE_heatmap.png"), width=sizeImg*1.5, height=sizeImg*1.5)
    graphics::par(oma=c(2,2,2,2))
    gplots::heatmap.2(cpm_gstats,
              trace="none",
              scale="row",
              labCol=dge$samples$Name,
              main = contrasko,
              xlab = "samples",
              ColSideColors = dge$samples$color,
              margins = c(12,12),
              Rowv = FALSE,
              dendrogram="col")
    grDevices::dev.off()
  }
}
