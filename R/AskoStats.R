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
#' @param glm_test tests for one or more coefficients in the linear model (likelihood ratio tests or empirical Bayes quasi-likelihood F-tests).
#' @param fit fitted linear model object.
#' @param contrast coefficient/contrast names tested.
#' @param ASKOlist list of data.frame contain condition, contrast and context informations made by asko3c.
#' @param dge large DGEList with normalized counts by GEnorm function.
#' @param data_list list contain all data and metadata (DGEList, samples descriptions, contrast, design and annotations).
#' @param parameters list that contains all arguments charged in Asko_start.
#' @return none
#'
#' @examples
#' \dontrun{
#'     AskoStats(glm_test, fit, contrast, ASKOlist, dge, data_list, parameters)
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
AskoStats <- function (glm_test, fit, contrast, ASKOlist, dge, data_list, parameters){
  study_dir = paste0(parameters$dir_path,"/",parameters$analysis_name,"/")
  asko_dir = paste0(study_dir, "DEanalysis/AskoTables/")
  image_dir = paste0(study_dir, "DEanalysis/DEimages/")
  table_dir = paste0(study_dir, "DEanalysis/DEtables/")
  norm_dir  = paste0(study_dir, "NormCountsTables/")

  contrasko<-ASKOlist$contrast$Contrast[row.names(ASKOlist$contrast)==contrast]   # to retrieve the name of contrast from Asko object
  contx1<-ASKOlist$contrast$context1[row.names(ASKOlist$contrast)==contrast]      # to retrieve the name of 1st context from Asko object
  contx2<-ASKOlist$contrast$context2[row.names(ASKOlist$contrast)==contrast]      # to retrieve the name of 2nd context from Asko object
  DETable=ASKOlist

  ASKO_stat<-glm_test$table
  ASKO_stat$contrast<-contrasko                                                   # addition of the contrast of the test
  ASKO_stat$gene<-row.names(ASKO_stat)                                            # addition of gene column = gene ID
  ASKO_stat$FDR<-stats::p.adjust(ASKO_stat$PValue, method=parameters$p_adj_method)       # computation of False Discovery Rate

  # Between context1 and context2 :
  ASKO_stat$Significance=0
  ASKO_stat$Significance[ASKO_stat$logFC <= -parameters$threshold_logFC & ASKO_stat$FDR <= parameters$threshold_FDR] = -1  # Significance values = -1 for down regulated genes
  ASKO_stat$Significance[ASKO_stat$logFC >= parameters$threshold_logFC  & ASKO_stat$FDR <= parameters$threshold_FDR] = 1   # Significance values =  1 for up regulated genes

  # addition of column "expression"
  ASKO_stat$Expression=NA
  ASKO_stat$Expression[ASKO_stat$Significance==-1]<-paste(contx1, contx2, sep="<")  # the value of attribute "Expression" is a string
  ASKO_stat$Expression[ASKO_stat$Significance==1]<-paste(contx1, contx2, sep=">")   # this attribute is easier to read the Significance
  ASKO_stat$Expression[ASKO_stat$Significance==0]<-paste(contx1, contx2, sep="=")   # of expression between two contexts

  if(parameters$Expression==TRUE){colg="Expression"}else{colg=NULL}
  if(parameters$logFC==TRUE){cola="logFC"}else{cola=NULL}
  if(parameters$FC==TRUE){colb="FC";ASKO_stat$FC <- 2^abs(ASKO_stat$logFC)}else{colb=NULL}
  if(parameters$Sign==TRUE){colc="Significance"}
  if(parameters$logCPM==TRUE){cold="logCPM"}else{cold=NULL}
  if(parameters$LR==TRUE){cole="LR"}else{cole=NULL}
  if(parameters$FDR==TRUE){colf="FDR"}else{colf=NULL}
  if(parameters$projectName!="DEprj" && stringr::str_replace_all(parameters$projectName, " ", "")!=""){
    colp="Project"
    ASKO_stat$Project<-parameters$projectName
    ASKO_stat$contrast2<-paste0(contx1, "vs", contx2)
    ASKO_stat$Test_id<-paste0(contx1, "vs", contx2, "_", rownames(ASKO_stat))
  }else{
    colp=NULL
    ASKO_stat$contrast2<-contrasko
    ASKO_stat$Test_id<-paste(contrasko, rownames(ASKO_stat), sep = "_")
  }

  # adding table "stat.table" to the ASKOlist and DETable
  ASKOlist$stat.table<-ASKO_stat[,c("Test_id",colp,"contrast","gene",cola,colb,"PValue",colg,colc,cold,cole,colf)]
  DETable$stat.table<-ASKO_stat[,c("gene","contrast2",colp,colg,colc,"PValue",cola,colb,colf,cold,cole)]

  grDevices::png(paste0(image_dir, contrast, "_Pval_Plot.png"),width=1000,height=500)
  graphics::par(mfrow = c(1,2))
  graphics::hist(ASKO_stat[,"PValue"], breaks=50, xlab = "p-value", main = "Raw p-values", xlim = c(0, 1))
  graphics::hist(ASKO_stat[,"FDR"], breaks=50, xlab = "p-value", main = "Adjusted p-values", xlim = c(0, 1))
  grDevices::dev.off()

  # Norm Mean Counts in DEtables
  #allElem<-rownames(which(data_list$contrast[contrast]!=0,arr.ind=TRUE)) #changé par la ligne suivante le 201220022 car allElem contenait les nom des conditions à la place des noms des contextes !!! Or il faut des noms de contexte pour que la fonction NormCountsMean fonctionne correctement
  allContext <- c(contx1,contx2)
  #for(contx in allElem){
  for(contx in allContext){
    mean<-NormCountsMean(fit, DETable, contx)
    #colnames(mean)<-paste0(contx,"_NormMeanCount") #remplacé car on veut bien ici les noms des conditions correspondantes aux conditions dans les contextes
    colnames(mean)<-paste0(ASKOlist$context$condition[ASKOlist$context$context==contx],"_NormMeanCount")
    DETable$stat.table<-cbind(DETable$stat.table, mean)
  }



  # CPM Norm Mean Count in DEtables
  moys<-utils::read.csv(paste0(norm_dir, parameters$analysis_name,"_CPM_NormMeanCounts.txt"), header=TRUE, sep="\t", row.names=1)
  rnames<-row.names(DETable$stat.table)                 # get Genes DE names
  allElem<-rownames(which(data_list$contrast[contrast]!=0,arr.ind=TRUE))
  CpmNMC<-as.matrix(moys[rnames,allElem])               # get cpm norm mean count for each genes DE
  rownames(CpmNMC)<-rnames
  colnames(CpmNMC)<-paste0(colnames(CpmNMC),"_CPMnormMeanCount")
  DETable$stat.table<-cbind(DETable$stat.table, CpmNMC) # merge the two matrix

  # Annotation Genes in DEtables
  if(is.null(parameters$annotation)==FALSE)
  {
    rnames<-row.names(DETable$stat.table)                # get Genes DE names
    annDE<-as.matrix(data_list$annot[rnames,])           # get annotations for each genes DE
    rownames(annDE)<-rnames
    colnames(annDE)<-colnames(data_list$annot)
    DETable$stat.table<-cbind(DETable$stat.table, annDE) # merge the two matrix
  }

  print(table(ASKO_stat$Expression))
  colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="gene"] <- paste("is", "gene", sep="@")                  # header formatting for askomics
  if(parameters$projectName!="DEprj" && stringr::str_replace_all(parameters$projectName, " ", "")!=""){
    colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="Project"] <- paste("from", "Project", sep="@")        # header formatting for askomics
  }
  colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="contrast"] <- paste("measured_in", "Contrast", sep="@") # header formatting for askomics
  o <- order(ASKOlist$stat.table$FDR)                                                                                   # ordering genes by FDR value
  ASKOlist$stat.table<-ASKOlist$stat.table[o,]

  o <- order(DETable$stat.table$FDR)                                                                                   # ordering genes by FDR value
  DETable$stat.table<-DETable$stat.table[o,]
  colnames(DETable$stat.table)[colnames(DETable$stat.table)=="contrast2"] <- "contrast"

  utils::write.table(ASKOlist$stat.table,paste0(asko_dir, "Asko_", contrasko, ".txt"), sep=parameters$sep, col.names = TRUE, row.names = FALSE, quote=FALSE)
  utils::write.table(DETable$stat.table,paste0(table_dir, contrasko, ".txt"), sep=parameters$sep, col.names = TRUE, row.names = FALSE, quote=FALSE)

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
