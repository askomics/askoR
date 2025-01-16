#' @title DEanalysis
#'
#' @description Genewise statistical tests for a given coefficient or contrast, with edgeR method.
#'
#' @param norm_GE large DGEList with normalized counts and data description.
#' @param data_list list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations).
#' @param asko_list list of data.frame contain condition, contrast and context informations made by asko3c.
#' @param parameters list that contains all arguments charged in Asko_start.
#' @return list (TestResults format class limma) contains for each contrast the significance expression (1/0/-1) for all gene.
#'
#' @import edgeR
#' @import limma
#'
#' @examples
#' \dontrun{
#'     sum_table<-DEanalysis(norm_GE, data_list, asko_list, parameters)
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
DEanalysis <- function(norm_GE, data_list, asko_list, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "DEanalysis/DEimages/")
  table_dir = paste0(study_dir, "DEanalysis/DEtables/")

  # for image size
  nsamples <- ncol(data_list$dge$counts)
  sizeImg=15*nsamples
  if(sizeImg < 480){ sizeImg=480 }

  # Checks Contrasts
  if(is.null(parameters$select_sample) & is.null(parameters$rm_sample)){
    c1<-unique(data_list$samples$condition)
    len1<-length(c1)
    c2<-rownames(data_list$contrast)
    len2<-length(c2)
    if(len1 > len2){
      cat("\n\n")
      stop("Contrast files must be contain all conditions (in rows).\n\n")
    }
    if(len1 < len2){
      cat("\n\n")
      stop("Too much condtions in contrast file!\n\n")
    }
    if(length(setdiff(c1,c2)) > 0 && length(setdiff(c2,c1)) > 0){
      cat("\n\n")
      stop("Erronate or unknown conditions names in contrast file!\n\n")
    }
  }

  # Estimate Common, Trended and Tagwise Dispersion for Negative Binomial GLMs
  if(parameters$glm_disp==TRUE)
  {
    normGEdisp <- estimateGLMCommonDisp(norm_GE, data_list$design)
    normGEdisp <- estimateGLMTrendedDisp(normGEdisp, data_list$design)
    normGEdisp <- estimateGLMTagwiseDisp(normGEdisp, data_list$design)
  }
  # Estimate Common, Trended and Tagwise Negative Binomial dispersions by weighted likelihood empirical Bayes
  else
  {
    normGEdisp <- estimateDisp(norm_GE, data_list$design)
  }
  grDevices::png(paste0(image_dir, "Biological_coefficient_of_variation.png"), width=sizeImg, height=sizeImg)
  plotBCV(normGEdisp)
  grDevices::dev.off()

  # Genewise Negative Binomial Generalized Linear Models
  if(parameters$glm=="lrt"){
    fit <- glmFit(normGEdisp, data_list$design, robust = TRUE)
  }
  # Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests
  else if(parameters$glm=="qlf"){
    fit <- glmQLFit(normGEdisp, data_list$design, robust = TRUE)
    grDevices::png(paste0(image_dir, parameters$analysis_name, "_quasi-likelihood_dispersion.png"), width=sizeImg, height=sizeImg)
    plotQLDisp(fit)
    grDevices::dev.off()
  }



  # data frame combine all status genes results for summary file
  sum<-data.frame(row.names = rownames(fit))
  sum2=list()
  # if only one contrast ask
  if(length(data_list$contrast)==1){

    contrast<-makeContrasts(contrasts = data_list$contrast, levels = data_list$design)
    colnames(contrast)<-colnames(data_list$contrast)

    # likelihood ratio tests for one or more coefficients in the linear model.
    if(parameters$glm=="lrt" & parameters$threshold_logFC==0){
      glm_test<-glmLRT(fit, contrast=contrast)
    }
    # similar to glmLRT except that it replaces likelihood ratio tests with empirical Bayes quasi-likelihood F-tests
    if(parameters$glm=="qlf" & parameters$threshold_logFC==0){
      glm_test<-glmQLFTest(fit, contrast=contrast)
    }
    if(parameters$threshold_logFC!=0){
      glm_test<-glmTreat(fit, contrast=contrast, lfc=parameters$threshold_logFC)
    }

    sum[,colnames(contrast)]<-decideTests(glm_test, adjust.method = parameters$p_adj_method, lfc=parameters$threshold_logFC, p.value=parameters$threshold_FDR)
    AskoStats(glm_test, fit, colnames(contrast), asko_list, normGEdisp, data_list, parameters)

    # display grahes (volcano or/and MD)
    if(parameters$plotMD==TRUE) { plot_expr(glm_test, normGEdisp, sum, colnames(contrast), "MD", parameters) }
    if(parameters$plotVO==TRUE) { plot_expr(glm_test, normGEdisp, sum, colnames(contrast), "VO", parameters) }
    if(parameters$glimMD==TRUE) { plot_glimma(glm_test, normGEdisp, sum, colnames(contrast), "MD", parameters) }
    if(parameters$glimVO==TRUE) { plot_glimma(glm_test, normGEdisp, sum, colnames(contrast), "VO", parameters) }
  }
  # for more than one contrast
  else{
    for (contrast in colnames(data_list$contrast)){

      if(parameters$glm=="lrt" & parameters$threshold_logFC==0){
        glm_test<-glmLRT(fit, contrast=data_list$contrast[,contrast])
      }
      # similar to glmLRT except that it replaces likelihood ratio tests with empirical Bayes quasi-likelihood F-tests
      if(parameters$glm=="qlf" & parameters$threshold_logFC==0){
        glm_test<-glmQLFTest(fit, contrast=data_list$contrast[,contrast])
      }
      if(parameters$threshold_logFC!=0){
        glm_test<-glmTreat(fit, contrast=data_list$contrast[,contrast], lfc=parameters$threshold_logFC)
      }
      sum[,contrast]<-decideTests(glm_test, adjust.method = parameters$p_adj_method, p.value=parameters$threshold_FDR)
      AskoStats(glm_test, fit, contrast, asko_list, normGEdisp, data_list, parameters)

      # display grahes (volcano or/and MD)
      if(parameters$plotMD==TRUE) { plot_expr(glm_test, normGEdisp, sum, contrast, "MD", parameters) }
      if(parameters$plotVO==TRUE) { plot_expr(glm_test, normGEdisp, sum, contrast, "VO", parameters) }
      if(parameters$glimMD==TRUE) { plot_glimma(glm_test, normGEdisp, sum, contrast, "MD", parameters) }
      if(parameters$glimVO==TRUE) { plot_glimma(glm_test, normGEdisp, sum, contrast, "VO", parameters) }
    }
  }

  # Create summary file with annotations (if available) and contrast value for each gene
  cat("\nCreate Summary file\n\n")
  sumFile<-paste0(table_dir,"Summary_DEresults.txt")
  if(is.null(data_list$annot)==FALSE)
  {
    rnames<-row.names(sum)                        # get Genes DE names
    annDE<-as.matrix(data_list$annot[rnames,])    # get annotations for each genes DE
    rownames(annDE)<-rnames
    colnames(annDE)<-colnames(data_list$annot)
    SumMat<-cbind(sum,annDE)                      # merge the two matrix

    utils::write.table(SumMat, file=sumFile, col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')
  }
  else
  {
    utils::write.table(sum, file=sumFile, col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')
  }

  # DE graph in contrasts
  tableau = data.frame(contrast=colnames(data_list$contrast), UP=0, DOWN=0)

  for (i in 1:ncol(data_list$contrast)) {
    tableau[i,2]<-sum(sum[,tableau$contrast[i]]==1)
    tableau[i,3]<-sum(sum[,tableau$contrast[i]]==-1)
  }

  tableau_final = tidyr::pivot_longer(tableau, cols = c("UP","DOWN"))
  comp_names <- c( `UP` = "UP in first condition", `DOWN` = "DOWN in first condition")

  ggplot2::ggplot(data=tableau_final, aes(x=tableau_final$contrast, y=tableau_final$value, fill=tableau_final$name)) +
    geom_bar(stat="identity")+
    geom_text(aes(label=tableau_final$value), vjust=-0.5,
              color="black", size=3.5)+
    scale_fill_brewer(palette="Paired", labels=comp_names)+
    ylim(0, max(tableau_final$value)+(0.1*max(tableau_final$value)))+
    #theme_minimal()+
    facet_wrap(~fct_rev(name), ncol=1)+
    labs(fill = "Differential\nExpression",title = "Number of DEGs in each contrast", subtitle="(DEGs = Diffentially Expressed Genes)")+
    ylab("Number of DEGs")+
    xlab("Contrast")
  ggplot2::ggsave(filename=paste0(image_dir,"DEGsnumber_InContrasts.png"),width=10+(0.5*ncol(data_list$contrast)), height=10)


  # reformate summary result table
  newMat <- as.data.frame(matrix(unlist(sum), nrow=nrow(sum)))
  rownames(newMat)<-rownames(sum)
  colnames(newMat)<-colnames(sum)

  return(newMat)
}
