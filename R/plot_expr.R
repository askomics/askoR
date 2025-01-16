#' @title plot_expr
#'
#' @description
#' Function to generate plots showing different aspects of differential expression results.
#' Two types of graphs can be created:
#' \itemize{
#'    \item Mean-Difference Plot of Expression Data (aka MA plot).
#'    \item Volcano plot for a specified coefficient/contrast of a linear model.
#' }
#'
#' @param fit fitted linear model object.
#' @param normGE large DGEList with normalized counts and data description.
#' @param resDE vector containing integer values of -1 to represent down-regulated
#' genes, 0 for no differential expression, and 1 for up-regulated genes.
#' @param contrast coefficient/contrast tested.
#' @param tplot type of plot selected for display.
#' @param parameters list that contains all arguments charged in Asko_start.
#' @return none.
#'
#' @examples
#' \dontrun{
#'    plot_expr(fit, normGE, resDE, contrast, "MD", parameters)  # smear plot
#'    plot_expr(fit, normGE, resDE, contrast, "VO", parameters)  # volcano plot
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
plot_expr <- function(fit, normGE, resDE, contrast, tplot, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "DEanalysis/DEimages/")



  # Mean-Difference Plot
  if(tplot=="MD"){
    grDevices::png(paste0(image_dir, contrast, "_MeanDifference_of_ExpressionData.png"))
    edgeR::plotMD.DGELRT(fit, xlab="Average log CPM", ylab="log-fold-change", main=paste0("MD plot - ", contrast),
                         cex=0.5, status=resDE[,contrast], values=c("-1","1"), col=c("blue","red"))
    grDevices::dev.off()
  }

  # Volcano plot
  if(tplot=="VO"){
    tglm<-fit$table
    tglm$FDR<-stats::p.adjust(tglm$PValue, method=parameters$p_adj_method)
    tglm$PValue[tglm$PValue==0]<-1e-300

    grDevices::png(paste0(image_dir, contrast, "_VolcanoPlot.png"))
    with(tglm, plot(tglm$logFC, -log10(tglm$PValue), pch=16, cex=0.5, xlim=c(min(tglm$logFC)-0.5, max(tglm$logFC)+0.5),
                    ylim=c(min(-log10(tglm$PValue))-0.5, max(-log10(tglm$PValue))+0.5),
                    main=paste0("Volcano plot - ", contrast), xlab="Log2FoldChange", ylab="-log10(pvalue)"))
    subset1 <- subset(tglm, tglm$FDR <= parameters$threshold_FDR & tglm$logFC >  parameters$threshold_logFC)
    with(subset1, graphics::points(subset1$logFC, -log10(subset1$PValue), pch=16, cex=0.5, col="red"))
    subset2 <- subset(tglm, tglm$FDR <= parameters$threshold_FDR & tglm$logFC < -parameters$threshold_logFC)
    with(subset2, graphics::points(subset2$logFC, -log10(subset2$PValue), pch=16, cex=0.5, col="blue"))
    grDevices::dev.off()
  }
}
