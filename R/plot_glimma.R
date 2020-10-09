#' @title plot_glimma
#'
#' @description
#' Use Glimma package for interactive visualization of results from differential
#' expression analyses. Two types of graphs can be created:
#' \itemize{
#'    \item Mean-Difference Plot of Expression Data (aka MA plot).
#'    \item Volcano plot for a specified coefficient/contrast of a linear model.
#' }
#'
#' @param fit, fitted linear model object.
#' @param normGE, large DGEList with normalized counts and data description.
#' @param resDE, vector containing integer values of -1 to represent down-regulated
#' genes, 0 for no differential expression, and 1 for up-regulated genes.
#' @param contrast, coefficient/contrast tested.
#' @param tplot, type of plot selected for display.
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @return none.
#'
#' @examples
#' \dontrun{
#'    plot_glimma(fit, normGE, resDE, contrast, "MD", parameters)  # smear plot
#'    plot_glimma(fit, normGE, resDE, contrast, "VO", parameters)  # volcano plot
#' }
#'
#' @export
plot_glimma <- function(fit, normGE, resDE, contrast, tplot, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "images/")

  # Mean-Difference Plot
  #--------------------------
  if(tplot=="MD"){
    if (is.null(normGE$samples$color)==TRUE){
      suppressWarnings(Glimma::glMDPlot(fit, status=resDE[,contrast], counts=normGE, group=normGE$samples$condition,
                                transform=TRUE, anno=NULL, launch=FALSE, main=contrast,
                                folder=paste0(image_dir, "Glimma_Plots"), html=paste0("MDPlot_",contrast)))
    }
    else{
      suppressWarnings(Glimma::glMDPlot(fit, status=resDE[,contrast], counts=normGE, group=normGE$samples$condition, transform=TRUE,
                                sample.cols=normGE$samples$color, anno=NULL, launch=FALSE, main=contrast,
                                folder=paste0(image_dir, "Glimma_Plots"), html=paste0("MDPlot_",contrast)))
    }
  }

  # Volcano plot
  #--------------------------
  if (tplot=="VO"){
    if (is.null(normGE$samples$color)==TRUE){
      Glimma::glXYPlot(x=fit$table$logFC, y=-log10(fit$table$PValue), status=resDE[,contrast], counts=normGE,
               group=normGE$samples$condition, xlab="Log2FoldChange", ylab="-log10(pvalue)", main=contrast,
               launch=FALSE, folder=paste0(image_dir, "Glimma_Plots"), html=paste0("Volcano_",contrast))
    }
    else{
      Glimma::glXYPlot(x=fit$table$logFC, y=-log10(fit$table$PValue), status=resDE[,contrast], counts=normGE, main=contrast,
               group=normGE$samples$condition, xlab="Log2FoldChange", ylab="-log10(pvalue)", launch=FALSE,
               sample.cols=normGE$samples$color, folder=paste0(image_dir, "Glimma_Plots"), html=paste0("Volcano_",contrast))
    }
  }
}
