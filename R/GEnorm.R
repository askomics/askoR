#' @title GEnorm
#'
#' @description Normalize counts
#' \itemize{
#'    \item Calculate normalization factors to scale the filtered library sizes.
#'    \item Plot different graphes to explore data before and after normalization.
#'    \item Optionally, write file with mean counts and normalized mean counts in Askomics format.
#' }
#'
#' @param filtered_GE large DGEList with filtered counts by GEfilt function.
#' @param parameters list that contains all arguments charged in Asko_start.
#' @param asko_list list of data.frame contain condition, contrast and context informations made by asko3c.
#' @param data_list list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations)
#' @return large DGEList with normalized counts and data descriptions.
#'
#' @examples
#' \dontrun{
#'     norm_GE<-GEnorm(filtered_counts, asko, parameters)
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
GEnorm <- function(filtered_GE, asko_list, data_list, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "DataExplore/")
  norm_dir  = paste0(study_dir, "NormCountsTables/")

  # for image size
  nsamples <- ncol(filtered_GE$counts)
  sizeImg=15*nsamples
  if(sizeImg < 480){ sizeImg=480 }

  # Normalization counts
  norm_GE<-edgeR::calcNormFactors(filtered_GE, method = parameters$normal_method)
  if(parameters$norm_factor == TRUE){
    utils::write.table(norm_GE$samples, file=paste0(norm_dir, parameters$analysis_name, "_normalization_factors.txt"), col.names=NA, row.names=TRUE, quote=FALSE, sep="\t", dec=".")
  }

  # boxplot log2(cpm) values after normalization
  logcpm_norm <- edgeR::cpm(norm_GE, log=TRUE)
  colnames(logcpm_norm)<-rownames(filtered_GE$samples)

  grDevices::png(paste0(image_dir,parameters$analysis_name,"_boxplot_logcpm_after_norm.png"), width=sizeImg, height=sizeImg)
  graphics::par(oma=c(1,1,1,1))
  graphics::boxplot(logcpm_norm,
                    col=filtered_GE$samples$color,
                    main="C. Log2(cpm) distribution after normalization",
                    cex.axis=0.8,
                    las=2,
                    ylab="Log2(cpm)")
  grDevices::dev.off()

  dScaleFactors <- norm_GE$samples$lib.size * norm_GE$samples$norm.factors
  normCounts <- t(t(filtered_GE$counts)/dScaleFactors)*mean(dScaleFactors)

  # File with normalized counts
  if (parameters$norm_counts == TRUE){
    utils::write.table(normCounts, file=paste0(norm_dir, parameters$analysis_name, "_NormCounts.txt"), col.names=NA, row.names=TRUE, quote=FALSE, sep="\t", dec=".")
  }


  # barplot counts after normalization
  grDevices::png(paste0(image_dir,parameters$analysis_name,"_barplot_SumCounts_after_norm.png"), width=sizeImg, height=sizeImg)
  graphics::par(oma=c(1,1,1,1),mar=c(7, 7, 4, 2), mgp=c(5,0.7,0))
  graphics::barplot(colSums(normCounts),
                    col=data_list$dge$samples$color,
                    main=paste0("C. Sum of normalized counts from ",nrow(normCounts)," transcripts"),
                    cex.axis=0.8,
                    las=2,
                    ylab="Counts sum",
                    xlab="Samples")
  grDevices::dev.off()

  # heatmap visualisation
  if(parameters$CompleteHeatmap==TRUE)
  {
    # heatmap cpm value per sample
    cpm_norm  <- edgeR::cpm(norm_GE, log=FALSE)
    cpmscale  <- scale(t(cpm_norm))
    tcpmscale <- t(cpmscale)

    d1 <- stats::dist(cpmscale,  method = parameters$distcluts, diag = FALSE, upper = FALSE)
    d2 <- stats::dist(tcpmscale, method = parameters$distcluts, diag = FALSE, upper = TRUE)
    hc <- stats::hclust(d1, method = parameters$hclust, members = NULL)
    hr <- stats::hclust(d2, method = parameters$hclust, members = NULL)
    my_palette <- grDevices::colorRampPalette(c("green","black","red"), interpolate = "linear")

    grDevices::png(paste0(image_dir,parameters$analysis_name,"_heatmap_CPMcounts_per_sample.png"), width=sizeImg*1.5, height=sizeImg*1.25)
    graphics::par(oma=c(2,1,2,2))
    gplots::heatmap.2(tcpmscale, Colv = stats::as.dendrogram(hc), Rowv = stats::as.dendrogram(hr), density.info="histogram",
                      trace = "none", dendrogram = "column", xlab = "samples", col = my_palette, labRow = FALSE,
                      cexRow = 0.1, cexCol = 1.25, ColSideColors = norm_GE$samples$color, margins = c(10,1),
                      main = paste0("CPM counts per sample\nGenes 1 to ",nrow(norm_GE)))
    grDevices::dev.off()

    # Normalized mean by conditions
    # heatmap mean counts per condition
    n_count <- NormCountsMean(norm_GE, ASKOlist = asko_list)
    countscale  <- scale(t(n_count))
    tcountscale <- t(countscale)

    d1 <- stats::dist(countscale,  method = parameters$distcluts, diag = FALSE, upper = FALSE)
    d2 <- stats::dist(tcountscale, method = parameters$distcluts, diag = FALSE, upper = TRUE)
    hc <- stats::hclust(d1, method = parameters$hclust, members = NULL)
    hr <- stats::hclust(d2, method = parameters$hclust, members = NULL)
    my_palette <- grDevices::colorRampPalette(c("green","black","red"), interpolate = "linear")

    grDevices::png(paste0(image_dir,parameters$analysis_name,"_heatmap_CPMmean_per_condition.png"), width=sizeImg*1.5, height=sizeImg*1.25)
    graphics::par(oma=c(2,1,2,2))
    gplots::heatmap.2(tcountscale, Colv = stats::as.dendrogram(hc), Rowv = stats::as.dendrogram(hr), density.info="histogram",
                      trace = "none", dendrogram = "column", xlab = "Condition", col = my_palette, labRow = FALSE,
                      cexRow = 0.1, cexCol = 1.5, ColSideColors = unique(norm_GE$samples$color), margins = c(10,1),
                      main = paste0("CPM counts per condition (mean)\nGenes 1 to ",nrow(norm_GE)))
    grDevices::dev.off()
  }

  # File with normalized counts in CPM by sample
  cpm_norm <- edgeR::cpm(norm_GE, log=FALSE)
  utils::write.table(cpm_norm, file=paste0(norm_dir, parameters$analysis_name, "_CPM_NormCounts.txt"), col.names=NA, row.names=TRUE, quote=FALSE, sep="\t", dec=".")
  # File with normalized mean counts in CPM, grouped by condition
  tempo<-as.data.frame(t(stats::aggregate(t(cpm_norm),list(data_list$samples$condition), mean)))
  colnames(tempo)<-tempo["Group.1",]
  meancpmDEGnorm<-tempo[-1,]
  utils::write.table(meancpmDEGnorm,paste0(norm_dir, parameters$analysis_name,"_CPM_NormMeanCounts.txt"), sep="\t", dec=".", row.names=TRUE, col.names=NA, quote=FALSE)

  return(norm_GE)
}
