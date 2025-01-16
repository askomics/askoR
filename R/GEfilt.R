#' @title GEfilt
#'
#' @description
#' \itemize{
#'    \item Filter genes according to cpm threshold value and replicated cpm option value.
#'    \item Plot different graphes to explore data before and after filtering.
#' }
#'
#' @param data_list list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations)
#' @param parameters list that contains all arguments charged in Asko_start
#' @return large DGEList with filtered counts and data descriptions.
#'
#' @examples
#' \dontrun{
#'     filtered_counts<-GEfilt(data, parameters)
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
GEfilt <- function(data_list, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "DataExplore/")

  # plot density before filtering
  cpm<-edgeR::cpm(data_list$dge)
  logcpm<-edgeR::cpm(data_list$dge, log=TRUE)
  colnames(logcpm)<-rownames(data_list$dge$samples)
  nsamples <- ncol(data_list$dge$counts)

  maxi<-c()
  for (i in seq(nsamples)){
    m=max(stats::density(logcpm[,i])$y)
    maxi<-c(maxi,m)
  }
  ymax<-round(max(maxi),1) + 0.02

  sizeImg=15*nsamples
  if(sizeImg < 480){ sizeImg=480 }
  btm=round((nsamples/6),0)+0.5

  grDevices::png(paste0(image_dir, parameters$analysis_name, "_densityGraph_RawData.png"), width=1024, height=1024)
  graphics::par(oma=c(2,2,2,0), mar=c(parameters$densbotmar,5,5,5))
  plot(stats::density(logcpm[,1]),
       col=as.character(data_list$dge$samples$color[1]),
       lwd=1,
       las=2,
       ylim=c(0,ymax),
       main="A. Raw data",
       xlab="Log-cpm")
  graphics::abline(v=0, lty=3)
  for (i in 2:nsamples){
    den<-stats::density(logcpm[,i])
    graphics::lines(den$x, col=as.character(data_list$dge$samples$color[i]), den$y, lwd=1)
  }

  graphics::legend("bottom", fill=data_list$dge$samples$color, bty="n", ncol=parameters$legendcol,
                   legend=rownames(data_list$dge$samples), xpd=TRUE, inset=-parameters$densinset)
  grDevices::dev.off()

  # plot density after filtering
  keep.exprs <- rowSums(cpm>parameters$threshold_cpm)>=parameters$replicate_cpm
  filtered_counts <- data_list$dge[keep.exprs,,keep.lib.sizes=FALSE]
  filtered_cpm<-edgeR::cpm(filtered_counts$counts, log=TRUE)

  maxi<-c()
  for (i in seq(nsamples)){
    m=max(stats::density(filtered_cpm[,i])$y)
    maxi<-c(maxi,m)
  }
  ymax<-round(max(maxi),1) + 0.02

  grDevices::png(paste0(image_dir,parameters$analysis_name,"_densityGraph_FilteredData.png"), width=1024, height=1024)
  graphics::par(oma=c(2,2,2,0), mar=c(parameters$densbotmar,5,5,5))
  plot(stats::density(filtered_cpm[,1]),
       col=as.character(data_list$dge$samples$color[1]),
       lwd=1,
       ylim=c(0,ymax),
       las=2,
       main="B. Filtered data",
       xlab="Log-cpm")
  graphics::abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- stats::density(filtered_cpm[,i])
    graphics::lines(den$x,col=as.character(data_list$dge$samples$color[i]), den$y, lwd=1)
  }
  graphics::legend("bottom", fill=data_list$dge$samples$color, bty="n", ncol=parameters$legendcol,
                   legend=rownames(data_list$dge$samples), xpd=TRUE, inset=-parameters$densinset)
  grDevices::dev.off()

  # histogram cpm values distribution before filtering
  grDevices::png(paste0(image_dir,parameters$analysis_name,"_barplot_logcpm_before_filtering.png"), width=sizeImg, height=sizeImg)
  graphics::hist(logcpm,
                 main= "A. Log2(cpm) distribution before filtering",
                 xlab = "log2(cpm)",
                 col = "grey")
  grDevices::dev.off()

  # histogram cpm values distribution after filtering
  grDevices::png(paste0(image_dir,parameters$analysis_name,"_barplot_logcpm_after_filtering.png"), width=sizeImg, height=sizeImg)
  graphics::hist(filtered_cpm,
                 main= "B. Log2(cpm) distribution after filtering",
                 xlab = "log2(cpm)",
                 col = "grey")
  grDevices::dev.off()

  # boxplot cpm values distribution before filtering
  grDevices::png(paste0(image_dir,parameters$analysis_name,"_boxplot_logcpm_before_filtering.png"), width=sizeImg, height=sizeImg)
  graphics::par(oma=c(1,1,1,1))
  graphics::boxplot(logcpm,
                    col=data_list$dge$samples$color,
                    main="A. Log2(cpm) distribution before filtering",
                    cex.axis=0.8,
                    las=2,
                    ylab="log2(cpm)")
  grDevices::dev.off()

  # boxplot cpm values distribution after filtering
  grDevices::png(paste0(image_dir,parameters$analysis_name,"_boxplot_logcpm_after_filtering.png"), width=sizeImg, height=sizeImg)
  graphics::par(oma=c(1,1,1,1))
  graphics::boxplot(filtered_cpm,
                    col=data_list$dge$samples$color,
                    main="B. Log2(cpm) distribution after filtering",
                    cex.axis=0.8,
                    las=2,
                    ylab="log2(cpm)")
  grDevices::dev.off()


  # barplot count distribution before filtering
  grDevices::png(paste0(image_dir,parameters$analysis_name,"_barplot_SumCounts_before_filtering.png"), width=sizeImg, height=sizeImg)
  graphics::par(oma=c(1,1,1,1),mar=c(7, 7, 4, 2), mgp=c(5,0.7,0))
  graphics::barplot(colSums(data_list$dge$counts),
                    col=data_list$dge$samples$color,
                    main=paste0("A. Sum of raw counts from ",nrow(data_list$dge$counts)," transcripts"),
                    cex.axis=0.8,
                    las=2,
                    ylab="Counts sum",
                    xlab="Samples")
  grDevices::dev.off()

  # barplot count distribution after filtering
  grDevices::png(paste0(image_dir,parameters$analysis_name,"_barplot_SumCounts_after_filtering.png"), width=sizeImg, height=sizeImg)
  graphics::par(oma=c(1,1,1,1),mar=c(7, 7, 4, 2), mgp=c(5,0.7,0))
  graphics::barplot(colSums(filtered_counts$counts),
                    col=data_list$dge$samples$color,
                    main=paste0("B. Sum of filtered counts from ",nrow(filtered_counts$counts)," transcripts"),
                    cex.axis=0.8,
                    las=2,
                    ylab="Counts sum",
                    xlab="Samples")
  grDevices::dev.off()

  return(filtered_counts)
}
