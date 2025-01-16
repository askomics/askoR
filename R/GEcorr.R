#' @title GEcorr
#'
#' @description
#' Plot some graphes to see data correlation :
#' \itemize{
#'    \item heatmap sample correslation
#'    \item MDS plots
#'    \item hierarchical clustering
#' }
#'
#' @param asko_norm large DGEList with normalized counts by GEnorm function.
#' @param parameters list that contains all arguments charged in Asko_start.
#' @return none
#'
#' @import ggfortify
#'
#' @examples
#' \dontrun{
#'     GEcorr(asko_norm,parameters)
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
GEcorr <- function(asko_norm, parameters){
  options(warn = -1)
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "DataExplore/")

  # for image size
  nsamples <- ncol(asko_norm$counts)
  sizeImg=15*nsamples
  if(sizeImg < 1024){ sizeImg=1024 }

  lcpm<-edgeR::cpm(asko_norm, log=TRUE)
  colnames(lcpm)<-rownames(asko_norm$samples)

  # Heatmap sample correlation
  cormat<-stats::cor(lcpm)
  color<-grDevices::colorRampPalette(c("black","red","yellow","white"),space="rgb")(35)
  grDevices::png(paste0(image_dir, parameters$analysis_name, "_heatmap_of_sample_correlation.png"), width=sizeImg, height=sizeImg)
  graphics::par(oma=c(4,2,4,1))
  stats::heatmap(cormat, col=color, symm=TRUE, RowSideColors=as.character(asko_norm$samples$color),
                 ColSideColors=as.character(asko_norm$samples$color), main="")
  graphics::title("Sample Correlation Matrix", adj=0.5, outer=TRUE)
  grDevices::dev.off()

  # Correlogram
  corr<-stats::cor(edgeR::cpm(asko_norm, log=FALSE))
  minCorr=min(corr)
  maxCorr=max(corr)
  grDevices::png(paste0(image_dir, parameters$analysis_name, "_correlogram.png"), width=sizeImg, height=sizeImg)
  corrplot::corrplot(corr, method="ellipse", type = "lower", tl.col = "black", tl.srt = 45, is.corr = FALSE, cl.lim=c(minCorr,maxCorr),col=RColorBrewer::brewer.pal(n=8, name="RdBu"))
  graphics::title("Sample Correlogram", adj=0.5)
  grDevices::dev.off()

  # MDS Plot
  mds <- stats::cmdscale(stats::dist(t(lcpm)),k=3, eig=TRUE)
  eigs<-round((mds$eig)*100/sum(mds$eig[mds$eig>0]),2)
  dfmds<-as.data.frame(mds$points)

  # Axe 1 and 2
  grDevices::png(paste0(image_dir, parameters$analysis_name, "_MDS_corr_axe1_2.png"), width=sizeImg*1.25, height=sizeImg*1.25)
  mds1<-ggplot2::ggplot(dfmds, ggplot2::aes(dfmds$V1, dfmds$V2, label=rownames(mds$points))) + ggplot2::labs(title="MDS Axes 1 and 2") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::theme(plot.margin=ggplot2::margin(20,30,20,15)) +
    ggplot2::geom_point(ggplot2::aes(color=as.character(asko_norm$samples$condition)),shape=17,size=15 ) +
    ggplot2::xlab(paste('dim 1 [', eigs[1], '%]')) + ggplot2::ylab(paste('dim 2 [', eigs[2], "%]")) +
    ggrepel::geom_label_repel(ggplot2::aes(label = rownames(mds$points),fill = factor(as.character(asko_norm$samples$condition))), color = 'white',size = 7) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face="bold",size=30),
      axis.text.x = ggplot2::element_text(face="bold",size=30),
      legend.text = ggplot2::element_text(size = 30),
      legend.title = ggplot2::element_text(size=30,face="bold"),
      strip.text.y = ggplot2::element_text(size=30, face="bold"),
      axis.title = ggplot2::element_text(size=30,face="bold"),
      plot.title = ggplot2::element_text(size=35)) +
    ggplot2::guides(color=FALSE, fill = ggplot2::guide_legend(title = "Condition", override.aes = ggplot2::aes(label = ""), ncol=parameters$legendcol)) +
    ggplot2::theme(legend.position="bottom",legend.direction="vertical",legend.margin=ggplot2::margin(5,5,5,5),legend.box.margin=ggplot2::margin(10,10,10,10))
  print(mds1)
  grDevices::dev.off()

  # Axe 2 and 3
  grDevices::png(paste0(image_dir, parameters$analysis_name, "_MDS_corr_axe2_3.png"), width=sizeImg*1.25, height=sizeImg*1.25)
  mds2<-ggplot2::ggplot(dfmds, ggplot2::aes(dfmds$V2, dfmds$V3, label = rownames(mds$points))) + ggplot2::labs(title="MDS Axes 2 and 3") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::theme(plot.margin=ggplot2::margin(20,30,20,15)) +
    ggplot2::geom_point(ggplot2::aes(color=as.character(asko_norm$samples$condition)),shape=17,size=15 ) +
    ggplot2::xlab(paste('dim 2 [', eigs[2], '%]')) + ggplot2::ylab(paste('dim 3 [', eigs[3], "%]")) +
    ggrepel::geom_label_repel(ggplot2::aes(label = rownames(mds$points),fill = factor(as.character(asko_norm$samples$condition))), color = 'white',size = 7) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face="bold",size=30),
      axis.text.x = ggplot2::element_text(face="bold",size=30),
      legend.text = ggplot2::element_text(size = 30),
      legend.title = ggplot2::element_text(size=30,face="bold"),
      strip.text.y = ggplot2::element_text(size=30, face="bold"),
      axis.title=ggplot2::element_text(size=30,face="bold"),
      plot.title = ggplot2::element_text(size=35)) +
    ggplot2::guides(color=FALSE, fill = ggplot2::guide_legend(title = "Condition", override.aes = ggplot2::aes(label = ""), ncol=parameters$legendcol)) +
    ggplot2::theme(legend.position="bottom",legend.direction="vertical",legend.margin=ggplot2::margin(5,5,5,5),legend.box.margin=ggplot2::margin(10,10,10,10))
  print(mds2)
  grDevices::dev.off()

  # Axe 1 and 3
  grDevices::png(paste0(image_dir, parameters$analysis_name, "_MDS_corr_axe1_3.png"), width=sizeImg*1.25, height=sizeImg*1.25)
  mds3<-ggplot2::ggplot(dfmds, ggplot2::aes(dfmds$V1, dfmds$V3, label = rownames(mds$points))) + ggplot2::labs(title="MDS Axes 1 and 3") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::theme(plot.margin=ggplot2::margin(20,30,20,15)) +
    ggplot2::geom_point(ggplot2::aes(color=as.character(asko_norm$samples$condition)),shape=17,size=15 ) +
    ggplot2::xlab(paste('dim 1 [', eigs[1], '%]')) + ggplot2::ylab(paste('dim 3 [', eigs[3], "%]")) +
    ggrepel::geom_label_repel(ggplot2::aes(label = rownames(mds$points),fill = factor(as.character(asko_norm$samples$condition))), color = 'white',size = 7) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face="bold",size=30),
      axis.text.x = ggplot2::element_text(face="bold",size=30),
      legend.text = ggplot2::element_text(size = 30),
      legend.title = ggplot2::element_text(size=30,face="bold"),
      strip.text.y = ggplot2::element_text(size=30, face="bold"),
      axis.title=ggplot2::element_text(size=30,face="bold"),
      plot.title = ggplot2::element_text(size=35)) +
    ggplot2::guides(color=FALSE, fill = ggplot2::guide_legend(title = "Condition",override.aes = ggplot2::aes(label = ""), ncol=parameters$legendcol)) +
    ggplot2::theme(legend.position="bottom",legend.direction="vertical",legend.margin=ggplot2::margin(5,5,5,5),legend.box.margin=ggplot2::margin(10,10,10,10))
  print(mds3)
  grDevices::dev.off()

  # PCA
  lcpm2=as.data.frame(t(lcpm))
  lcpm3=lcpm2
  lcpm3$cond=asko_norm$samples$condition

  # Axe 1 and 2
  grDevices::png(paste0(image_dir, parameters$analysis_name, "_PCA_axe1_2.png"), width=sizeImg*1.25, height=sizeImg*1.25)
  p=ggplot2::autoplot(stats::prcomp(lcpm2), data = lcpm3, size=15,colour="cond",show.legend = FALSE, x=1, y=2,shape=17)
  test=p + ggrepel::geom_label_repel(ggplot2::aes(label = rownames(lcpm3),fill = factor(lcpm3$cond)), color = 'white',size = 7) +
    ggplot2::labs(title="PCA Axes 1 and 2") + ggplot2::theme(
      plot.margin=ggplot2::margin(20,30,20,15),
      axis.text.y = ggplot2::element_text(face="bold",size=30),
      axis.text.x = ggplot2::element_text(face="bold",size=30),
      legend.text = ggplot2::element_text(size = 30),
      legend.title = ggplot2::element_text(size=30,face="bold"),
      strip.text.y = ggplot2::element_text(size=30, face="bold"),
      axis.title=ggplot2::element_text(size=30,face="bold"),
      plot.title = ggplot2::element_text(size=35, hjust=0.5)) +
    ggplot2::guides(colour=FALSE, fill = ggplot2::guide_legend(title = "Condition",override.aes = ggplot2::aes(label = ""), ncol=parameters$legendcol)) +
    ggplot2::theme(legend.position="bottom",legend.direction="vertical",legend.margin=ggplot2::margin(5,5,5,5),legend.box.margin=ggplot2::margin(10,10,10,10))
  print(test)
  grDevices::dev.off()

  # Axe 2 and 3
  grDevices::png(paste0(image_dir, parameters$analysis_name, "_PCA_axe2_3.png"), width=sizeImg*1.25, height=sizeImg*1.25)
  p=ggplot2::autoplot(stats::prcomp(lcpm2), data = lcpm3, size=15,colour="cond",show.legend = FALSE, x=2, y=3,shape=17)
  test=p + ggrepel::geom_label_repel(ggplot2::aes(label = rownames(lcpm3),fill = factor(lcpm3$cond)), color = 'white',size = 7) +
    ggplot2::labs(title="PCA Axes 2 and 3") + ggplot2::theme(
      plot.margin=ggplot2::margin(20,30,20,15),
      axis.text.y = ggplot2::element_text(face="bold",size=30),
      axis.text.x = ggplot2::element_text(face="bold",size=30),
      legend.text = ggplot2::element_text(size = 30),
      legend.title = ggplot2::element_text(size=30,face="bold"),
      strip.text.y = ggplot2::element_text(size=30, face="bold"),
      axis.title=ggplot2::element_text(size=30,face="bold"),
      plot.title = ggplot2::element_text(size=35, hjust=0.5)) +
    ggplot2::guides(colour=FALSE, fill = ggplot2::guide_legend(title = "Condition",override.aes = ggplot2::aes(label = ""), ncol=parameters$legendcol)) +
    ggplot2::theme(legend.position="bottom",legend.direction="vertical",legend.margin=ggplot2::margin(5,5,5,5),legend.box.margin=ggplot2::margin(10,10,10,10))
  print(test)
  grDevices::dev.off()

  # Axe 1 and 3
  grDevices::png(paste0(image_dir, parameters$analysis_name, "_PCA_axe1_3.png"), width=sizeImg*1.25, height=sizeImg*1.25)
  p=ggplot2::autoplot(stats::prcomp(lcpm2), data = lcpm3, size=15,colour="cond",show.legend = FALSE, x=1, y=3,shape=17)
  test=p + ggrepel::geom_label_repel(ggplot2::aes(label = rownames(lcpm3),fill = factor(lcpm3$cond)), color = 'white',size = 7) +
    ggplot2::labs(title="PCA Axes 1 and 3") + ggplot2::theme(
      plot.margin=ggplot2::margin(20,30,20,15),
      axis.text.y = ggplot2::element_text(face="bold",size=30),
      axis.text.x = ggplot2::element_text(face="bold",size=30),
      legend.text = ggplot2::element_text(size = 30),
      legend.title = ggplot2::element_text(size=30,face="bold"),
      strip.text.y = ggplot2::element_text(size=30, face="bold"),
      axis.title=ggplot2::element_text(size=30,face="bold"),
      plot.title = ggplot2::element_text(size=35, hjust=0.5)) +
    ggplot2::guides(colour=FALSE, fill = ggplot2::guide_legend(title = "Condition",override.aes = ggplot2::aes(label = ""), ncol=parameters$legendcol)) +
    ggplot2::theme(legend.position="bottom",legend.direction="vertical",legend.margin=ggplot2::margin(5,5,5,5),legend.box.margin=ggplot2::margin(10,10,10,10))
  print(test)
  grDevices::dev.off()

  # hierarchical clustering
  # mat.dist <- stats::dist(t(edgeR::cpm(asko_norm)), method = parameters$distcluts)
  mat.dist <- stats::dist(t(lcpm), method = parameters$distcluts)
  clustering <- stats::hclust(mat.dist, method=parameters$hclust)
  grDevices::png(paste0(image_dir, parameters$analysis_name, "_hclust.png"), width=sizeImg, height=sizeImg)
  graphics::par(oma=c(1,1,1,1))
  plot(clustering,
       main = 'Distances Correlation\nHierarchical clustering', sub="",
       xlab=paste0("Samples\nmethod hclust: ",parameters$hclust),
       hang = -1)
  j<-grDevices::dev.off()
}
