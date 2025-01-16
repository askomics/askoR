#' @title GeneInfo_OnList
#'
#' @description
#' Produce a heatmap of the expression of the genes included in the list with information on DE status of each gene and resume all informations on the genes in a table
#' \itemize{
#'    \item Table with gene expression, DE status, clustering, and gene description
#'    \item Heatmap of gene expression and DE status
#' }
#'
#' @param list list contain all the genes you want to get information on
#' @param resDEG data frame contains for each contrast the significance expression (1/0/-1) for all genes coming from DEanalysis function.
#' @param data list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations)
#' @param parameters list that contains all arguments charged in Asko_start.
#' @param title name of the gene list
#' @param clustering data frame with clusters of each gene produced by ClustAndGO function
#' @param conditions list of the conditions you want to see in graph and table
#' @param contrasts list of the conditions you want to see in graph and table
#' @return none
#'
#' @import tidyverse
#' @import Rgraphviz
#' @import gghalves
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#'    GeneInfo_OnList(list, resDEG, data, parameters, title)
#'    # OR
#'    GeneInfo_OnList(list, resDEG, data, parameters, title, clustering)
#'    # OR
#'    GeneInfo_OnList(list, resDEG, data, parameters, title, clustering, conditions)
#'    # OR
#'    GeneInfo_OnList(list, resDEG, data, parameters, title, clustering, conditions, contrasts)
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
GeneInfo_OnList<-function(list, resDEG, data, parameters, title, clustering=NULL, conditions=NULL, contrasts=NULL){
  study_dir  = paste0(parameters$dir_path, "/", parameters$analysis_name, "/")
  input_path = paste0(parameters$dir_path, "/input/")
  norm_dir = paste0(study_dir, "NormCountsTables/")

  list_dir = paste0(study_dir, "GeneListExplore/")
  if(dir.exists(list_dir)==FALSE){
    dir.create(list_dir)
    cat("\n\nDirectory: ",list_dir," created\n")
  }

  img_InfosOnGenes_dir = paste0(list_dir, "/", title, "/")
  if(dir.exists(img_InfosOnGenes_dir)==FALSE){
    dir.create(img_InfosOnGenes_dir)
    cat("\n\nDirectory: ",img_InfosOnGenes_dir," created\n")
  }

  # import normalized MEAN counts in CPM
  moys<-utils::read.csv(paste0(norm_dir, parameters$analysis_name,"_CPM_NormMeanCounts.txt"), header=TRUE, sep="\t", row.names=1)
  if(is.null(conditions) == FALSE) {
    moys = moys[,colnames(moys) %in% conditions]
  }

  if(is.null(contrasts) == FALSE) {
    resDEG = resDEG[,colnames(resDEG) %in% contrasts]
  }

  moys = moys[rowSums(moys[])>0,]

  # for image size
  nsamples <- ncol(moys)
  sizeImg=15*nsamples
  if(sizeImg < 1024){ sizeImg=1024 }

  # Merge normalized MEAN counts in CPM and DE status from resDEG
  totalData = merge(moys,resDEG,by="row.names")
  rownames(totalData)=totalData[,1]
  totalData=totalData[,-1]

  if(is.null(clustering)==FALSE){
    names(clustering)="CoExpression_Cluster"
    `%notin%` <- Negate(`%in%`)
    clustNotDE=data.frame(Row.names=rownames(totalData[rownames(totalData) %notin% rownames(clustering),]))
    clustNotDE$CoExpression_Cluster="NOT DE"
    rownames(clustNotDE)=clustNotDE$'Row.names'
    clustNotDE2 = data.frame(clustNotDE[,-1])
    rownames(clustNotDE2)=rownames(clustNotDE)
    names(clustNotDE2)="CoExpression_Cluster"
    clust2=rbind(clustNotDE2,clustering)
    totalData = merge(totalData,clust2,by="row.names")
    rownames(totalData)=totalData[,1]
    totalData=totalData[,-1]
  }

  if(is.null(data$annot)==FALSE)
  {
    rnames<-rownames(totalData)                        # get Genes DE names
    annDE<-as.matrix(data$annot[rnames,])    # get annotations for each genes DE
    rownames(annDE)<-rnames
    colnames(annDE)<-colnames(data$annot)
    totalData<-cbind(totalData,annDE)                      # merge the two matrix
  }

  totalData=totalData[rownames(totalData) %in% list,]
  utils::write.table(totalData,paste0(img_InfosOnGenes_dir, title, "_SummaryGeneList.txt"), sep="\t", dec=".", row.names = TRUE, col.names = NA)

  if(is.null(conditions) == FALSE | is.null(contrasts) == FALSE) {
    suff = "_Subset"
  }
  else{
    suff = "_Complete"
  }

  n = ncol(moys)
  totalDataLOG = totalData[, seq_len(n)]
  mat = as.matrix(totalDataLOG)
  mat_scaled = t(apply(mat, 1, scale)) # same as : t(scale(t(mat)))
  colnames(mat_scaled)=colnames(mat)

  min = min(mat_scaled)
  max = max(mat_scaled)

  if (is.null(clustering)==FALSE){
    graphTitle = "Clusters"
  }
  else{
    graphTitle=""
  }

  col_fun = circlize::colorRamp2(c(min, 0, max), c("green", "white", "red"))

  hc = ComplexHeatmap::rowAnnotation("DE Status in contrasts" = as.matrix(totalData[,(n+1):(n+ncol(resDEG))]),simple_anno_size = unit(0.5, "cm"),gp=grid::gpar(pch=1,col="white",lwd = 4),col = list("DE Status in contrasts" = c("-1" = "green", "0" = "lightgrey", "1" = "red")),
                                     annotation_legend_param = list(
                                       at = c(-1, 0, 1),
                                       legend_height = unit(4, "cm"),
                                       title_position = "topleft",
                                       legend_side = "bottom", direction ="horizontal")
  )
  ComplexHeatmap::ht_opt("TITLE_PADDING" = unit(c(7, 7), "points"))
  if (is.null(clustering)==FALSE){
    ht_list = ComplexHeatmap::Heatmap((mat_scaled), name = "Expression \n(Scaled CPM)",
                                      heatmap_legend_param = list(
                                        #at = c(-2, 0, 2),
                                        legend_height = unit(4, "cm"),
                                        title_position = "topleft", direction = "horizontal"
                                      ),
                                      col = col_fun,
                                      row_split = totalData$CoExpression_Cluster,
                                      row_title_gp = grid::gpar(fill = grDevices::grey.colors(0.5), col="white", font = 2, fontsize=10),
                                      row_title_rot = 0,
                                      show_row_dend = FALSE,
                                      show_column_names = TRUE,
                                      show_column_dend = FALSE,
                                      row_names_side = "left",
                                      column_order = sort(colnames(mat)),
                                      row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
                                      right_annotation = hc,
                                      width=ncol(mat_scaled)*unit(20,"mm"),
                                      height=nrow(mat_scaled)*unit(5,"mm")
    )
  }
  else{
    ht_list = ComplexHeatmap::Heatmap((mat_scaled), name = "Expression \n(Scaled CPM)",
                                      heatmap_legend_param = list(
                                        #at = c(min, 0, max),
                                        legend_height = unit(4, "cm"),
                                        title_position = "topleft", direction = "horizontal"
                                      ),
                                      col = col_fun,
                                      show_row_dend = TRUE,
                                      show_column_names = TRUE,
                                      show_column_dend = FALSE,
                                      row_names_side = "left",
                                      column_order = sort(colnames(mat)),
                                      row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
                                      right_annotation = hc,
                                      width=ncol(mat_scaled)*unit(20,"mm"),
                                      height=nrow(mat_scaled)*unit(5,"mm")
    )
  }
  if (nrow(mat_scaled)<15){
    grDevices::png(paste0(img_InfosOnGenes_dir, title, suff, "_heatmap.png"), width=sizeImg*0.8, height=nrow(mat_scaled)*40)
  }
  else if (nrow(mat_scaled)>14 & nrow(mat_scaled)<30){
    grDevices::png(paste0(img_InfosOnGenes_dir, title, suff, "_heatmap.png"), width=sizeImg*0.8, height=nrow(mat_scaled)*25)
  }
  else if (nrow(mat_scaled)>29 & nrow(mat_scaled)<100){
    grDevices::png(paste0(img_InfosOnGenes_dir, title, suff, "_heatmap.png"), width=sizeImg*0.8, height=nrow(mat_scaled)*20)
  }
  else{
    grDevices::png(paste0(img_InfosOnGenes_dir, title, suff, "_heatmap.png"), width=sizeImg*0.8, height=nrow(mat_scaled)*17)
  }
  ComplexHeatmap::draw(ht_list, row_title = graphTitle, column_title_gp = grid::gpar(font=2, fontsize=15), heatmap_legend_side = "bottom",column_title = paste0("Expression and DE status \n on genes from list ", title, ""))
  grDevices::dev.off()

  if (length(list)!=nrow(totalData)){
    cat(paste0("WARNING: ", length(list)-nrow(totalData), " genes from the list (",length(list)," genes) were eliminated at the filtering step due to a very low expression level or because gene name is not valid. "))
  }

}
