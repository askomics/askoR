#' @title ClustAndGO
#'
#' @description
#' Clusterize genes with same profile, proceed to GO-enrichment on clusters, search for contrasts enriched with genes of specific clusters, and identify intersections of DE list genes
#' \itemize{
#'    \item Graphs of clusters (heatmap and boxplot)
#'    \item Expression profiles in each cluster
#'    \item GO enrichments in each cluster
#'    \item Files with gene description of each significant enriched GO
#'    \item Over-representation of genes of each cluster in each contrast
#'    \item Intersections of DE list genes in each cluster
#' }
#'
#' @param asko_norm large DGEList with normalized counts by GEnorm function.
#' @param resDEG data frame contains for each contrast the significance expression (1/0/-1) for all genes coming from DEanalysis function.
#' @param parameters list that contains all arguments charged in Asko_start.
#' @param data list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations)
#' @param list gene list of interest if you want to apply ClustAndGO function on a specific gene list
#' @param title name of the gene list if you want to apply ClustAndGO function on a specific gene list
#' @return data frame with clusters of each gene
#'
#' @import topGO
#' @import goSTAG
#' @import tidyverse
#' @import Rgraphviz
#' @import gghalves
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#'    clust<-ClustAndGO(asko_norm, resDEG, parameters, data)
#'    # OR
#'    clust<-ClustAndGO(asko_norm, resDEG, parameters, data, list, title)
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
ClustAndGO <- function(asko_norm, resDEG, parameters, data, list=NULL, title=NULL){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  norm_dir = paste0(study_dir, "NormCountsTables/")

  CLUST_dir  = paste0(study_dir, "/Clustering/")
  if(dir.exists(CLUST_dir)==FALSE){
    dir.create(CLUST_dir)
    cat("\n\nDirectory: ",CLUST_dir," created\n")
  }

  input_path = paste0(parameters$dir_path, "/input/")

  if (is.null(list) == TRUE){
    img_Clustering_dir = paste0(CLUST_dir, "OnDEgenes/")
    if(dir.exists(img_Clustering_dir)==FALSE){
      dir.create(img_Clustering_dir)
      cat("Directory: ",img_Clustering_dir," created\n")
    }
  }
  else {
    img_Clustering_dir = paste0(CLUST_dir, title,"/")
    if(dir.exists(img_Clustering_dir)==FALSE){
      dir.create(img_Clustering_dir)
      cat("Directory: ",img_Clustering_dir," created\n")
    }
  }

  # for image size
  nsamples <- ncol(asko_norm$counts)
  sizeImg=15*nsamples
  if(sizeImg < 1024){ sizeImg=1024 }

  # import normalized MEAN counts in CPM
  moys<-utils::read.csv(paste0(norm_dir, parameters$analysis_name,"_CPM_NormMeanCounts.txt"), header=TRUE, sep="\t", row.names=1)
  moys = as.matrix(moys)

  # import normalized counts (all samples) in CPM
  object=utils::read.csv(paste0(norm_dir, parameters$analysis_name,"_CPM_NormCounts.txt"), header=TRUE, sep="\t", row.names=1)

  if (is.null(list) == TRUE){
    # keep only DE genes in at least "coseq_ContrastsThreshold" contrasts
    resDEG2=resDEG
    resDEG2[resDEG2== -1] <- 1
    object=object[which(rowSums(resDEG2)>=1),]
  }
  else {
    object=object[rownames(object) %in% list,]
    resDEG2 = resDEG[rownames(resDEG) %in% list,]
    resDEG2[resDEG2== -1] <- 1
  }

  if (parameters$coseq_data == 'LogScaledData'){
    object = log2(object+1)
    object = t(apply(object, 1, scale))
  }

  cat("Number of differentially expressed genes kept : ")
  print(nrow(object))

  conds=asko_norm$samples$condition

  ## run coseq ##
  if (length(parameters$coseq_ClustersNb)== 1){
    if (parameters$coseq_ClustersNb > 25){
      stop("TOO MANY CLUSTERS : Please set parameters$coseq_ClustersNb to default or under 26")
    }
  }
  else{
    if (length(parameters$coseq_ClustersNb) > 25) {
      stop("TOO MANY CLUSTERS : Please set parameters$coseq_ClustersNb to default or under a vector length of 26")
    }
  }

  if (parameters$coseq_data == 'LogScaledData' & parameters$coseq_transformation != 'none' ){
    stop("WRONG TRANSFORMATION CHOSEN  : You try to cluster data already transformed in log2+1. So you don't want to cluster transformed expression profiles (as recommended by coseq creators). Please set parameters$coseq_transformation to 'none' or parameters$coseq_data to 'ExpressionProfiles' ")
  }

  if (length(parameters$coseq_ClustersNb)== 1 & parameters$coseq_model=="kmeans"){
    coexpr=coseq::coseq(object, K=2:25, model = parameters$coseq_model, transformation = parameters$coseq_transformation,normFactors = "none", seed = 12345)
    clust=as.data.frame(coseq::clusters(coexpr, K=parameters$coseq_ClustersNb))
    names(clust)=c("clusters(coexpr)")
  }
  else{
    coexpr=coseq::coseq(object, K=parameters$coseq_ClustersNb, model = parameters$coseq_model, transformation = parameters$coseq_transformation,normFactors = "none", seed = 12345)
    clust=as.data.frame(coseq::clusters(coexpr))
    names(clust)=c("clusters(coexpr)")

    cat("\nSummary of CoSeq\n")
    print(summary(coexpr))
  }

  GeneToClusters<-merge(clust,moys,by="row.names")

  if (parameters$coseq_data == 'LogScaledData'){
    img_transfo_dir = paste0(img_Clustering_dir,parameters$coseq_model,"_OnLog2ScaledData_",length(unique(clust$`clusters(coexpr)`)),"clusters/")
    if(dir.exists(img_transfo_dir)==FALSE){
      dir.create(img_transfo_dir)
      cat("Directory: ",img_transfo_dir," created\n")
    }
  }
  else{
    img_transfo_dir = paste0(img_Clustering_dir,parameters$coseq_model,"_",parameters$coseq_transformation,"_",length(unique(clust$`clusters(coexpr)`)),"clusters/")
    if(dir.exists(img_transfo_dir)==FALSE){
      dir.create(img_transfo_dir)
      cat("Directory: ",img_transfo_dir," created\n")
    }
  }

  tempGeneToClusters = GeneToClusters
  rownames(tempGeneToClusters) = GeneToClusters$Row.names
  tempGeneToClusters = tempGeneToClusters[,-1]
  GeneToClustersSummary<-merge(tempGeneToClusters,resDEG,by="row.names")
  rownames(GeneToClustersSummary) = GeneToClustersSummary$Row.names
  GeneToClustersSummary = GeneToClustersSummary[,-1]

  if(is.null(data$annot)==FALSE)
  {
    rnames<-row.names(GeneToClustersSummary)                        # get Genes DE names
    annDE<-as.matrix(data$annot[rnames,])    # get annotations for each genes DE
    rownames(annDE)<-rnames
    colnames(annDE)<-colnames(data$annot)
    GeneToClustersSummary<-cbind(GeneToClustersSummary,annDE)                      # merge the two matrix

    utils::write.table(GeneToClustersSummary,paste0(img_transfo_dir, parameters$analysis_name,"_ClusteringSUMMARY_",parameters$coseq_model,"_",parameters$coseq_transformation,".txt"),sep="\t",dec=".",row.names = TRUE,col.names = NA)
  }
  else
  {
    utils::write.table(GeneToClustersSummary,paste0(img_transfo_dir, parameters$analysis_name,"_ClusteringSUMMARY_",parameters$coseq_model,"_",parameters$coseq_transformation,".txt"),sep="\t",dec=".",row.names = TRUE,col.names = NA)
  }

  ## Global graphs ##

  # Boxplots (scaled expression)
  GeneToClustersScaled=GeneToClusters
  GeneToClustersScaled=GeneToClustersScaled[,-2]
  rownames(GeneToClustersScaled)=GeneToClustersScaled$Row.names
  GeneToClustersScaled=GeneToClustersScaled[,-1]
  GeneToClustersScaled=t(apply(as.matrix(GeneToClustersScaled), 1, scale))
  colnames(GeneToClustersScaled)=colnames(GeneToClusters[,-c(seq_len(2))])

  final=data.frame()
  n=as.numeric(ncol(GeneToClustersScaled))
  for (i in seq_len(n)) {
    BDD <- data.frame(gene=rownames(GeneToClustersScaled))
    BDD$cluster=GeneToClusters$`clusters(coexpr)`
    BDD$expression=GeneToClustersScaled[,i]
    BDD$sample=colnames(GeneToClusters[i+2])
    final=rbind(final,BDD)
  }
  lab=c()
  for (x in unique(final$cluster)){
    lab=c(lab,paste0("Cluster ",x," (",nrow(GeneToClusters[GeneToClusters$`clusters(coexpr)`==x,])," genes)"))
  }
  names(lab)<-unique(final$cluster)
  ggplot2::ggplot(final, aes(x=final$sample, y=final$expression, fill=final$sample)) + geom_boxplot() +
    stat_summary(fun=mean, geom="line", aes(group=1), colour="red")+
    stat_summary(fun=mean, geom="point", colour="red")+
    facet_wrap(~final$cluster, labeller = as_labeller(lab))+
    theme_bw()+
    theme(strip.text.x = element_text(size=12),
          axis.text.x =element_blank(),
          axis.text.y=element_text(size=12),
          axis.ticks = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=15),
          legend.title = element_text(size=15,face="bold"),
          legend.text = element_text(size=12))+
    scale_y_continuous(name="Scaled expression")+
    scale_fill_discrete(name="Experimental \nconditions")
  if (length(unique(final$cluster)) > 3 & length(unique(final$cluster)) <= 6){
    ggplot2::ggsave(filename=paste0(img_transfo_dir, parameters$analysis_name, "_Boxplots_ScaledCPM_",parameters$coseq_model,"_",parameters$coseq_transformation,".png"),width=12,height=8)
  }
  else if (length(unique(final$cluster)) <= 3) {
    ggplot2::ggsave(filename=paste0(img_transfo_dir, parameters$analysis_name, "_Boxplots_ScaledCPM_",parameters$coseq_model,"_",parameters$coseq_transformation,".png"),width=12,height=4)
  }
  else {
    ggplot2::ggsave(filename=paste0(img_transfo_dir, parameters$analysis_name, "_Boxplots_ScaledCPM_",parameters$coseq_model,"_",parameters$coseq_transformation,".png"),width=12,height=12)
  }

  #Global expression profiles with probability (red genes are under proba 0.8)
  proba = 0.80
  for (thresh in proba){
    p <- coseq::plot(coexpr, graphs="profiles", K=length(unique(clust$`clusters(coexpr)`)), threshold=thresh)
    p2 <- p$profiles + ggtitle(paste0("Red genes are affiliated to the cluster with a probability lower than ",thresh))

    b <- coseq::plot(coexpr, graphs="probapost_barplots", K=length(unique(clust$`clusters(coexpr)`)), threshold=thresh)
    b2 <- b$probapost_barplots + ggtitle(paste0("Threshold probability = ",thresh)) +scale_fill_manual(values = c("black", "red"),name="Max\nconditional\nprobability")

    cowplot::plot_grid(p2, b2, ncol = 2, nrow = 1,rel_widths = c(1,1))
    ggplot2::ggsave(filename=paste0(img_transfo_dir, parameters$analysis_name, "_ClusterProbabilities_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Threshold",thresh,".png"),width=16,height=8)
  }


  # Heatmap on ScaledCPM
  m=n+2
  mat = as.matrix(GeneToClusters[, 3:m])
  mat_scaled = t(apply(mat, 1, scale))
  colnames(mat_scaled)=colnames(mat)
  rownames(mat_scaled)=GeneToClusters[,1]
  cluster=GeneToClusters[,2]

  if (parameters$coseq_HeatmapOrderSample==TRUE){
    ComplexHeatmap::ht_opt("TITLE_PADDING" = unit(c(8.5, 8.5), "points"))
    ht_list = ComplexHeatmap::Heatmap(t(mat_scaled), cluster_rows = FALSE, column_order=order(cluster), name = "Scaled CPM expression", column_split = cluster,
                                      heatmap_legend_param = list(title_position = "topcenter",legend_direction = "horizontal"),
                                      col = viridis::viridis(100),
                                      show_column_names = FALSE,
                                      column_title_gp = grid::gpar(fill = grDevices::grey.colors(0.5), col="white", font = 2, fontsize=15),
                                      row_gap = unit(2, "mm"), column_gap = unit(2, "mm")
    )
    grDevices::png(paste0(img_transfo_dir, parameters$analysis_name, "_Heatmap_ScaledCPM_",parameters$coseq_model,"_",parameters$coseq_transformation,"_MySampleOrder.png"), width=sizeImg*1.75, height=sizeImg/4*1.25)
    ComplexHeatmap::draw(ht_list,column_title_gp = grid::gpar(font=2, fontsize=20), heatmap_legend_side = "bottom",column_title = paste0("Heatmap on Clusters (parameters : ",parameters$coseq_model," and ",parameters$coseq_transformation," transformation)"))
    grDevices::dev.off()
  }
  else{
    ComplexHeatmap::ht_opt("TITLE_PADDING" = unit(c(8.5, 8.5), "points"))
    ht_list = ComplexHeatmap::Heatmap(t(mat_scaled),column_order=order(cluster), name = "Scaled CPM expression",column_split = cluster,
                                      heatmap_legend_param = list(title_position = "topcenter",legend_direction = "horizontal"),
                                      col = viridis::viridis(100),
                                      show_column_names = FALSE,
                                      column_title_gp = grid::gpar(fill = grDevices::grey.colors(0.5), col="white",font=2, fontsize=15),
                                      row_gap = unit(2, "mm"), column_gap = unit(2, "mm")
    )
    grDevices::png(paste0(img_transfo_dir, parameters$analysis_name, "_Heatmap_ScaledCPM_",parameters$coseq_model,"_",parameters$coseq_transformation,".png"), width=sizeImg*1.75, height=sizeImg/4*1.25)
    ComplexHeatmap::draw(ht_list,column_title_gp = grid::gpar(font=2, fontsize=20), heatmap_legend_side = "bottom",column_title = paste0("Heatmap on Clusters (parameters : ",parameters$coseq_model," and ",parameters$coseq_transformation," transformation)"))
    grDevices::dev.off()
  }

  ## Graphs for each cluster ##

  # import data and create vectors for color and cluster
  uniqClust=unique(GeneToClusters$`clusters(coexpr)`)

  GoCoul=c("palegreen", "skyblue", "lightsalmon", "thistle", "tan", "pink", "aquamarine", "violetred", "darkorange", "yellow", "mediumpurple", "wheat","palegreen", "skyblue", "lightsalmon", "thistle", "tan", "pink", "aquamarine", "violetred", "darkorange", "yellow", "mediumpurple", "wheat","palegreen")

  if (is.null(list) == TRUE){
    # create file with matrix of DE genes and cluster for each gene
    resDEG3=resDEG2[which(rowSums(resDEG2)>=1),]
    # delete contrasts with no DE genes
    if (ncol(resDEG2)==1){resDEG3 = as.matrix(resDEG3, ncol=2)}
    resDEG3=resDEG3[,(apply(resDEG3,2,sum)!=0)]
    if (ncol(resDEG2)==1){resDEG3 = as.matrix(resDEG3, ncol=2)}
  }
  else {
    if (ncol(resDEG2)==1){resDEG3 = as.matrix(resDEG3, ncol=2)}
    resDEG3=resDEG2[,(apply(resDEG2,2,sum)!=0)]
  }


  rownames(GeneToClusters)=GeneToClusters$`Row.names`
  ForContrast<-merge(resDEG3,GeneToClusters,by="row.names")

  if (is.null(list) == TRUE){
    FileForContrast=data.frame()
    for (z in 2:(ncol(resDEG3)+1)){
      tab <- data.frame(cluster=uniqClust)
      tab$contrast <- colnames(ForContrast[z])
      tab$TotalGenesInContrast <- sum(ForContrast[,z])
      tab$GenesOfContrastInCluster <- 0
      for (y in tab$cluster){
        ligne=which(tab$cluster==y)
        if (length(which((ForContrast$`clusters(coexpr)`==y) & (ForContrast[,z]=="1"))) >= 1) {
          tab$GenesOfContrastInCluster[ligne] <- length(which((ForContrast$`clusters(coexpr)`==y) & (ForContrast[,z]=="1")))
        }
      }
      tab$ObservedProportion <- paste0(round((tab$GenesOfContrastInCluster * 100 / tab$TotalGenesInContrast),1),"%")
      tab$ExpectedProportion <- 0
      tab$ChiTest <- ""
      FileForContrast=rbind(FileForContrast,tab)
    }
  }

  for (clustered in uniqClust){
    img_CLUST_dir = paste0(img_transfo_dir,"Cluster_",clustered,"/")
    if(dir.exists(img_CLUST_dir)==FALSE){
      dir.create(img_CLUST_dir)
      cat("Directory: ",img_CLUST_dir," created\n")
    }

    # Upset On each cluster
    cols=ncol(resDEG3)+1
    ForUpset = ForContrast[which(ForContrast$`clusters(coexpr)`==clustered),2:cols]
    if (ncol(resDEG2)==1){ForUpset = as.matrix(ForUpset, ncol=2)}
    ForUpset = ForUpset[,(apply(ForUpset,2,sum)!=0)]
    ForUpset = data.frame(ForUpset)

    if (ncol(ForUpset)>1) {
      grDevices::png(paste0(img_CLUST_dir,parameters$analysis_name,"_UpSet_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_",clustered,".png"), width=1600, height=1024, units = "px")
      print(UpSetR::upset(data=ForUpset, sets=rev(colnames(ForUpset)), nsets=ncol(ForUpset), keep.order=TRUE, att.color ="black" ,sets.bar.color=GoCoul[clustered],point.size = 5, line.size = 1.5, nintersects=NA, text.scale = 2))
      grid::grid.text(paste0("All differentially expressed genes (up+down) in cluster ",clustered), x=0.65, y=0.95, gp=grid::gpar(fontsize=20))
      grDevices::dev.off()
    }

    # Scaled expression of each condition in the cluster
    if (nrow(GeneToClusters[GeneToClusters$`clusters(coexpr)`==clustered,])<=750) {alph=1} else {alph=0.2}

    TabTempo<-final[which(final$cluster==clustered),]
    ggplot2::ggplot(TabTempo,aes(x=TabTempo$sample, y=TabTempo$expression))+
      geom_jitter(color = "gray50", alpha = alph, size = 1.5, show.legend = FALSE) +
      geom_violin(fill = GoCoul[clustered], alpha = 0.75, show.legend = FALSE) +
      stat_boxplot(geom = "errorbar", width = 0.15) +
      geom_boxplot(outlier.size = 0, width = 0.2, alpha = 0.75, show.legend = FALSE) +
      theme_bw() +
      labs(title = paste0("Scaled Expression of Cluster ",clustered, "\n (",nrow(GeneToClusters[GeneToClusters$`clusters(coexpr)`==clustered,])," genes)"), x="", y="Scaled Expression") +
      theme(legend.position = "none",
            axis.text.x =element_text(size=10,angle=90),
            axis.text.y=element_text(size=10),
            axis.ticks = element_blank(),
            plot.title = element_text(face="bold",size=15),
            axis.title.x=element_text(size=12),
            axis.title.y=element_text(size=12))
    ggplot2::ggsave(filename=paste0(img_CLUST_dir,parameters$analysis_name,"_ScaledExpression_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_",clustered,".png"),width=10, height=10)

    if (is.null(list) == TRUE){
      # Genes of each contrast in cluster and significance (Chi2) of enrichment of the cluster in each contrast
      FileForContrast$ExpectedProportion[FileForContrast$cluster==clustered] = length(which(GeneToClusters[,2]==clustered)) / nrow(resDEG3)
      proportion = length(which(GeneToClusters[,2]==clustered)) / nrow(resDEG3)
      pr=1-proportion

      for (a in which(FileForContrast$cluster==clustered)){
        b=FileForContrast$GenesOfContrastInCluster[a]
        d=FileForContrast$TotalGenesInContrast[a] - FileForContrast$GenesOfContrastInCluster[a]
        obs1=c(b,d)
        obs2=FileForContrast$GenesOfContrastInCluster[a]/FileForContrast$TotalGenesInContrast[a]
        proba=c(proportion,pr)
        if (stats::chisq.test(obs1,p=proba)$p.value<0.001 & obs2>=proportion) {
          FileForContrast$ChiTest[a]<-"***"
          FileForContrast$ObservedProportion[a]<-paste0(FileForContrast$ObservedProportion[a],FileForContrast$ChiTest[a])
        }
        else if (stats::chisq.test(obs1,p=proba)$p.value>=0.001 & stats::chisq.test(obs1,p=proba)$p.value<0.01  & obs2>=proportion){
          FileForContrast$ChiTest[a]<-"**"
          FileForContrast$ObservedProportion[a]<-paste0(FileForContrast$ObservedProportion[a],FileForContrast$ChiTest[a])
        }
        else if (stats::chisq.test(obs1,p=proba)$p.value>=0.01 & stats::chisq.test(obs1,p=proba)$p.value<0.05  & obs2>=proportion){
          FileForContrast$ChiTest[a]<-"*"
          FileForContrast$ObservedProportion[a]<-paste0(FileForContrast$ObservedProportion[a],FileForContrast$ChiTest[a])
        }
      }

      TabTempo<-FileForContrast[FileForContrast$cluster==clustered,]
      ggplot2::ggplot(TabTempo, aes(x=TabTempo$contrast, y=TabTempo$GenesOfContrastInCluster)) +
        coord_flip()+
        geom_col(fill=GoCoul[clustered])+
        theme_classic()+
        geom_text(aes(label=TabTempo$ObservedProportion), position=position_stack(0.5),color="black")+
        scale_y_reverse()+
        labs(title = paste0("DE Genes in contrasts for cluster ",clustered, "\n (",length(which(GeneToClusters[,2]==clustered))," genes in the cluster)"), x="Contrasts", y="Number of genes") +
        scale_x_discrete(position = "top")+
        theme(
          axis.text.y = element_text(face="bold",size=10),
          axis.text.x = element_text(face="bold",size=10),
          axis.title.x=element_text(face="bold",size=12),
          axis.title.y=element_blank(),
          legend.title = element_text(size=12,face="bold"),
          plot.title = element_text(face="bold",size=15),
          legend.text = element_text(size=12),
          panel.background = element_rect(colour = "black", size=0.5, fill=NA))
      ggplot2::ggsave(filename=paste0(img_CLUST_dir,parameters$analysis_name,"_GenesInContrasts_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_",clustered,".png"),width=10, height = 8)
    }


    # GO enrichment in the cluster for MF, CC, and BP category (if annotation file is provided)
    if(is.null(parameters$geneID2GO_file)==FALSE){
      img_GOtoGene_dir = paste0(img_CLUST_dir,"SignificantGO/")
      if(dir.exists(img_GOtoGene_dir)==FALSE){
        dir.create(img_GOtoGene_dir)
        cat("Directory: ",img_GOtoGene_dir," created\n")
      }

      geneID2GO <- readMappings(file = paste0(input_path,parameters$geneID2GO_file))
      geneNames <- names(geneID2GO)
      geneList <- factor(as.integer(geneNames %in% GeneToClusters[which(GeneToClusters$`clusters(coexpr)`==clustered),1]))
      names(geneList) <- geneNames

      if(nrow(GeneToClusters[which(GeneToClusters$`clusters(coexpr)`==clustered),])==0){
        cat("\n -> No DE genes found!\n")
        next
      }

      if(sum(levels(geneList)==1)==0){
        cat("\n -> No DE genes with GO annotation!\n")
        next
      }

      GO=NULL

      listOnto <- c("MF","BP","CC")
      for(ontology in listOnto){
        GOdata <- methods::new("topGOdata",
                               nodeSize = parameters$GO_min_num_genes,
                               ontology = ontology,
                               allGenes = geneList,
                               annot = annFUN.gene2GO,
                               gene2GO = geneID2GO)

        resultTest <- topGO::runTest(GOdata, algorithm = parameters$GO_algo, statistic = parameters$GO_stats)
        resGenTab <- topGO::GenTable(GOdata, numChar = 1000,statisticTest = resultTest, orderBy = "statisticTest", topNodes=length(graph::nodes(graph(GOdata))) )
        resGenTab$Ratio = as.numeric(as.numeric(resGenTab$Significant)/as.numeric(resGenTab$Expected))
        resGenTab$GO_cat <- ontology

        if (is.null(parameters$annotation)==FALSE){
          annot<-utils::read.csv(paste0(input_path, parameters$annotation), header = TRUE, row.names = 1, sep = '\t', quote = "")
        }

        myterms = as.character(resGenTab$GO.ID[as.numeric(resGenTab$statisticTest)<=parameters$GO_threshold])

        if (length(myterms) != "0"){
          cat("\nAskoR is saving one file per enriched GO-term in cluster ", clustered, " (category ", ontology, ").\n")
          mygenes <- genesInTerm(GOdata, myterms)
          noms=names(mygenes)
          nomss=stringr::str_replace(noms,":","_")
          for (z in seq_len(length(mygenes))){
            listes=mygenes[[z]][mygenes[[z]] %in% GeneToClusters[which(GeneToClusters$`clusters(coexpr)`==clustered),1] == TRUE]
            GOtab <- data.frame(Gene=listes)
            GOtab$Gene_cluster = clustered
            rownames(GOtab)=GOtab$Gene
            if (is.null(parameters$annotation)==FALSE){
              GOtab = merge(GOtab, annot, by="row.names")
              GOtab = GOtab[,-1]
              GOtab = GOtab[,seq_len(3)]
              colnames(GOtab)[3] <- "Gene_description"
              rownames(GOtab)=GOtab$Gene
            }
            else{
              GOtab$Gene_description = "No annotation file provided"
            }
            GOtab = merge(GOtab, resDEG, by="row.names")
            GOtab = GOtab[,-1]
            rownames(GOtab)=GOtab$Gene
            GOtab = merge(GOtab, moys, by="row.names")
            GOtab = GOtab[,-1]
            GOtab$GO_ID = noms[z]
            GOtab$GO_term = resGenTab[which(resGenTab$GO.ID==noms[z]),2]
            GOtab$GO_cat = resGenTab[which(resGenTab$GO.ID==noms[z]),8]
            utils::write.table(GOtab,paste0(img_GOtoGene_dir, ontology, "_", nomss[z],".txt"), sep="\t", dec=".", row.names = FALSE, col.names = TRUE, quote=FALSE)
          }
        }

        if(ontology == "MF"){
          TabCompl<-resGenTab
          resGenTab[resGenTab=="< 1e-30"]<-"1.0e-30"


          if(nrow(resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold & resGenTab$Significant >= parameters$GO_min_sig_genes,])!=0){
            maxi<-parameters$GO_max_top_terms
            TabSigCompl<-resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold & resGenTab$Significant >= parameters$GO_min_sig_genes,]
            if(maxi > nrow(TabSigCompl)){ maxi<-nrow(TabSigCompl) }
            TabSigCompl<-TabSigCompl[seq_len(maxi),]
          }else{
            cat("\n\n->Cluster ",clustered," - ontology: ",ontology," - No enrichment can pe performed - there are no feasible GO terms!\n\n")
            TabSigCompl<-as.data.frame(stats::setNames(replicate(8,numeric(0), simplify = FALSE),c("GO.ID","Term","Annotated","Significant","Expected","statisticTest","Ratio","GO_cat") ))
          }
        }else{
          TabCompl=rbind(TabCompl,resGenTab)
          resGenTab[resGenTab=="< 1e-30"]<-"1.0e-30"

          if(nrow(resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold & resGenTab$Significant >= parameters$GO_min_sig_genes,])!=0){
            maxi<-parameters$GO_max_top_terms
            tempSig<-resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold & resGenTab$Significant >= parameters$GO_min_sig_genes,]
            if(maxi > nrow(tempSig)){ maxi<-nrow(tempSig) }
            TabSigCompl=rbind(TabSigCompl,tempSig[seq_len(maxi),])
          }else{
            cat("\n\n->Cluster ",clustered," - ontology: ",ontology," - No enrichment can pe performed - there are no feasible GO terms!\n\n")
          }
        }

        if (ontology == "BP"){
          goCat= "Biological Process"
        }
        if (ontology == "CC"){
          goCat= "Cellular Component"
        }
        if (ontology == "MF"){
          goCat= "Molecular Function"
        }

        ## Bargraph in each GO cat separately (ratio, pval, and number of genes)
        if(exists("TabSigCompl")==TRUE){
          if(nrow(TabSigCompl[TabSigCompl$GO_cat==ontology,])>=1){
            TabTempo<-TabSigCompl[TabSigCompl$GO_cat==ontology,]
            ggplot2::ggplot(TabTempo, aes(x=stringr::str_wrap(TabTempo$Term, 55), y=TabTempo$Ratio,fill=-1*log10(as.numeric(TabTempo$statisticTest)))) +
              coord_flip()+
              geom_col()+
              theme_classic()+
              geom_text(aes(label=TabTempo$Significant), position=position_stack(0.5),color="white")+
              scale_fill_gradient(name="-log10pval",low=GoCoul[clustered],high=paste0(GoCoul[clustered],"4"))+
              scale_y_reverse()+
              labs(title = paste0("GO Enrichment in cluster ",clustered, " (", goCat, " category)", "\n (",length(which(geneList==1)), " annotated genes among the ",length(which(GeneToClusters[,2]==clustered))," in the cluster)"), x="GOterm", y="Ratio Significant / Expected") +
              scale_x_discrete(position = "top")+
              theme(
                axis.text.y = element_text(face="bold",size=10),
                axis.text.x = element_text(face="bold",size=10),
                axis.title.x=element_text(face="bold",size=12),
                axis.title.y=element_blank(),
                legend.title = element_text(size=12,face="bold"),
                plot.title = element_text(face="bold",size=15),
                legend.text = element_text(size=12),
                panel.background = element_rect(colour = "black", size=0.5, fill=NA))
            ggplot2::ggsave(filename=paste0(img_CLUST_dir,parameters$analysis_name,"_GOEnrichment_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_",clustered,"_", ontology,".png"),width=10, height = 8)
          }
        }
      }

      TabCompl<-TabCompl[TabCompl$Significant > 0,]
      utils::write.table(TabCompl, file=paste0(img_CLUST_dir,parameters$analysis_name,"_GOEnrichmentTable_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_",clustered, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

      ## Dotplot of all GO cat
      if(exists("TabSigCompl")==TRUE){
        if(nrow(TabSigCompl)>=1){
          if (parameters$GO_max_top_terms > 10) {
            TabSigCompl$Term = stringr::str_trunc(TabSigCompl$Term, 67)
          }else{
            TabSigCompl$Term = stringr::str_trunc(TabSigCompl$Term, 137)
          }
          comp_names <- c( `MF` = "Molecular Function", `BP` = "Biological Process", `CC` = "Cellular Component")
          coul <- c(`MF` = "green4", `BP` = "red", `CC` = "blue")
          comp_names2 <- c(`MF` = "MF", `BP` = "BP", `CC` = "CC")

          TabSigCompl$Term = factor(TabSigCompl$Term, levels = unique(TabSigCompl$Term))
          minR=(min(TabSigCompl$Ratio)+max(TabSigCompl$Ratio))/4
          minP=(min(as.numeric(TabSigCompl$statisticTest))+max(as.numeric(TabSigCompl$statisticTest)))/4

          # Ratio Graph
          ggplot2::ggplot(TabSigCompl, aes(x=TabSigCompl$Ratio, y=TabSigCompl$Term, size=TabSigCompl$Significant, color=TabSigCompl$GO_cat)) +
            geom_point(alpha=1) +
            labs(title = paste0("GO Enrichment for Cluster ",clustered, "\n (",length(which(geneList==1)), " annotated genes among the ",length(which(GeneToClusters[,2]==clustered))," in the cluster)"), x="Ratio Significant / Expected", y="GOterm") +
            scale_color_manual(values=coul,labels=comp_names,name="GO categories") +
            facet_grid(TabSigCompl$GO_cat~., scales="free", space = "free",labeller = as_labeller(comp_names2)) +
            scale_size_continuous(name="Number of genes") + scale_x_continuous(expand = expansion(add = minR)) +
            scale_y_discrete(labels = function(x) stringr::str_wrap(x, 70)) +
            theme_linedraw() + theme(
              panel.background = element_rect(fill = "grey93", colour = "grey93", size = 0.5, linetype = "solid"),
              panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
              panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
              axis.text.y = element_text(face="bold",size=8),
              axis.text.x = element_text(face="bold",size=10),
              legend.title = element_text(size=9,face="bold"),
              plot.title = element_text(face="bold",size=10),
              legend.text = element_text(size=9),
              strip.text.y = element_text(size=12, face="bold"))
          ggplot2::ggsave(filename=paste0(img_CLUST_dir,parameters$analysis_name,"_Ratio_BUBBLESgraph_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_",clustered, ".png"),width=10, height=10)

          # Pvalue Graph
          ggplot2::ggplot(TabSigCompl, aes(x=as.numeric(TabSigCompl$statisticTest), y=TabSigCompl$Term, size=TabSigCompl$Significant, color=TabSigCompl$GO_cat)) +
            geom_point(alpha=1) + labs(title = paste0("GO Enrichment for Cluster ",clustered, "\n (",length(which(geneList==1)), " annotated genes among the ",length(which(GeneToClusters[,2]==clustered))," in the cluster)"),x="Pvalue",y="GOterm")+
            scale_color_manual(values=coul,labels=comp_names,name="GO categories")+
            facet_grid(TabSigCompl$GO_cat~., scales="free", space = "free",labeller = as_labeller(comp_names2))+
            scale_size_continuous(name="Number of genes") + scale_x_continuous(expand = expansion(add = minP)) +
            scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 70)) + theme_linedraw() +
            scale_x_reverse()+
            theme(
              panel.background = element_rect(fill = "grey90", colour = "grey90", size = 0.5, linetype = "solid"),
              panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
              panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
              axis.text.y = element_text(face="bold", size=rel(0.75)),
              axis.text.x = element_text(face="bold", size=rel(0.75), angle=45, hjust=1),
              axis.title = element_text(face="bold", size=rel(0.75)),
              legend.title = element_text(size=rel(0.75), face="bold"),
              plot.title = element_text(face="bold", size=rel(1), hjust=1),
              legend.text = element_text(size=rel(0.5)))
          ggplot2::ggsave(filename=paste0(img_CLUST_dir,parameters$analysis_name,"_Pvalue_BUBBLESgraph_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_",clustered, ".png"),width=10, height=10)
        }
      }else{
        cat("\n\nToo few results to display the graph.\n\n")
      }
    }
  }
  return(clust)
}

