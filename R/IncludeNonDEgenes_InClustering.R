#' @title IncludeNonDEgenes_InClustering
#'
#' @description
#' Add a cluster with the genes that are not DE in the analysis to the ClustAndGO analysis
#' \itemize{
#'    \item Graphs of clusters (heatmap and boxplot) created through ClustAndGO function with an additionnal cluster représenting NON DE genes
#'    \item Expression profile of NON DE genes
#'    \item GO enrichments on NON DE genes
#'    \item Files with gene description of each significant enriched GO
#' }
#'
#' @param data list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations)
#' @param asko_norm large DGEList with normalized counts by GEnorm function.
#' @param resDEG data frame contains for each contrast the significance expression (1/0/-1) for all genes coming from DEanalysis function.
#' @param parameters list that contains all arguments charged in Asko_start.
#' @param clustering data frame with clusters of each gene produced by ClustAndGO function
#' @return none
#'
#' @import tidyverse
#' @import Rgraphviz
#' @import gghalves
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#'    IncludeNonDEgenes_InClustering(data, asko_norm, resDEG, parameters, clustering)
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
IncludeNonDEgenes_InClustering <- function(data, asko_norm, resDEG, parameters, clustering){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  input_path = paste0(parameters$dir_path, "/input/")
  norm_dir = paste0(study_dir, "NormCountsTables/")
  img_Clustering_dir = paste0(study_dir, "Clustering/OnDEgenes/")
  if(dir.exists(img_Clustering_dir)==FALSE){
    dir.create(img_Clustering_dir)
    cat("Directory: ",img_Clustering_dir," created\n")
  }


  if (parameters$coseq_data == 'LogScaledData'){
    img_transfo_dir = paste0(img_Clustering_dir,parameters$coseq_model,"_OnLog2ScaledData_",length(unique(clustering$`clusters(coexpr)`)),"clusters/")
    if(dir.exists(img_transfo_dir)==FALSE){
      dir.create(img_transfo_dir)
      cat("Directory: ",img_transfo_dir," created\n")
    }
  }
  else{
    img_transfo_dir = paste0(img_Clustering_dir,parameters$coseq_model,"_",parameters$coseq_transformation,"_",length(unique(clustering$`clusters(coexpr)`)),"clusters/")
    if(dir.exists(img_transfo_dir)==FALSE){
      dir.create(img_transfo_dir)
      cat("Directory: ",img_transfo_dir," created\n")
    }
  }

  img_CLUST_dir = paste0(img_transfo_dir,"NOT_DE/")
  if(dir.exists(img_CLUST_dir)==FALSE){
    dir.create(img_CLUST_dir)
    cat("Directory: ",img_CLUST_dir," created\n")
  }


  # for image size
  nsamples <- ncol(asko_norm$counts)
  sizeImg=15*nsamples
  if(sizeImg < 1024){ sizeImg=1024 }

  # import normalized MEAN counts in CPM
  moys<-utils::read.csv(paste0(norm_dir, parameters$analysis_name,"_CPM_NormMeanCounts.txt"), header=TRUE, sep="\t", row.names=1)

  # import GeneToClusters
  GeneToClusters<-utils::read.csv(paste0(img_transfo_dir, parameters$analysis_name,"_ClusteringSUMMARY_",parameters$coseq_model,"_",parameters$coseq_transformation,".txt"), header=TRUE, sep="\t", row.names=1)
  nbCond = length(unique(asko_norm$samples$condition))
  GeneToClusters = GeneToClusters[,seq_len(1+nbCond)]

  `%notin%` <- Negate(`%in%`)
  moysNotDE = moys[rownames(moys) %notin% rownames(GeneToClusters),]

  moys2=data.frame(Rnames=rownames(moysNotDE))
  moys2$clusters.coexpr.= 100
  rownames(moys2)=rownames(moysNotDE)
  moysNotDE=merge(moys2,moysNotDE,by="row.names")
  moysNotDE=moysNotDE[,-1]
  rownames(moysNotDE)=rownames(moys2)
  moysNotDE=moysNotDE[,-1]
  GeneToClusters=rbind(GeneToClusters,moysNotDE)

  GeneToClustersSummary = GeneToClusters
  GeneToClustersSummary[,1] = gsub(100, "NOT DE", GeneToClusters[,1])

  tempGeneToClusters = GeneToClustersSummary
  rownames(tempGeneToClusters) = rownames(GeneToClustersSummary)
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

    utils::write.table(GeneToClustersSummary, paste0(img_transfo_dir, parameters$analysis_name,"_ClusteringSUMMARY_WithNonDEgenes_",parameters$coseq_model,"_",parameters$coseq_transformation,".txt"),sep="\t",dec=".",row.names = TRUE,col.names = NA)
  }
  else
  {
    utils::write.table(GeneToClustersSummary, paste0(img_transfo_dir, parameters$analysis_name,"_ClusteringSUMMARY_WithNonDEgenes_",parameters$coseq_model,"_",parameters$coseq_transformation,".txt"),sep="\t",dec=".",row.names = TRUE,col.names = NA)
  }


  # Boxplots (scaled expression)
  GeneToClustersScaled=GeneToClusters
  GeneToClustersScaled=GeneToClustersScaled[,-1]
  GeneToClustersScaled=t(apply(as.matrix(GeneToClustersScaled), 1, scale))
  colnames(GeneToClustersScaled)=colnames(GeneToClusters[,-1])

  final=data.frame()
  n=as.numeric(ncol(GeneToClustersScaled))
  for (i in seq_len(n)) {
    BDD <- data.frame(gene=rownames(GeneToClustersScaled))
    BDD$cluster=GeneToClusters$clusters.coexpr.
    BDD$expression=GeneToClustersScaled[,i]
    BDD$sample=colnames(GeneToClusters[i+1])
    final=rbind(final,BDD)
  }
  lab=c()
  for (x in unique(final$cluster)){
    lab=c(lab,paste0("Cluster ",x," (",nrow(GeneToClusters[GeneToClusters$clusters.coexpr.==x,])," genes)"))
  }
  names(lab)<-unique(final$cluster)

  lab[[length(lab)]] = paste0("NOT DE (",nrow(GeneToClusters[GeneToClusters$clusters.coexpr.==100,])," genes)")

  ggplot(final,aes(x=sample, y=expression,fill=sample))+geom_boxplot()+
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
    ggsave(filename=paste0(img_transfo_dir, parameters$analysis_name, "_Boxplots_ScaledCPM_WithNonDEgenes_",parameters$coseq_model,"_",parameters$coseq_transformation,".png"),width=12,height=8)
  }
  else if (length(unique(final$cluster)) <= 3) {
    ggsave(filename=paste0(img_transfo_dir, parameters$analysis_name, "_Boxplots_ScaledCPM_WithNonDEgenes_",parameters$coseq_model,"_",parameters$coseq_transformation,".png"),width=12,height=4)
  }
  else {
    ggsave(filename=paste0(img_transfo_dir, parameters$analysis_name, "_Boxplots_ScaledCPM_WithNonDEgenes_",parameters$coseq_model,"_",parameters$coseq_transformation,".png"),width=12,height=12)
  }

  # Heatmap on ScaledCPM
  m=n+1
  mat = as.matrix(GeneToClusters[, 2:m])
  mat_scaled = t(apply(mat, 1, scale))
  colnames(mat_scaled)=colnames(mat)
  rownames(mat_scaled)=rownames(GeneToClusters)
  cluster=GeneToClusters[,1]

  if (parameters$coseq_HeatmapOrderSample==TRUE){
    ComplexHeatmap::ht_opt("TITLE_PADDING" = unit(c(8.5, 8.5), "points"))
    ht_list = ComplexHeatmap::Heatmap(t(mat_scaled),cluster_rows = FALSE,column_order=order(cluster), name = "Scaled CPM expression",column_split = cluster,
                                      heatmap_legend_param = list(title_position = "topcenter",legend_direction = "horizontal"),
                                      col = viridis::viridis(100),
                                      show_column_names = FALSE,
                                      column_title = c(c(seq_len(length(lab)-1)),"NOT DE"),
                                      column_title_gp = grid::gpar(fill = grDevices::grey.colors(0.5), col="white", font = 2, fontsize=15),
                                      row_gap = unit(2, "mm"), column_gap = unit(2, "mm")
    )
    grDevices::png(paste0(img_transfo_dir, parameters$analysis_name, "_Heatmap_ScaledCPM_WithNonDEgenes_",parameters$coseq_model,"_",parameters$coseq_transformation,"_MySampleOrder.png"), width=sizeImg*1.75, height=sizeImg/4*1.25)
    ComplexHeatmap::draw(ht_list,column_title_gp = grid::gpar(font=2, fontsize=20), heatmap_legend_side = "bottom",column_title = paste0("Heatmap on Clusters (parameters : ",parameters$coseq_model," and ",parameters$coseq_transformation," transformation)"))
    grDevices::dev.off()
  }
  else{
    ComplexHeatmap::ht_opt("TITLE_PADDING" = unit(c(8.5, 8.5), "points"))
    ht_list = ComplexHeatmap::Heatmap(t(mat_scaled),column_order=order(cluster), name = "Scaled CPM expression",column_split = cluster,
                                      heatmap_legend_param = list(title_position = "topcenter",legend_direction = "horizontal"),
                                      col = viridis::viridis(100),
                                      show_column_names = FALSE,
                                      column_title = c(c(seq_len(length(lab)-1)),"NOT DE"),
                                      column_title_gp = grid::gpar(fill = grDevices::grey.colors(0.5), col="white",font=2, fontsize=15),
                                      row_gap = unit(2, "mm"), column_gap = unit(2, "mm")
    )
    grDevices::png(paste0(img_transfo_dir, parameters$analysis_name, "_Heatmap_ScaledCPM_WithNonDEgenes_",parameters$coseq_model,"_",parameters$coseq_transformation,".png"), width=sizeImg*1.75, height=sizeImg/4*1.25)
    ComplexHeatmap::draw(ht_list, column_title_gp = grid::gpar(font=2, fontsize=20),heatmap_legend_side = "bottom",column_title = paste0("Heatmap on Clusters (parameters : ",parameters$coseq_model," and ",parameters$coseq_transformation," transformation)"))
    grDevices::dev.off()
  }

  # import data and create vectors for color and cluster
  clustered = 100

  GoCoul="gray"

  # Scaled expression of NON DE genes
  if (nrow(GeneToClusters[GeneToClusters$clusters.coexpr.==clustered,])<=750) {alph=1} else {alph=0.2}

  ggplot(final[which(final$cluster==clustered),],aes(x=sample, y=expression))+
    geom_jitter(color = "gray50", alpha = alph, size = 1.5, show.legend = FALSE) +
    geom_violin(fill = GoCoul, alpha = 0.75, show.legend = FALSE) +
    stat_boxplot(geom = "errorbar", width = 0.15) +
    geom_boxplot(outlier.size = 0, width = 0.2, alpha = 0.75, show.legend = FALSE) +
    theme_bw() +
    labs(title = paste0("Scaled Expression of Cluster NOT DE \n(",nrow(GeneToClusters[GeneToClusters$clusters.coexpr.==clustered,])," genes)"), x="", y="Scaled Expression") +
    theme(legend.position = "none",
          axis.text.x =element_text(size=10,angle=90),
          axis.text.y=element_text(size=10),
          axis.ticks = element_blank(),
          plot.title = element_text(face="bold",size=15),
          axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12))
  ggsave(filename=paste0(img_CLUST_dir,parameters$analysis_name,"_ScaledExpression_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_NOT_DE.png"),width=10, height=10)



  # GO enrichment in the cluster for MF, CC, and BP category
  if(is.null(parameters$geneID2GO_file)==FALSE){

    img_GOtoGene_dir = paste0(img_CLUST_dir,"SignificantGO/")
    if(dir.exists(img_GOtoGene_dir)==FALSE){
      dir.create(img_GOtoGene_dir)
      cat("Directory: ",img_GOtoGene_dir," created\n")
    }

    geneID2GO <- readMappings(file = paste0(input_path,parameters$geneID2GO_file))
    geneNames <- names(geneID2GO)
    geneList <- factor(as.integer(geneNames %in% rownames(GeneToClusters[which(GeneToClusters$clusters.coexpr.==clustered),])))

    names(geneList) <- geneNames

    if(nrow(GeneToClusters[which(GeneToClusters$clusters.coexpr.==clustered),])==0){
      cat("\n -> No DE genes found!\n")
    }

    if(sum(levels(geneList)==1)==0){
      cat("\n -> No DE genes with GO annotation!\n")
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

      resultTest <- runTest(GOdata, algorithm = parameters$GO_algo, statistic = parameters$GO_stats)
      resGenTab <- GenTable(GOdata, numChar = 1000,statisticTest = resultTest, orderBy = "statisticTest", topNodes=length(graph::nodes(graph(GOdata))) )
      resGenTab$Ratio = as.numeric(as.numeric(resGenTab$Significant)/as.numeric(resGenTab$Expected))
      resGenTab$GO_cat <- ontology

      if (is.null(parameters$annotation)==FALSE){
        annot<-utils::read.csv(paste0(input_path, parameters$annotation), header = TRUE, row.names = 1, sep = '\t', quote = "")
      }

      myterms = as.character(resGenTab$GO.ID[as.numeric(resGenTab$statisticTest)<=parameters$GO_threshold])


      if (length(myterms) != "0"){
        cat("\nAskoR is saving one file per enriched GO-term in cluster NOT DE (category ", ontology, ").\n")
        mygenes <- genesInTerm(GOdata, myterms)
        noms=names(mygenes)
        nomss=stringr::str_replace(noms,":","_")
        for (z in seq_len(length(mygenes))){
          listes=mygenes[[z]][mygenes[[z]] %in% rownames(GeneToClusters[which(GeneToClusters$clusters.coexpr.==clustered),]) == TRUE]
          GOtab <- data.frame(Gene=listes)
          GOtab$Gene_cluster = "NOT DE"
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
          cat("\n\n->Cluster NOT DE - ontology: ",ontology," - No enrichment can pe performed - there are no feasible GO terms!\n\n")
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
          cat("\n\n->Cluster NOT DE - ontology: ",ontology," - No enrichment can pe performed - there are no feasible GO terms!\n\n")
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
          ggplot(TabTempo, aes(x=stringr::str_wrap(TabTempo$Term, 55), y=TabTempo$Ratio,fill=-1*log10(as.numeric(TabTempo$statisticTest)))) +
            coord_flip()+
            geom_col()+
            theme_classic()+
            geom_text(aes(label=TabTempo$Significant), position=position_stack(0.5),color="white")+
            scale_fill_gradient(name="-log10pval",low=GoCoul,high=paste0(GoCoul,"4"))+
            scale_y_reverse()+
            labs(title = paste0("GO Enrichment in cluster NOT DE (", goCat, " category)", "\n (",length(which(geneList==1)), " annotated genes among the ",length(which(GeneToClusters[,1]==clustered))," in the cluster)"), x="GOterm", y="Ratio Significant / Expected") +
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
          ggsave(filename=paste0(img_CLUST_dir,parameters$analysis_name,"_GOEnrichment_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_NOT_DE_", ontology,".png"),width=10, height = 8)
        }
      }
    }

    TabCompl<-TabCompl[TabCompl$Significant > 0,]
    utils::write.table(TabCompl, file=paste0(img_CLUST_dir,parameters$analysis_name,"_GOEnrichmentTable_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_NOT_DE.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

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
        ggplot(TabSigCompl, aes(x=TabSigCompl$Ratio, y=TabSigCompl$Term, size=TabSigCompl$Significant, color=TabSigCompl$GO_cat)) +
          geom_point(alpha=1) +
          labs(title = paste0("GO Enrichment for Cluster NOT DE \n(",length(which(geneList==1)), " annotated genes among the ",length(which(GeneToClusters[,1]==clustered))," in the cluster)"), x="Ratio Significant / Expected", y="GOterm") +
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
        ggsave(filename=paste0(img_CLUST_dir,parameters$analysis_name,"_Ratio_BUBBLESgraph_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_NOT_DE.png"),width=10, height=10)

        ggplot2::ggplot(TabSigCompl, aes(x=as.numeric(TabSigCompl$statisticTest), y=TabSigCompl$Term, size=TabSigCompl$Significant, color=TabSigCompl$GO_cat)) +
          geom_point(alpha=1) + labs(title = paste0("GO Enrichment for Cluster NOT DE \n(",length(which(geneList==1)), " annotated genes among the ",length(which(GeneToClusters[,1]==clustered))," in the cluster)"), x="Pvalue", y="GOterm") +
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
