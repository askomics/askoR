#' @title GOenrichment
#'
#' @description Perform GO enrichment analysis with topGO package.
#' This package provides tools for testing GO terms while accounting for
#' the topology of the GO graph. Different test statistics and different
#' methods for eliminating local similarities and dependencies between GO
#' terms can be implemented and applied.
#'
#' @param resDEG, data frame contains for each contrast the significance expression (1/0/-1) for all gene.
#' @param data_list, list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations).
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @return none.
#'
#' @import topGO
#' @import goSTAG
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#'    GOenrichment(resDEG, data_list, parameters)
#' }
#'
#' @export
GOenrichment<-function(resDEG, data_list, parameters){
  study_dir  = paste0(parameters$dir_path, "/", parameters$analysis_name, "/")
  input_path = paste0(parameters$dir_path, "/input/")
  img_go_dir = paste0(study_dir, "GO_images/")
  if(dir.exists(img_go_dir)==FALSE){
    dir.create(img_go_dir)
    cat("\n\nDirectory: ",img_go_dir," created\n")
  }

  if(is.null(parameters$GO)==TRUE){ return(NULL) }

  # Get GO annotations
  geneID2GO <- readMappings(file = paste0(input_path,parameters$geneID2GO_file))
  geneNames <- names(geneID2GO)

  for(contrast in colnames(data_list$contrast)){
    if(parameters$GO == "both"){
      print(contrast)
      geneSelected <- rownames(resDEG[apply(as.matrix(resDEG[,contrast]), 1, function(x) all(x!=0)),])
      titlename<-"all differentially expressed genes (up+down)"
    }else if(parameters$GO == "up"){
      geneSelected<-rownames(resDEG[apply(as.matrix(resDEG[,contrast]), 1, function(x) all(x==1)),])
      titlename<-"genes expressed UP"
    }else if(parameters$GO == "down"){
      geneSelected<-rownames(resDEG[apply(as.matrix(resDEG[,contrast]), 1, function(x) all(x==-1)),])
      titlename<-"genes expressed DOWN"
    }else{
      cat("\nBad value for GO parameters : autorized values are both, up, down or NULL.\n")
      return(NULL)
    }
    geneList <- factor(as.integer(geneNames %in% geneSelected))
    names(geneList) <- geneNames

    if(length(geneSelected)==0){
      cat("\nContrast:",contrast,"-> No DE genes found!\n")
      next
    }

    if(sum(levels(geneList)==1)==0){
      cat("\nContrast:",contrast,"-> No DE genes with GO annotation!\n")
      next
    }

    listOnto <- c("MF","BP","CC")
    for(ontology in listOnto){
      cat("\nContrast :",contrast," et ontologie :",ontology,"\n")
      GOdata <- methods::new("topGOdata",
                    nodeSize = parameters$GO_min_num_genes,
                    ontology = ontology,
                    allGenes = geneList,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)
      print(GOdata)

      resultTest <- runTest(GOdata, algorithm = parameters$GO_algo, statistic = parameters$GO_stats)
      print(resultTest)

      resGenTab <- GenTable(GOdata, numChar = 1000000, statisticTest = resultTest, orderBy = "statisticTest", topNodes=length(graph::nodes(graph(GOdata))) )
      resGenTab$Ratio = as.numeric(as.numeric(resGenTab$Significant)/as.numeric(resGenTab$Expected))
      resGenTab$GO_cat <- ontology

      if(ontology == "MF"){
        TabCompl<-resGenTab
        resGenTab[resGenTab=="< 1e-30"]<-"1.0e-30"

        if(nrow(resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold,])!=0){
          maxi<-parameters$GO_max_top_terms
          TabSigCompl<-resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold,]
          if(maxi > nrow(TabSigCompl)){ maxi<-nrow(TabSigCompl) }
          TabSigCompl<-TabSigCompl[1:maxi,]
        }else{
          cat("\n\n->",contrast," - ontology: ",ontology," - No enrichment can pe performed - there are no feasible GO terms!\n\n")
        }

      }else{
        TabCompl=rbind(TabCompl,resGenTab)
        resGenTab[resGenTab=="< 1e-30"]<-"1.0e-30"

        if(nrow(resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold,])!=0){
          maxi<-parameters$GO_max_top_terms
          tempSig<-resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold,]
          if(maxi > nrow(tempSig)){ maxi<-nrow(tempSig) }
          TabSigCompl=rbind(TabSigCompl,tempSig[1:maxi,])
        }else{
          cat("\n\n->",contrast," - ontology: ",ontology," - No enrichment can pe performed - there are no feasible GO terms!\n\n")
        }

      }
    }
    TabCompl<-TabCompl[TabCompl$Significant > 0,]
    utils::write.table(TabCompl, file=paste0(img_go_dir, parameters$analysis_name, "_", contrast, "_Complet_GOenrichment.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

    if(exists("TabSigCompl")==TRUE){
      if(nrow(TabSigCompl)>=1){
        if (parameters$GO_max_top_terms > 10) {
          TabSigCompl$Term = stringr::str_trunc(TabSigCompl$Term, 67)
        }else{
          TabSigCompl$Term = stringr::str_trunc(TabSigCompl$Term, 137)
        }
        # Graph for one contrast
        comp_names <- c( `MF` = "Molecular Function", `BP` = "Biological Process", `CC` = "Cellular Component")
        coul <- c(`MF` = "green4", `BP` = "red", `CC` = "blue")
        comp_names2 <- c(`MF` = "MF", `BP` = "BP", `CC` = "CC")

        TabSigCompl$Term = factor(TabSigCompl$Term, levels = unique(TabSigCompl$Term))
        minR=(min(TabSigCompl$Ratio)+max(TabSigCompl$Ratio))/4
        minP=(min(as.numeric(TabSigCompl$statisticTest))+max(as.numeric(TabSigCompl$statisticTest)))/4

        # Ratio Graph
        ggplot(TabSigCompl, aes(x=TabSigCompl$Ratio, y=TabSigCompl$Term, size=TabSigCompl$Significant, color=TabSigCompl$GO_cat)) +
          geom_point(alpha=1) + labs(title = paste0("GO Enrichment for contrast\n",contrast), x="Ratio Significant/Expected", y="GOterm") +
          scale_color_manual(values=coul,labels=comp_names,name="GO categories") +
          facet_grid(GO_cat~., scales="free", space = "free",labeller = as_labeller(comp_names2)) +
          scale_size_continuous(name="Number of genes") + scale_x_continuous(expand = expansion(add = minR)) +
          scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 70)) + theme_linedraw() +
          theme(
            panel.background = element_rect(fill = "grey90", colour = "grey90", size = 0.5, linetype = "solid"),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
            axis.text.y = element_text(face="bold", size=rel(0.75)),
            axis.text.x = element_text(face="bold", size=rel(0.75)),
            axis.title = element_text(face="bold", size=rel(0.75)),
            legend.title = element_text(size=rel(0.75), face="bold"),
            plot.title = element_text(face="bold", size=rel(1), hjust=1),
            legend.text = element_text(size=rel(0.5)))
        ggsave(filename=paste0(img_go_dir,contrast,"_Ratio_BUBBLESgraph.png"), width=7, height=7)

        # Pvalue Graph
        ggplot(TabSigCompl, aes(x=as.numeric(TabSigCompl$statisticTest), y=TabSigCompl$Term, size=TabSigCompl$Significant, color=TabSigCompl$GO_cat)) +
          geom_point(alpha=1) + labs(title = paste0("GO Enrichment for contrast\n",contrast),x="Pvalue",y="GOterm")+
          scale_color_manual(values=coul,labels=comp_names,name="GO categories")+
          facet_grid(TabSigCompl$GO_cat~., scales="free", space = "free",labeller = as_labeller(comp_names2))+
          scale_size_continuous(name="Number of genes") + scale_x_continuous(expand = expansion(add = minP)) +
          scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 70)) + theme_linedraw() +
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
        ggsave(filename=paste0(img_go_dir,contrast,"_Pvalue_BUBBLESgraph.png"), width=7, height=7)
      }
    }else{
      cat("\n\nToo few results to display the graph.\n\n")
    }
  }
}
