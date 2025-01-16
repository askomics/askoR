#' @title GOenrichment
#'
#' @description Perform GO enrichment analysis with topGO package.
#' This package provides tools for testing GO terms while accounting for
#' the topology of the GO graph. Different test statistics and different
#' methods for eliminating local similarities and dependencies between GO
#' terms can be implemented and applied.
#'
#' @param resDEG data frame contains for each contrast the significance expression (1/0/-1) for all gene.
#' @param data_list list contain all data and metadata (DGEList, samples descriptions, contrast, design and annotations).
#' @param parameters list that contains all arguments charged in Asko_start.
#' @param list gene list of interest if you want to apply GOenrichment function on a specific gene list.
#' @param title name of the gene list if you want to apply GOenrichment function on a specific gene list.
#' @return none.
#'
#' @import topGO
#' @import goSTAG
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#'    GOenrichment(resDEG, data_list, parameters)
#'    # OR
#'    GOenrichment(resDEG, data_list, parameters, list, title)
#'    # OR
#'    GOenrichment(resDEG, data_list, parameters, list=NULL, title=NULL)
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
GOenrichment<-function(resDEG, data_list, parameters, list=NULL, title=NULL){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  input_path = paste0(parameters$dir_path, "/input/")
  norm_dir = paste0(study_dir, "NormCountsTables/")
  GO_dir = paste0(study_dir, "/GOenrichment/")
  if(dir.exists(GO_dir)==FALSE){
    dir.create(GO_dir)
    cat("\n\nDirectory: ",GO_dir," created\n")
  }

  # Get GO annotations
  geneID2GO <- topGO::readMappings(file = paste0(input_path,parameters$geneID2GO_file))
  geneNames <- names(geneID2GO)

  if (is.null(list) == FALSE){
    img_go_dir = paste0(GO_dir,title,"/")
    if(dir.exists(img_go_dir)==FALSE){
      dir.create(img_go_dir)
      cat("Directory: ",img_go_dir," created\n")
    }
    GeneListName = title
    geneSelected = list
  }
  else if (parameters$GO == "up") {
    img_go_dir = paste0(GO_dir,"UP_DEgenes/")
    if(dir.exists(img_go_dir)==FALSE){
      dir.create(img_go_dir)
      cat("\n\nDirectory: ",img_go_dir," created\n")
    }
    GeneListName = colnames(data_list$contrast)
  }
  else if (parameters$GO == "down") {
    img_go_dir = paste0(GO_dir,"DOWN_DEgenes/")
    if(dir.exists(img_go_dir)==FALSE){
      dir.create(img_go_dir)
      cat("\n\nDirectory: ",img_go_dir," created\n")
    }
    GeneListName = colnames(data_list$contrast)
  }
  else {
    img_go_dir = paste0(GO_dir,"TOTAL_DEgenes/")
    if(dir.exists(img_go_dir)==FALSE){
      dir.create(img_go_dir)
      cat("\n\nDirectory: ",img_go_dir," created\n")
    }
    GeneListName = colnames(data_list$contrast)
  }

  for(contrast in GeneListName){

    # transfo resDEG vector to dataframe si on a un seul contraste
    if (ncol(resDEG)==1) {
      row = rownames(resDEG)
      resDEG = data.frame(x=row, contrast=resDEG[,1])
      colnames(resDEG) = c("NA",contrast)
      rownames(resDEG) = row
    }

    if (is.null(list) == TRUE){
      if(is.null(parameters$GO)==TRUE){ return(NULL) }


      if(parameters$GO == "both"){
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

    }

    geneList <- factor(as.integer(geneNames %in% geneSelected))
    names(geneList) <- geneNames

    img_GOtoGene_dir = paste0(img_go_dir, contrast,"_SignificantGO/")
    if(dir.exists(img_GOtoGene_dir)==FALSE){
      dir.create(img_GOtoGene_dir)
      cat("Directory: ",img_GOtoGene_dir," created\n")
    }

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

      resultTest <- topGO::runTest(GOdata, algorithm = parameters$GO_algo, statistic = parameters$GO_stats)

      resGenTab <- topGO::GenTable(GOdata, numChar = 1000000, statisticTest = resultTest, orderBy = "statisticTest", topNodes=length(graph::nodes(graph(GOdata))) )
      resGenTab$Ratio = as.numeric(as.numeric(resGenTab$Significant)/as.numeric(resGenTab$Expected))
      resGenTab$GO_cat <- ontology

      if (is.null(parameters$annotation)==FALSE){
        annot<-utils::read.csv(paste0(input_path, parameters$annotation), header = TRUE, row.names = 1, sep = '\t', quote = "")
      }

      # import normalized MEAN counts in CPM
      moys<-utils::read.csv(paste0(norm_dir, parameters$analysis_name,"_CPM_NormMeanCounts.txt"), header=TRUE, sep="\t", row.names=1)
      moys = as.matrix(moys)

      # Create files of genes for each enrichied GO
      myterms = as.character(resGenTab$GO.ID[as.numeric(resGenTab$statisticTest)<=parameters$GO_threshold])

      if (length(myterms) != "0"){
        cat("\nAskoR is saving one file per enriched GO-term (category ", ontology, ").\n")
        mygenes <- genesInTerm(GOdata, myterms)
        noms=names(mygenes)
        nomss=stringr::str_replace(noms,":","_")
        for (z in seq_len(length(mygenes))){
          listes=mygenes[[z]][mygenes[[z]] %in% geneSelected == TRUE]
          GOtab <- data.frame(Gene=listes)
          #GOtab$Gene_cluster = clustered
          rownames(GOtab)=GOtab$Gene
          if (is.null(parameters$annotation)==FALSE){
            GOtab = merge(GOtab, annot, by="row.names")
            GOtab = GOtab[,-1]
            GOtab = GOtab[,seq_len(2)]
            colnames(GOtab)[2] <- "Gene_description"
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
          utils::write.table(GOtab,paste0(img_GOtoGene_dir, ontology, "_", nomss[z],".txt"), sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)
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
          cat("\n\n->",contrast," - ontology: ",ontology," - No enrichment can pe performed - there are no feasible GO terms!\n\n")
          TabSigCompl<-as.data.frame(stats::setNames(replicate(8,numeric(0), simplify=FALSE),c("GO.ID","Term","Annotated","Significant","Expected","statisticTest","Ratio","GO_cat") ))
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
          cat("\n\n->",contrast," - ontology: ",ontology," - No enrichment can pe performed - there are no feasible GO terms!\n\n")
        }

      }

      ## Bargraph in each GO cat separately (ratio, pval, and number of genes)
      GoCoul="gray"

      if (is.null(list) == FALSE){
        GraphTitle0 = paste0("GO Enrichment (",ontology, " category)", "\n for list ", contrast, "\n (",length(which(geneList==1)), " annotated genes among ",length(geneSelected)," genes)")
      }
      else{
        GraphTitle0 = paste0("GO Enrichment (",ontology, " category)", "\n for contrast ", contrast, "\n (",length(which(geneList==1)), " annotated genes among ",length(geneSelected)," genes)")
      }

      if(exists("TabSigCompl")==TRUE){
        if(nrow(TabSigCompl[TabSigCompl$GO_cat==ontology,])>=1){
          TabOntology<-TabSigCompl[TabSigCompl$GO_cat==ontology,]
          ggplot2::ggplot(TabOntology, aes(x=stringr::str_wrap(TabOntology$Term, 55), y=TabOntology$Ratio,fill=-1*log10(as.numeric(TabOntology$statisticTest)))) +
            coord_flip()+
            geom_col()+
            theme_classic()+
            geom_text(aes(label=TabOntology$Significant), position=position_stack(0.5),color="white")+
            scale_fill_gradient(name="-log10pval",low=GoCoul,high=paste0(GoCoul,"4"))+
            scale_y_reverse()+
            labs(title = GraphTitle0, x="GOterm", y="Ratio Significant / Expected") +
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
          ggplot2::ggsave(filename=paste0(img_go_dir,contrast,"_",ontology,"_GOgraph.png"),width=10, height = 8)
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

        if (parameters$GO == "both"){
          GraphTitleList = paste0("GO Enrichment for list\n",contrast, "\n (Total DE)", "\n (",length(which(geneList==1)), " annotated genes among ",length(geneSelected)," genes)")
          GraphTitleContrast = paste0("GO Enrichment for contrast\n",contrast, "\n (Total DE)", "\n (",length(which(geneList==1)), " annotated genes among ",length(geneSelected)," genes)")
        }

        if (parameters$GO == "up"){
          GraphTitleList = paste0("GO Enrichment for list\n",contrast, "\n (UP DE)", "\n (",length(which(geneList==1)), " annotated genes among ",length(geneSelected)," genes)")
          GraphTitleContrast = paste0("GO Enrichment for contrast\n",contrast, "\n (UP DE)", "\n (",length(which(geneList==1)), " annotated genes among ",length(geneSelected)," genes)")
        }

        if (parameters$GO == "down"){
          GraphTitleList = paste0("GO Enrichment for list\n",contrast, "\n (DOWN DE)", "\n (",length(which(geneList==1)), " annotated genes among ",length(geneSelected)," genes)")
          GraphTitleContrast = paste0("GO Enrichment for contrast\n",contrast, "\n (DOWN DE)", "\n (",length(which(geneList==1)), " annotated genes among ",length(geneSelected)," genes)")
        }

        if (is.null(list) == FALSE){
          GraphTitle = GraphTitleList
        }
        else{
          GraphTitle = GraphTitleContrast
        }

        # Ratio Graph
        ggplot2::ggplot(TabSigCompl, aes(x=TabSigCompl$Ratio, y=TabSigCompl$Term, size=TabSigCompl$Significant, color=TabSigCompl$GO_cat)) +
          geom_point(alpha=1) +
          labs(title = GraphTitle, x="Ratio Significant/Expected", y="GOterm")+
          scale_color_manual(values=coul,labels=comp_names,name="GO categories") +
          facet_grid(TabSigCompl$GO_cat~., scales="free", space = "free",labeller = as_labeller(comp_names2)) +
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
        ggplot2::ggsave(filename=paste0(img_go_dir,contrast,"_Ratio_BUBBLESgraph.png"), width=7, height=7)

        # Pvalue Graph
        ggplot2::ggplot(TabSigCompl, aes(x=as.numeric(TabSigCompl$statisticTest), y=TabSigCompl$Term, size=TabSigCompl$Significant, color=TabSigCompl$GO_cat)) +
          geom_point(alpha=1) + labs(title = GraphTitle,x="Pvalue",y="GOterm")+
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
        ggplot2::ggsave(filename=paste0(img_go_dir,contrast,"_Pvalue_BUBBLESgraph.png"), width=7, height=7)
      }
    }else{
      cat("\n\nToo few results to display the graph.\n\n")
    }
  }
}
