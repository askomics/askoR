#' @title VD
#'
#' @description Plot Venn Diagram to compare different contrast
#'
#' @param resDEG data frame contains for each contrast the significance expression (1/0/-1) for all gene.
#' @param asko_list list of data.frame contain condition, contrast and context informations made by asko3c.
#' @param parameters list that contains all arguments charged in Asko_start.
#' @return none.
#'
#' @examples
#' \dontrun{
#'    VD(resDEG, parameters, asko_list)
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
VD <- function(resDEG, parameters, asko_list){
  options(warn = -1)
  # check parameters
  if(is.null(parameters$VD)==TRUE){ return(NULL) }
  if(is.null(parameters$compaVD)==TRUE ){ ## enleve car pb version R 4.4.1: || parameters$compaVD==""
    cat("compaVD parameter must not be empty!")
    return(NULL)
  }

  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  venn_dir = paste0(study_dir, "VennDiagrams/")
  if(dir.exists(venn_dir)==FALSE){
    dir.create(venn_dir)
    cat("\n\nDirectory: ",venn_dir," created\n")
  }

  # don't write log file
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

  cat("Created VennDiagrams ")
  if(parameters$VD == "all"){
    cat("for all differentially expressed genes.\n")
    for(comparaison in parameters$compaVD){
      compa<-limma::strsplit2(comparaison, "-")
      nbCompa<-length(compa)
      input<-list()

      if(nbCompa > 4) {
        cat("Comp : ",comparaison," - Accepts up to 4 comparisons!\n")
        next
      }

      for(n in compa){
        listDEG<-rownames(resDEG[apply(as.matrix(resDEG[,n]), 1, function(x) all(x!=0)),])
        nameDEG<-gsub("vs", "/", n)
        input[[nameDEG]]<-listDEG
      }

      color <- RColorBrewer::brewer.pal(nbCompa, parameters$palette)
      VennDiagram::venn.diagram(input,
                                main="All differentially expressed genes (up+down)",
                                filename=paste0(venn_dir, "/", comparaison, "_all.png"),
                                imagetype = "png",
                                main.cex = 1,
                                cat.cex = 0.8,
                                cex = 0.8,
                                fill = color,
                                category.names = labels(input),
                                col=0,euler.d = FALSE,scaled=FALSE)
    }
  }
  else if(parameters$VD == "both"){
    cat("for genes expressed UP and DOWN.\n")
    for(comparaison in parameters$compaVD){
      compa<-limma::strsplit2(comparaison, "-")
      nbCompa<-length(compa)
      input<-list()

      if(nbCompa != 2) {
        cat("Comp : ",comparaison," - Accepts only 2 comparisons!\n")
        next
      }

      for(n in compa){
        listUp<-rownames(resDEG[apply(as.matrix(resDEG[,n]), 1, function(x) all(x==1)),])
        listDown<-rownames(resDEG[apply(as.matrix(resDEG[,n]), 1, function(x) all(x==-1)),])
        nameUp<-gsub("vs", " > ", n)
        nameDown<-gsub("vs", " < ", n)
        input[[nameUp]]<-listUp
        input[[nameDown]]<-listDown
      }
      venn<-VennDiagram::venn.diagram(input, main="Genes expressed \"UP\" and \"DOWN\"",
                                      filename=paste0(venn_dir, "/", comparaison, "_mixed.png"),
                                      imagetype = "png",
                                      main.cex = 1,
                                      cat.cex = 0.8,
                                      cex=0.8,
                                      cat.dist = c(-0.4,-0.4,0.1,0.1),
                                      cat.col = c( "red1","royalblue1", "red3", "royalblue4"),
                                      category.names = labels(input),
                                      col=c( "red1","royalblue1", "red3", "royalblue4"),
                                      euler.d = FALSE,
                                      scaled=FALSE)
    }
  }
  else if(parameters$VD == "up"){
    cat("for genes expressed UP.\n")
    for(comparaison in parameters$compaVD){
      compa<-limma::strsplit2(comparaison, "-")
      nbCompa<-length(compa)
      input<-list()

      if(nbCompa > 4) {
        cat("Comp : ",comparaison," - Accepts up to 4 comparisons!\n")
        next
      }

      for(n in compa){
        listDEG<-rownames(resDEG[apply(as.matrix(resDEG[,n]), 1, function(x) all(x==1)),])
        nameDEG<-gsub("vs", " > ", n)
        input[[nameDEG]]<-listDEG
      }

      color <- RColorBrewer::brewer.pal(nbCompa, parameters$palette)
      VennDiagram::venn.diagram(input,
                                main="Genes expressed \"UP\"",
                                filename=paste0(venn_dir, "/", comparaison, "_up.png"),
                                imagetype = "png",
                                main.cex = 1,
                                cat.cex = 0.8,
                                cex = 0.8,
                                fill = color,
                                category.names = labels(input),
                                col=0,euler.d = FALSE,scaled=FALSE)
    }
  }
  else if(parameters$VD == "down"){
    cat("for genes expressed DOWN.\n")
    for(comparaison in parameters$compaVD){
      compa<-limma::strsplit2(comparaison, "-")
      nbCompa<-length(compa)
      input<-list()

      if(nbCompa > 4) {
        cat("Comp : ",comparaison," - Accepts up to 4 comparisons!\n")
        next
      }

      for(n in compa){
        listDEG<-rownames(resDEG[apply(as.matrix(resDEG[,n]), 1, function(x) all(x==-1)),])
        nameDEG<-gsub("vs", " < ", n)
        input[[nameDEG]]<-listDEG
      }

      color <- RColorBrewer::brewer.pal(nbCompa, parameters$palette)
      VennDiagram::venn.diagram(input,
                                main="Genes expressed \"DOWN\"",
                                filename=paste0(venn_dir, "/", comparaison, "_down.png"),
                                imagetype = "png",
                                main.cex = 1,
                                cat.cex = 0.8,
                                cex = 0.8,
                                fill = color,
                                category.names = labels(input),
                                col=0,euler.d = FALSE,scaled=FALSE)
    }
  }
}
