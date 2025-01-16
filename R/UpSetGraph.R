#' @title UpSetGraph
#'
#' @description Generate upsetR graphs.
#'
#' @param resDEG data frame contains for each contrast the significance expression (1/0/-1) for all gene.
#' @param data_list list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations).
#' @param parameters list that contains all arguments charged in Asko_start.
#' @return none
#'
#' @examples
#' \dontrun{
#'    UpSetGraph(sumDEG, data_list, parameters)
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
UpSetGraph <- function(resDEG, data_list, parameters){
  options(warn = -1)
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "UpsetGraphs/")
  if(dir.exists(image_dir)==FALSE){
    dir.create(image_dir)
    cat("\n\nDirectory: ",image_dir," created\n")
  }

  # Global UpsetR
  if(is.null(parameters$upset_basic)==FALSE){

    # created directory
    global_dir = paste0(image_dir, "Global_upset/")
    if(dir.exists(global_dir)==FALSE){
      dir.create(global_dir)
      cat("\n\nDirectory: ",global_dir," created\n\n")
    }

    if (parameters$upset_basic == "all"){
      cat("\nCreated global upset charts for all differentially expressed genes.")

      # verify empty groups
      if(sum(colSums(abs(resDEG)))==0){
        warning("Each group consists of none observation. Do you need to verify these empty groups?", immediate.=TRUE, call.=FALSE)
      }
      else{
        if(ncol(resDEG)<=6){tsc=2}else{tsc=1.45}
        # reoder columns by colsums value
        ordDEG<-abs(resDEG)
        ordDEG<-ordDEG[,order(colSums(-ordDEG, na.rm=TRUE))]
        sets<-colnames(resDEG)

        # all genes differentially expressed
        grDevices::png(paste0(global_dir, parameters$analysis_name,"_UpSetR_allDEG.png"), width=1280, height=1024)
        print(UpSetR::upset(data=ordDEG, sets=rev(sets), nsets=ncol(ordDEG), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA, text.scale = tsc))
        grid::grid.text("All differentially expressed genes (up+down)", x=0.65, y=0.95, gp=grid::gpar(fontsize=26))
        grDevices::dev.off()
      }
    }
    else if(parameters$upset_basic == "up"){
      cat("\nCreated global upset charts for genes expressed UP.")
      # table with Down Expressed Genes
      upDEG<-resDEG
      upDEG[upDEG==-1]<-0
      colnames(upDEG)<-gsub("vs"," > ",colnames(upDEG))

      # verify empty groups
      if(sum(colSums(abs(upDEG)))==0){
        warning("Each group consists of none observation. Do you need to verify these empty groups?", immediate.=TRUE, call.=FALSE)
      }
      else{
        # reoder columns by colsums value
        ordDEG<-upDEG[,order(colSums(-upDEG, na.rm=TRUE))]
        sets<-colnames(upDEG)

        # record upsetR graph for Down Expressed Genes
        if(ncol(ordDEG)<=6){tsc=2}else{tsc=1.45}
        grDevices::png(paste0(global_dir, parameters$analysis_name,"_UpSetR_upDEG.png"), width=1280, height=1024)
        print(UpSetR::upset(data=ordDEG, sets=rev(sets), nsets=ncol(ordDEG), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA, text.scale = tsc))
        grid::grid.text("Genes expressed \"UP\"", x=0.65, y=0.95, gp=grid::gpar(fontsize=26))
        grDevices::dev.off()
      }
    }
    else if(parameters$upset_basic == "down"){
      cat("\nCreated global upset charts for genes expressed DOWN.")
      # table with Up Expressed Genes
      downDEG<-resDEG
      downDEG[downDEG==1]<-0
      downDEG[downDEG==-1]<-1
      colnames(downDEG)<-gsub("vs"," < ",colnames(downDEG))

      # verify empty groups
      if(sum(colSums(abs(downDEG)))==0){
        warning("Each group consists of none observation. Do you need to verify these empty groups?", immediate.=TRUE, call.=FALSE)
      }
      else{
        # reoder columns by colsums value
        ordDEG<-downDEG[,order(colSums(-downDEG, na.rm=TRUE))]
        sets<-colnames(downDEG)

        # record upsetR graph for Up Expressed Genes
        if(ncol(ordDEG)<=6){tsc=2}else{tsc=1.45}
        grDevices::png(paste0(global_dir, parameters$analysis_name,"_UpSetR_downDEG.png"), width=1280, height=1024)
        print(UpSetR::upset(data=downDEG, sets=rev(colnames(downDEG)), nsets=ncol(downDEG), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA, text.scale = tsc))
        grid::grid.text("Genes expressed \"DOWN\"", x=0.65, y=0.95, gp=grid::gpar(fontsize=26))
        grDevices::dev.off()
      }
    }
    else if(parameters$upset_basic == "mixed"){
      cat("\nCreated global upset charts for genes expressed distinctly UP and DOWN.")
      # table with Up Expressed Genes
      upDEG<-resDEG
      upDEG[upDEG==-1]<-0
      colnames(upDEG)<-gsub("vs"," > ",colnames(upDEG))

      # table with Down Expressed Genes
      downDEG<-resDEG
      downDEG[downDEG==1]<-0
      downDEG[downDEG==-1]<-1
      colnames(downDEG)<-gsub("vs"," < ",colnames(downDEG))

      # table mixed up and down
      mixDEG<-cbind(upDEG,downDEG)
      sets<-as.vector(rbind(colnames(upDEG),colnames(downDEG)))
      metadata<-as.data.frame(cbind(c(colnames(upDEG),colnames(downDEG)),c(rep("UP",ncol(upDEG)),rep("DOWN",ncol(downDEG)))))
      names(metadata)<-c("sets", "SENS")

      # verify empty groups
      if(sum(colSums(abs(mixDEG)))==0){
        warning("Each group consists of none observation. Do you need to verify these empty groups?", immediate.=TRUE, call.=FALSE)
      }
      else{
        # reoder columns by colsums value
        ordDEG<-mixDEG[,order(colSums(-mixDEG, na.rm=TRUE))]

        # record upsetR graph for Up and Down Expressed Genes
        if(ncol(ordDEG)<=6){tsc=2}else{tsc=1.45}
        grDevices::png(paste0(global_dir, parameters$analysis_name,"_UpSetR_mixedDEG.png"), width=1280, height=1024)
        print(UpSetR::upset(data=ordDEG, sets=rev(sets), nsets=ncol(ordDEG), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA,
                            text.scale = tsc, set.metadata = list(data = metadata, plots = list(list(type = "matrix_rows",
                                                                                                     column = "SENS", colors = c(UP = "#FF9999", DOWN = "#99FF99"), alpha = 0.5)))))
        grid::grid.text("Genes expressed \"UP\" and \"DOWN\"", x=0.65, y=0.95, gp=grid::gpar(fontsize=26))
        grDevices::dev.off()
      }
    }
  }

  # Multiple graphs UpSetR
  if(is.null(parameters$upset_type)==TRUE && is.null(parameters$upset_list)==FALSE){
    warning("For Subset chart: upset_type must be not empty.\n", immediate.=TRUE, call.=FALSE)
  }
  else if(is.null(parameters$upset_type)==FALSE && is.null(parameters$upset_list)==TRUE){
    warning("For Subset chart: upset_list must be not empty.\n", immediate.=TRUE, call.=FALSE)
  }
  else if(is.null(parameters$upset_type)==FALSE && is.null(parameters$upset_list)==FALSE){
    # created directory
    subset_dir = paste0(image_dir, "Subset_upset/")
    if(dir.exists(subset_dir)==FALSE){
      dir.create(subset_dir)
      cat("\n\nDirectory: ",subset_dir," created\n")
    }

    cat("\nCreated upset charts for each element in \"upset_list\":",parameters$upset_type," expressed genes.\n")

    for(comparaison in parameters$upset_list){
      compa<-as.vector(limma::strsplit2(comparaison, "-"))
      cat("    -> Subset:",compa,"\n")

      if (parameters$upset_type == "all"){
        # verify empty groups
        if(sum(colSums(abs(resDEG[,compa])))==0){
          warning("Each group consists of none observation. Do you need to verify these empty groups?", immediate.=TRUE, call.=FALSE)
        }
        else{
          # reoder columns by colsums value
          ordDEG<-abs(resDEG[,compa])
          ordDEG<-ordDEG[,order(colSums(-ordDEG, na.rm=TRUE))]

          # record upsetR graph for all Differentially Expressed Genes
          if(length(compa)<=6){tsc=2}else{tsc=1.45}
          grDevices::png(paste0(subset_dir,comparaison,"_allDEG.png"), width=1280, height=1024)
          print(UpSetR::upset(data=ordDEG, sets=rev(compa), nsets=length(compa), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA, text.scale = tsc))
          grid::grid.text("All differentially expressed genes (up+down)", x=0.65, y=0.95, gp=grid::gpar(fontsize=26))
          grDevices::dev.off()
        }
      }
      else if(parameters$upset_type == "up"){
        # table with Down Expressed Genes
        upDEG<-resDEG
        upDEG[upDEG==-1]<-0
        colnames(upDEG)<-gsub("vs"," > ",colnames(upDEG))
        compa<-gsub("vs"," > ", compa)

        # verify empty groups
        if(sum(colSums(abs(upDEG[,compa])))==0){
          warning("Each group consists of none observation. Do you need to verify these empty groups?", immediate.=TRUE, call.=FALSE)
        }
        else{
          # reoder columns by colsums value
          ordDEG<-upDEG[,compa]
          ordDEG<-ordDEG[,order(colSums(-ordDEG, na.rm=TRUE))]

          # record upsetR graph for Down Expressed Genes
          if(length(compa)<=6){tsc=2}else{tsc=1.45}
          grDevices::png(paste0(subset_dir,comparaison,"_upDEG.png"), width=1280, height=1024)
          print(UpSetR::upset(data=ordDEG, sets=rev(compa), nsets=length(compa), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA, text.scale = tsc))
          grid::grid.text("Genes expressed \"UP\"", x=0.65, y=0.95, gp=grid::gpar(fontsize=26))
          grDevices::dev.off()
        }
      }
      else if(parameters$upset_type == "down"){
        # table with Up Expressed Genes
        downDEG<-resDEG
        downDEG[downDEG==1]<-0
        downDEG[downDEG==-1]<-1
        colnames(downDEG)<-gsub("vs"," < ",colnames(downDEG))
        newcompa<-gsub("vs"," < ",compa)

        # verify empty groups
        if(sum(colSums(abs(downDEG[,newcompa])))==0){
          warning("Each group consists of none observation. Do you need to verify these empty groups?", immediate.=TRUE, call.=FALSE)
        }
        else{
          # reoder columns by colsums value
          ordDEG<-downDEG[,newcompa]
          ordDEG<-ordDEG[,order(colSums(-ordDEG, na.rm=TRUE))]

          # record upsetR graph for Down Expressed Genes
          if(length(newcompa)<=6){tsc=2}else{tsc=1.45}
          grDevices::png(paste0(subset_dir,comparaison,"_downDEG.png"), width=1280, height=1024)
          print(UpSetR::upset(data=ordDEG, sets=rev(newcompa), nsets=length(newcompa), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA, text.scale = tsc))
          grid::grid.text("Genes expressed \"DOWN\"", x=0.65, y=0.95, gp=grid::gpar(fontsize=26))
          grDevices::dev.off()
        }
      }
      # mixed up and down expressed genes
      else if(parameters$upset_type == "mixed"){
        # table with Up Expressed Genes
        upDEG<-resDEG
        upDEG[upDEG==-1]<-0
        colnames(upDEG)<-gsub("vs"," > ",colnames(upDEG))
        compa1<-gsub("vs"," > ",compa)

        # table with Down Expressed Genes
        downDEG<-resDEG
        downDEG[downDEG==1]<-0
        downDEG[downDEG==-1]<-1
        colnames(downDEG)<-gsub("vs"," < ",colnames(downDEG))
        compa2<-gsub("vs"," < ",compa)

        # table mixed up and down
        mixDEG<-cbind(upDEG,downDEG)
        metadata<-as.data.frame(cbind(c(colnames(upDEG),colnames(downDEG)),c(rep("UP",ncol(upDEG)),rep("DOWN",ncol(downDEG)))))
        names(metadata)<-c("sets", "SENS")
        sets<-as.vector(rbind(colnames(upDEG[,compa1]),colnames(downDEG[,compa2])))

        # verify empty groups
        if(sum(colSums(abs(mixDEG[,sets])))==0){
          warning("Each group consists of none observation. Do you need to verify these empty groups?", immediate.=TRUE, call.=FALSE)
        }
        else{
          # reoder columns by colsums value
          ordDEG<-mixDEG[,sets]
          ordDEG<-ordDEG[,order(colSums(-ordDEG, na.rm=TRUE))]

          # record upsetR graph for Up and Down Expressed Genes
          if(length(sets)<=6){tsc=2}else{tsc=1.25}
          grDevices::png(paste0(subset_dir,comparaison,"_mixedDEG.png"), width=1280, height=1024)
          print(UpSetR::upset(data=ordDEG, sets=rev(sets), nsets=length(sets), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA,
                              text.scale = tsc, set.metadata = list(data = metadata,
                                                                    plots = list(list(type = "matrix_rows",column = "SENS", colors = c(UP = "#FF9999", DOWN = "#99FF99"), alpha = 0.5)))))
          grid::grid.text("Genes expressed \"UP\" and \"DOWN\"", x=0.65, y=0.95, gp=grid::gpar(fontsize=26))
          grDevices::dev.off()
        }
      }
    }
  }
}
