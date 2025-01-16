#' @title loadData
#'
#' @description
#' Function to load :
#' \itemize{
#'   \item Count data (one of this below):
#'   \itemize{
#'      \item count matrix : 1 file with all counts for each samples/conditions or multiple
#'      \item list of files : 1 file of count per conditions, files names contained in sample file
#'   }
#'   \item Metatdata :
#'   \itemize{
#'      \item sample file : file describing the samples and the experimental design
#'      \item contrast file : matrix which specifies which comparisons you would like to make between the samples
#'      \item (optional) annotation file : functional/genomic annotation for each genes
#'      \item (optional) GO terms annotations files : GO annotations for each genes
#'   }
#' }
#' Three output directory will be create :
#' \itemize{
#'   \item images : contains all the images created when running this pipeline with the exception of the upsetR and Venn graphs.
#'   \item vennDiagram : contain all venn diagrams created by VD function
#'   \item UpSetR_graphs : contain all upset graphs created by UpSetGraph function
#'   \item Askomics : files compatible with Askomics Software
#' }
#'
#' @param parameters list that contains all arguments charged in Asko_start
#' @return list contain all data and metadata (DGEList, samples descriptions, contrast, design and annotations)
#'
#' @examples
#' \dontrun{
#'    data<-loadData(parameters)
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
loadData <- function(parameters){

  # Folders for output files
  cat("\nCreated directories:\n")
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  if(dir.exists(study_dir)==FALSE){ dir.create(study_dir) }
  cat("\t",study_dir,"\n")

  explo_dir = paste0(study_dir,"DataExplore/")
  if(dir.exists(explo_dir)==FALSE){ dir.create(explo_dir) }
  cat("\t",explo_dir,"\n")

  norm_dir = paste0(study_dir,"NormCountsTables/")
  if(dir.exists(norm_dir)==FALSE){ dir.create(norm_dir) }
  cat("\t",norm_dir,"\n")

  de_dir = paste0(study_dir,"DEanalysis/")
  if(dir.exists(de_dir)==FALSE){ dir.create(de_dir) }
  cat("\t",de_dir,"\n")

  image_dir = paste0(de_dir, "DEimages/")
  if(dir.exists(image_dir)==FALSE){ dir.create(image_dir) }
  cat("\t",image_dir,"\n")

  table_dir = paste0(de_dir, "DEtables/")
  if(dir.exists(table_dir)==FALSE){ dir.create(table_dir) }
  cat("\t",table_dir,"\n")

  asko_dir = paste0(de_dir, "AskoTables/")
  if(dir.exists(asko_dir)==FALSE){ dir.create(asko_dir) }
  cat("\t",asko_dir,"\n")

  if(is.null(parameters$VD)==FALSE){
    venn_dir = paste0(study_dir, "VennDiagrams/")
    if(dir.exists(venn_dir)==FALSE){ dir.create(venn_dir) }
    cat("\t",venn_dir,"\n")
  }

  if((is.null(parameters$upset_basic)==FALSE) || (is.null(parameters$upset_list)==FALSE && is.null(parameters$upset_type)==FALSE)){
    upset_dir = paste0(study_dir, "UpsetGraphs/")
    if(dir.exists(upset_dir)==FALSE){ dir.create(upset_dir) }
    cat("\t",upset_dir,"\n")
  }

  # Management of input files
  input_path = paste0(parameters$dir_path, "/input/")

  # Sample file
  sample_path<-paste0(input_path, parameters$sample_file)
  samples<-utils::read.csv(sample_path, header=TRUE, sep="\t", row.names=1)

  # Selecting some sample (select_sample parameter)
  if(is.null(parameters$select_sample)==FALSE){
    if(parameters$regex==TRUE){
      selected<-c()
      for(sel in parameters$select_sample){
        select<-grep(sel, rownames(samples))
        if(is.null(selected)){selected=select}else{selected<-append(selected, select)}
      }
      samples<-samples[selected,]
    }
    else{ samples<-samples[parameters$select_sample,] }
  }

  # Deleting some samples (rm_sample parameter)
  if(is.null(parameters$rm_sample)==FALSE){
    if(parameters$regex==TRUE){
      for(rm in parameters$rm_sample){
        removed<-grep(rm, rownames(samples))
        if(length(removed!=0)){samples<-samples[-removed,]}
      }
    }
    else{
      for (rm in parameters$rm_sample) {
        rm2<-match(rm, rownames(samples))
        samples<-samples[-rm2,]
      }
    }
  }

  # Conditions and colors
  if(is.null(samples$color)==TRUE){
    condition<-unique(samples$condition)
    if(length(condition)<3){ color=c("#FF9999","#99CCFF") }
    else{ color<-grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(length(condition)) }
    samples$color<-NA
    j=0
    for(name in condition){
      j=j+1
      samples$color[samples$condition==name]<-color[j]
    }
  }
  else if(typeof(samples$colors)!="character"){
    color<-as.character(unlist(samples$color))
    samples$color<-color
  }

  # Count file(s).
  # Two possibilities:
  #     - 1 count file per condition
  #     - a matrix of count for all samples/conditions
  #
  # Multiple count files, 1 per conditions
  if(is.null(parameters$fileofcount)){
    cat("\nFiles of counts:\n")
    print(samples$file)
    cat("\nSamples:\n")
    print(rownames(samples))

    # creates a DGEList object from a table of counts
    dge<-edgeR::readDGE(paste0(input_path,samples$file), labels=rownames(samples), columns=c(parameters$col_genes,parameters$col_counts), header=TRUE, comment.char="#")
    dge<-edgeR::DGEList(counts=dge$counts, samples=samples)
  }
  # Matrix file with all counts for all conditions
  else{
    cat("\nSamples:\n")
    print(rownames(samples))
    count_path<-paste0(input_path, parameters$fileofcount)
    if(grepl(".csv", parameters$fileofcount)==TRUE){
      count<-utils::read.csv(count_path, header=TRUE, sep = "\t", row.names = parameters$col_genes, comment.char="#")
    }
    else{
      count<-utils::read.table(count_path, header=TRUE, sep = "\t", row.names = parameters$col_genes, comment.char="#")
    }

    # If you ask for some samples were removed for analysis
    select_counts<-row.names(samples)
    countT<-count[,select_counts]


    # Creates a DGEList object from a table of counts
    dge<-edgeR::DGEList(counts=countT, samples=samples)
  }

  # Experimental design
  Group<-factor(samples$condition)
  cat("\nConditions :\n")
  print(Group)
  designExp<-stats::model.matrix(~0+Group)
  rownames(designExp) <- row.names(samples)
  colnames(designExp) <- levels(Group)

  # Contrast for DE analysis
  contrast_path<-paste0(input_path, parameters$contrast_file)
  contrastab<-utils::read.table(contrast_path, sep="\t", header=TRUE, row.names = 1, comment.char="#", stringsAsFactors = FALSE)

  # Verify if some colunms will be not use for analysis
  rmcol<-list()
  for(condition_name in row.names(contrastab)){
    test<-match(condition_name, colnames(designExp),nomatch = 0)
    if(test==0){
      rm<-grep("0", contrastab[condition_name,], invert = TRUE)
      if(is.null(rmcol)){rmcol=rm}else{rmcol<-append(rmcol, rm)}
    }
  }
  # If it's the case then it delete them
  if (length(rmcol)> 0){
    rmcol<-unlist(rmcol)
    rmcol<-unique(rmcol)
    contrastab<-subset(contrastab, select=-rmcol)
  }
  contrastab2<-as.data.frame(contrastab[colnames(designExp),])
  colnames(contrastab2)<-colnames(contrastab)
  rownames(contrastab2)<-colnames(designExp)

  # Sort contrast table if more than one contrast in contrastab
  if(length(contrastab2)>1){
    ord<-match(colnames(designExp),row.names(contrastab2), nomatch = 0)
    contrast_table<-contrastab2[ord,]
  }
  else{
    contrast_table<-contrastab2
  }

  cat("\nContrasts:\n")
  print(contrast_table)

  # Format contrast : convert "+" to "1" and - to "-1"
  colnum<-c()
  for(contrast in colnames(contrast_table)){
    set_cond1<-row.names(contrast_table)[contrast_table[,contrast]=="+"]
    set_cond2<-row.names(contrast_table)[contrast_table[,contrast]=="-"]
    if(length(set_cond1)!=length(set_cond2)){
      contrast_table[,contrast][contrast_table[,contrast]=="+"]=signif(1/length(set_cond1),digits = 2)
      contrast_table[,contrast][contrast_table[,contrast]=="-"]=signif(-1/length(set_cond2),digits = 2)
    }
    else {
      contrast_table[,contrast][contrast_table[,contrast]=="+"]=1
      contrast_table[,contrast][contrast_table[,contrast]=="-"]=-1
    }
    contrast_table[,contrast]<-as.numeric(contrast_table[,contrast])
  }
  data<-list("dge"=dge, "samples"=samples, "contrast"=contrast_table, "design"=designExp)

  # Annnotation file
  if(is.null(parameters$annotation)==FALSE){
    annot<-utils::read.csv(paste0(input_path, parameters$annotation), header = TRUE, row.names = 1, sep = '\t', quote = "")
    data[["annot"]]=annot
  }
  return(data)
}
