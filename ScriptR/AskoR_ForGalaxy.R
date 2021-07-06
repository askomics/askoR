#' @title Asko_start
#'
#' @description Initialize and Scans parameters from command line in a python-like style:
#' \itemize{
#'    \item declare options, their flags, types, default values and help messages,
#'    \item read the arguments passed to the R script and parse them according to what has been declared in step 1.
#' }
#' Parameters can be called by their names as declared in opt object.
#'
#' @return List of parameters that contains all arguments.
#'
#' @examples
#'    parameters <- Asko_start()
#'    parameters$threshold_cpm <- 1  # Set parameters threshold cpm to new value
#'
#' @note All parameters were describe in README documentation
#'
#' @export
Asko_start <- function(){
  # Loading libraries in silent mode (only error messages will be displayed)
  pkgs<-c("limma","statmod","edgeR","VennDiagram","RColorBrewer","UpSetR","grid","topGO","ggfortify","gghalves","tidyverse",
          "Rgraphviz","ggplot2","ggrepel","gplots","stringr","optparse","goSTAG","Glimma","ComplexHeatmap","cowplot","circlize","corrplot")
  for(p in pkgs) suppressWarnings(suppressMessages(library(p, quietly=TRUE, character.only=TRUE, warn.conflicts=FALSE)))

  # Specify desired options in a list
  option_list = list(
    optparse::make_option(c("-o", "--out"), type="character", default="AskoRanalysis",dest="analysis_name",
                          help="output directory name [default= %default]", metavar="character"),
    optparse::make_option(c("-d", "--dir"), type="character", default=".",dest="dir_path",
                          help="data directory path [default= %default]", metavar="character"),
    optparse::make_option(c("-O", "--org"), type="character", default="Asko", dest="organism",
                          help="Organism name [default= %default]", metavar="character"),
    optparse::make_option(c("-f", "--fileofcount"), type="character", default=NULL, dest="fileofcount",
                          help="file of counts [default= %default]", metavar="character"),
    optparse::make_option(c("-G", "--col_genes"), type="integer", default=1, dest="col_genes",
                          help="col of genes ids in count files [default= %default]", metavar="integer"),
    optparse::make_option(c("-C", "--col_counts"), type="integer", default=7,dest="col_counts",
                          help="col of counts in count files [default= %default (featureCounts) ]", metavar="integer"),
    optparse::make_option(c("-t", "--sep"), type="character", default="\t", dest="sep",
                          help="field separator for count files or count matrix [default= %default]", metavar="character"),
    optparse::make_option(c("-a", "--annotation"), type="character", default=NULL, dest="annotation",
                          help="annotation file [default= %default]", metavar="character"),
    optparse::make_option(c("-s", "--sample"), type="character", default="Samples.txt", dest="sample_file",
                          help="Samples file [default= %default]", metavar="character"),
    optparse::make_option(c("-c", "--contrasts"), type="character", default="Contrasts.txt",dest="contrast_file",
                          help="Contrasts file [default= %default]", metavar="character"),
    optparse::make_option(c("-k", "--mk_context"), type="logical", default=FALSE,dest="mk_context",
                          help="generate automatically the context names [default= %default]", metavar="logical"),
    optparse::make_option(c("--palette"), type="character", default="Set2", dest="palette",
                          help="Color palette (ggplot)[default= %default]", metavar="character"),
    optparse::make_option(c("-R", "--regex"), type="logical", default=FALSE, dest="regex",
                          help="use regex when selecting/removing samples [default= %default]", metavar="logical"),
    optparse::make_option(c("-S", "--select"), type="character", default=NULL, dest="select_sample",
                          help="selected samples [default= %default]", metavar="character"),
    optparse::make_option(c("-r", "--remove"), type="character", default=NULL, dest="rm_sample",
                          help="removed samples [default= %default]", metavar="character"),
    optparse::make_option(c("--th_cpm"), type="double", default=0.5, dest="threshold_cpm",
                          help="CPM's threshold [default= %default]", metavar="double"),
    optparse::make_option(c("--rep"), type="integer", default=3, dest="replicate_cpm",
                          help="Minimum number of replicates [default= %default]", metavar="integer"),
    optparse::make_option(c("--th_FDR"), type="double", default=0.05, dest="threshold_FDR",
                          help="FDR threshold [default= %default]", metavar="double"),
    optparse::make_option(c("-n", "--normalization"), type="character", default="TMM", dest="normal_method",
                          help="normalization method (TMM/RLE/upperquartile/none) [default= %default]", metavar="character"),
    optparse::make_option(c("--adj"), type="character", default="fdr", dest="p_adj_method",
                          help="p-value adjust method (holm/hochberg/hommel/bonferroni/BH/BY/fdr/none) [default= %default]", metavar="character"),
    optparse::make_option("--glm", type="character", default="qlf", dest="glm",
                          help="GLM method (lrt/qlf) [default= %default]", metavar="character"),
    optparse::make_option("--glmDisp", type="logical", default=FALSE, dest="glm_disp",
                          help="Estimate Common, Trended and Tagwise Negative Binomial dispersions GLMs (TRUE/FALSE) [default= %default]", metavar="logical"),
    optparse::make_option(c("--lfc"), type="logical", default=TRUE, dest="logFC",
                          help="logFC in the summary table [default= %default]", metavar="logical"),
    optparse::make_option(c("--th_lfc"), type="double", default=0, dest="threshold_logFC",
                          help="logFC threshold [default= %default]", metavar="double"),
    optparse::make_option("--CompleteHm", type="logical", default=FALSE, dest="CompleteHeatmap",
                          help="generation of the normalized expression heatmap on ALL genes (TRUE/FALSE) [default= %default]", metavar="logical"),
    optparse::make_option("--fc", type="logical", default=TRUE, dest="FC",
                          help="FC in the summary table [default= %default]", metavar="logical"),
    optparse::make_option(c("--lcpm"), type="logical", default=FALSE, dest="logCPM",
                          help="logCPm in the summary table [default= %default]", metavar="logical"),
    optparse::make_option("--fdr", type="logical", default=TRUE, dest="FDR",
                          help="FDR in the summary table [default= %default]", metavar="logical"),
    optparse::make_option("--lr", type="logical", default=FALSE, dest="LR",
                          help="LR in the summary table [default= %default]", metavar="logical"),
    optparse::make_option(c("--sign"), type="logical", default=TRUE, dest="Sign",
                          help="Significance (1/0/-1) in the summary table [default= %default]", metavar="logical"),
    optparse::make_option(c("--expr"), type="logical", default=TRUE, dest="Expression",
                          help="Significance expression in the summary table [default= %default]", metavar="logical"),
    optparse::make_option(c("--mc"), type="logical", default=TRUE, dest="mean_counts",
                          help="Mean counts in the summary table [default= %default]", metavar="logical"),
    optparse::make_option(c("--dclust"), type="character", default="euclidean", dest="distcluts",
                          help="The distance measure to be used : euclidean, maximum, manhattan, canberra, binary or minkowski [default= %default]", metavar="character"),
    optparse::make_option(c("--hclust"), type="character", default="complete", dest="hclust",
                          help="The agglomeration method to be used : ward.D, ward.D2, single, complete, average, mcquitty, median or centroid [default= %default]", metavar="character"),
    optparse::make_option(c("--hm"), type="logical", default=TRUE, dest="heatmap",
                          help="generation of the expression heatmap [default= %default]", metavar="logical"),
    optparse::make_option(c("--nh"), type="integer", default="50", dest="numhigh",
                          help="number of genes in the heatmap [default= %default]", metavar="integer"),
    optparse::make_option(c("--norm_factor"), type="logical", default=FALSE, dest="norm_factor",
                          help="generate file with normalize factor for each condition/sample [default= %default]", metavar="logical"),
    optparse::make_option(c("--norm_counts"), type="logical", default=FALSE, dest="norm_counts",
                          help="Generate files with mormalized counts [default= %default]", metavar="logical"),
    optparse::make_option(c("--VD"), type="character", default=NULL, dest="VD",
                          help="Plot VennDiagram, precise type of comparison: all, down, up or both [default=%default]", metavar = "character"),
    optparse::make_option(c("--compaVD"), type="character", default=NULL, dest="compaVD",
                          help="Contrast comparison list to display in VennDiagram [default= %default]", metavar="character"),
    optparse::make_option(c("--GO"), type="character", default=NULL, dest="GO",
                          help="GO enrichment analysis for gene expressed 'up', 'down' or 'both', NULL for no GO enrichment. [default= %default]", metavar="character"),
    optparse::make_option(c("--ID2GO"), type="character", default=NULL, dest="geneID2GO_file",
                          help="GO annotation file [default= %default]", metavar="character"),
    optparse::make_option(c("--GO_threshold"), type="numeric", default="0.05", dest="GO_threshold",
                          help="the significant threshold used to filter p-values [default=%default]", metavar="double"),
    optparse::make_option(c("--GO_min_num_genes"), type="integer", default="10", dest="GO_min_num_genes",
                          help="the minimum number of genes for each GO terms [default=%default]", metavar="integer"),
    optparse::make_option(c("--GO_min_sig_genes"), type="integer", default="0", dest="GO_min_sig_genes",
                          help="the minimum number of significant gene(s) behind the enriched GO-term [default=%default]", metavar="integer"),
    optparse::make_option(c("--GO_max_top_terms"), type="integer", default="10", dest="GO_max_top_terms",
                          help="the maximum number of GO terms plot [default=%default]", metavar="integer"),
    optparse::make_option(c("--GO_algo"), type="character", default="weight01", dest="GO_algo",
                          help="algorithms which are accessible via the runTest function: shown by the whichAlgorithms() function, [default=%default]", metavar="character"),
    optparse::make_option(c("--GO_stats"), type="character", default="fisher", dest="GO_stats",
                          help = "statistical tests which are accessible via the runTest function: shown by the whichTests() function, [default=%default]", metavar = "character"),
    optparse::make_option(c("--Ratio_threshold"), type="numeric", default="0", dest="Ratio_threshold",
                          help="the minimum ratio value to display GO in graph [default=%default]", metavar="double"),
    optparse::make_option(c("--plotMD"),type="logical", default=FALSE, dest="plotMD", metavar="logical",
                          help="Mean-Difference Plot of Expression Data (aka MA plot) [default= %default]"),
    optparse::make_option(c("--plotVO"),type="logical", default=FALSE, dest="plotVO", metavar="logical",
                          help="Volcano plot for a specified coefficient/contrast of a linear model [default= %default]"),
    optparse::make_option(c("--glimMD"),type="logical", default=FALSE, dest="glimMD", metavar="logical",
                          help="Glimma - Interactif Mean-Difference Plot of Expression Data (aka MA plot) [default= %default]"),
    optparse::make_option(c("--glimVO"),type="logical", default=FALSE, dest="glimVO", metavar="logical",
                          help="Glimma - Interactif Volcano plot for a specified coefficient/contrast of a linear model [default= %default]"),
    optparse::make_option(c("--dens_bottom_mar"), type="integer", default="20", dest="densbotmar", metavar="integer",
                          help="Set bottom margin of density plot to help position the legend [default= %default]"),
    optparse::make_option(c("--dens_inset"), type="double", default="0.45", dest="densinset", metavar="double",
                          help="Set position the legend in bottom density graphe [default= %default]"),
    optparse::make_option(c("--legend_col"), type="integer", default="6", dest="legendcol", metavar="integer",
                          help="Set numbers of column for density plot legends [default= %default]"),
    optparse::make_option(c("--upset_basic"),type="character", default=NULL, dest="upset_basic",
                          help="Display UpSetR charts for all contrasts, precise type of comparison: all, down, up, mixed [default=%default].", metavar = "character"),
    optparse::make_option(c("--upset_type"),type="character", default=NULL, dest="upset_type",
                          help="Display UpSetR charts for list of contrasts, precise type of comparison: all, down, up, mixed [default=%default].", metavar = "character"),
    optparse::make_option(c("--upset_list"), type="character", default=NULL, dest="upset_list",
                          help="Contrast comparison list to display in UpSetR chart. See documentation. [default=%default]", metavar="character"),
    optparse::make_option("--coseq_data", type="character", default="ExpressionProfiles", dest="coseq_data",
                          help="Coseq data (ExpressionProfiles, LogScaledData) [default= %default]", metavar="character"),
    optparse::make_option("--coseq_model", type="character", default="kmeans", dest="coseq_model",
                          help="Coseq model (kmeans, Normal) [default= %default]", metavar="character"),
    optparse::make_option("--coseq_transformation", type="character", default="clr", dest="coseq_transformation",
                          help="Coseq tranformation (voom, logRPKM, arcsin, logit, logMedianRef, logclr, clr, alr, ilr, none) [default= %default]", metavar="character"),
    optparse::make_option("--coseq_ClustersNb", type="double", default=2:25, dest="coseq_ClustersNb",
                          help="Coseq : number of clusters desired (2:25 (auto), number from 2 to 25) [default= %default]", metavar="double"),
    optparse::make_option("--coseq_HeatmapOrderSample", type="logical", default=FALSE, dest="coseq_HeatmapOrderSample",
                          help="Choose TRUE if you prefer keeping your sample order than clusterizing samples in heatmap  [default= %default]", metavar="logical")
  )
  # Get command line options
  opt_parser = optparse::OptionParser(option_list=option_list)
  parameters = optparse::parse_args(opt_parser)

  if(is.null(parameters$rm_sample) == FALSE ) {
    stringr::str_replace_all(parameters$rm_sample, " ", "")
    parameters$rm_sample<-limma::strsplit2(parameters$rm_sample, ",")
  }

  if(is.null(parameters$select_sample) == FALSE ) {
    stringr::str_replace_all(parameters$select_sample, " ", "")
    parameters$select_sample<-limma::strsplit2(parameters$select_sample, ",")
  }

  return(parameters)
}

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
#' @param parameters, list that contains all arguments charged in Asko_start
#' @return data, list contain all data and metadata (DGEList, samples descriptions, contrast, design and annotations)
#'
#' @examples
#' \dontrun{
#'     parameters<-Asko_start()
#'     data<-loadData(parameters)
#' }
#'
#' @export
loadData <- function(parameters){

  # Folders for output files
  #---------------------------------------------------------
  cat("\nCreated directories:\n")
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  if(dir.exists(study_dir)==FALSE){ dir.create(study_dir) }
  cat("\t",study_dir,"\n")

  explo_dir = paste0(study_dir,"DataExplore/")
  if(dir.exists(explo_dir)==FALSE){ dir.create(explo_dir) }
  cat("\t",explo_dir,"\n")

  de_dir = paste0(study_dir,"DEanalysis/")
  if(dir.exists(de_dir)==FALSE){ dir.create(de_dir) }
  cat("\t",de_dir,"\n")

  image_dir = paste0(de_dir, "Images/")
  if(dir.exists(image_dir)==FALSE){ dir.create(image_dir) }
  cat("\t",image_dir,"\n")

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
  #---------------------------------------------------------
  input_path = "/import/"

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
  #----------------------------------------------------
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
  #---------------------------------------------------------
  Group<-factor(samples$condition)
  cat("\nConditions :\n")
  print(Group)
  designExp<-stats::model.matrix(~0+Group)
  rownames(designExp) <- row.names(samples)
  colnames(designExp) <- levels(Group)

  # Contrast for DE analysis
  #---------------------------------------------------------
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

#' @title asko3c
#'
#' @description Create contrast/condition/context file in format readable by Askomics Software.
#'
#' @param data_list, list contain all data and metadata (DGEList, samples descriptions, contrast, design and annotations)
#' @param parameters, list that contains all arguments charged in Asko_start
#' @return asko, list of data.frame contain condition, contrast and context information
#'
#' @examples
#' \dontrun{
#'     parameters <- Asko_start()
#'     data<-loadData(parameters)
#'     asko<-asko3c(data, parameters)
#' }
#'
#' @export
asko3c <- function(data_list, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")

  # Askomics directory
  asko_dir = paste0(study_dir, "DEanalysis/AskoTables/")
  asko<-list()

  # Condition
  #---------------
  condition<-unique(data_list$samples$condition)                       # retrieval of different condition's names
  col1<-which(colnames(data_list$samples)=="condition")                # determination of number of the column "condition"
  colcol<-which(colnames(data_list$samples)=="color")
  if(is.null(parameters$fileofcount)){
    col2<-which(colnames(data_list$samples)=="file")                   # determination of number of the column "replicate"
    column_name<-colnames(data_list$samples[,c(-col1,-col2,-colcol)])  # retrieval of column names needful to create the file condition
  }else{column_name<-colnames(data_list$samples[,c(-col1,-colcol)])}
  condition_asko<-data.frame(row.names=condition)                      # initialization of the condition's data frame

  for (name in column_name){                                           # for each experimental factor :
    condition_asko$n<-NA                                               # initialization of new column in the condition's data frame
    colnames(condition_asko)[colnames(condition_asko)=="n"]<-name      # to rename the new column with the name of experimental factor
    for(condition_name in condition){                                  # for each condition's names
      condition_asko[condition_name,name]<-as.character(unique(data_list$samples[data_list$samples$condition==condition_name, name]))
    }                                                                  # filling the condition's data frame
  }

  # Contrast + Context
  #--------------------------
  i=0
  contrast_asko<-data.frame(row.names = colnames(data_list$contrast))           # initialization of the contrast's data frame
  contrast_asko$Contrast<-NA                                                    # all columns are created et initialized with
  contrast_asko$context1<-NA                                                    # NA values
  contrast_asko$context2<-NA                                                    #

  list_context<-list()                                                          # initialization of context and condition lists
  list_condition<-list()                                                        # will be used to create the context data frame
  if(parameters$mk_context==TRUE){
    for (contrast in colnames(data_list$contrast)){                             # for each contrast :
      i=i+1                                                                       # contrast data frame will be filled line by line
      set_cond1<-row.names(data_list$contrast)[data_list$contrast[,contrast]>0]   # retrieval of 1st set of condition's names implicated in a given contrast
      set_cond2<-row.names(data_list$contrast)[data_list$contrast[,contrast]<0]   # retrieval of 2nd set of condition's names implicated in a given contrast
      set_condition<-colnames(condition_asko)                                     # retrieval of names of experimental factor

      if(length(set_cond1)==1){complex1=FALSE}else{complex1=TRUE}                        # to determine if we have complex contrast (multiple conditions
      if(length(set_cond2)==1){complex2=FALSE}else{complex2=TRUE}                        # compared to multiple conditions) or not
      if(complex1==FALSE && complex2==FALSE){                                             # Case 1: one condition against one condition
        contrast_asko[i,"context1"]<-set_cond1                                    # filling contrast data frame with the name of the 1st context
        contrast_asko[i,"context2"]<-set_cond2                                    # filling contrast data frame with the name of the 2nd context
        contrast_name<-paste(set_cond1,set_cond2, sep = "vs")                     # creation of contrast name by associating the names of contexts
        contrast_asko[i,"Contrast"]<-contrast_name                                # filling contrast data frame with contrast name
        list_context<-append(list_context, set_cond1)                             #
        list_condition<-append(list_condition, set_cond1)                         # adding respectively to the lists "context" and "condition" the context name
        list_context<-append(list_context, set_cond2)                             # and the condition name associated
        list_condition<-append(list_condition, set_cond2)                         #
      }
      if(complex1==FALSE && complex2==TRUE){                                             # Case 2: one condition against multiple condition
        contrast_asko[i,"context1"]<-set_cond1                                    # filling contrast data frame with the name of the 1st context
        list_context<-append(list_context, set_cond1)                             # adding respectively to the lists "context" and "condition" the 1st context
        list_condition<-append(list_condition, set_cond1)                         # name and the condition name associated
        l=0
        # "common_factor" will contain the common experimental factors shared by
        common_factor=list()                                                      # conditions belonging to the complex context
        for (param_names in set_condition){                                       # for each experimental factor
          facteur<-unique(c(condition_asko[,param_names]))                        # retrieval of possible values for the experimental factor
          l=l+1                                                                   #
          for(value in facteur){                                                  # for each possible values
            verif<-unique(stringr::str_detect(set_cond2, value))                           # verification of the presence of values in each condition contained in the set
            if(length(verif)==1 && verif==TRUE){common_factor[l]<-value}          # if verif contains only TRUE, value of experimental factor
          }                                                                       # is added as common factor
        }
        if(length(common_factor)>1){                                              # if there are several common factor
          common_factor<-toString(common_factor)                                  # the list is converted to string
          contx<-stringr::str_replace(common_factor,", ","")
          contx<-stringr::str_replace_all(contx, "NULL", "")}else{contx<-common_factor}    # and all common factor are concatenated to become the name of context
        contrast_asko[i,"context2"]<-contx                                        # filling contrast data frame with the name of the 2nd context
        contrast_name<-paste(set_cond1,contx, sep = "vs")                         # concatenation of context names to make the contrast name
        contrast_asko[i,"Contrast"]<-contrast_name                                # filling contrast data frame with the contrast name
        for(j in length(set_cond2)){                                            # for each condition contained in the complex context (2nd):
          list_context<-append(list_context, contx)                               # adding condition name with the context name associated
          list_condition<-append(list_condition, set_cond2[j])                    # to their respective list
        }
      }
      if(complex1==TRUE && complex2==FALSE){                                             # Case 3: multiple conditions against one condition
        contrast_asko[i,"context2"]<-set_cond2                                    # filling contrast data frame with the name of the 2nd context
        list_context<-append(list_context, set_cond2)                             # adding respectively to the lists "context" and "condition" the 2nd context
        list_condition<-append(list_condition, set_cond2)                         # name and the 2nd condition name associated
        l=0
        # "common_factor" will contain the common experimental factors shared by
        common_factor=list()                                                      # conditions belonging to the complex context
        for (param_names in set_condition){                                       # for each experimental factor:
          facteur<-unique(c(condition_asko[,param_names]))                        # retrieval of possible values for the experimental factor
          l=l+1
          for(value in facteur){                                                  # for each possible values:
            verif<-unique(stringr::str_detect(set_cond1, value))                           # verification of the presence of values in each condition contained in the set
            if(length(verif)==1 && verif==TRUE){common_factor[l]<-value}          # if verif contains only TRUE, value of experimental factor
          }                                                                       # is added as common factor
        }
        if(length(common_factor)>1){                                              # if there are several common factor
          common_factor<-toString(common_factor)                                  # the list is converted to string
          contx<-stringr::str_replace(common_factor,", ","")
          contx<-stringr::str_replace_all(contx, "NULL", "")}else{contx<-common_factor}    # and all common factor are concatenated to become the name of context
        contrast_asko[i,"context1"]<-contx                                        # filling contrast data frame with the name of the 1st context
        contrast_name<-paste(contx,set_cond2, sep = "vs")                         # concatenation of context names to make the contrast name
        contrast_asko[i,"Contrast"]<-contrast_name                                # filling contrast data frame with the contrast name
        for(j in length(set_cond1)){                                            # for each condition contained in the complex context (1st):
          list_context<-append(list_context, contx)                               # adding condition name with the context name associated
          list_condition<-append(list_condition, set_cond1[j])                    # to their respective list
        }
      }
      if(complex1==TRUE && complex2==TRUE){                                             # Case 4: multiple conditions against multiple conditions
        m=0                                                                       #
        n=0                                                                       #
        common_factor1=list()                                                     # list of common experimental factors shared by conditions of the 1st context
        common_factor2=list()                                                     # list of common experimental factors shared by conditions of the 2nd context
        w=1
        for (param_names in set_condition){                                       # for each experimental factor:
          print(w)
          w=w+1
          facteur<-unique(c(condition_asko[,param_names]))                        # retrieval of possible values for the experimental factor

          for(value in facteur){                                                  # for each possible values:
            verif1<-unique(stringr::str_detect(set_cond1, value))                          # verification of the presence of values in each condition contained in the 1st context
            verif2<-unique(stringr::str_detect(set_cond2, value))                          # verification of the presence of values in each condition contained in the 2nd context

            if(length(verif1)==1 && verif1==TRUE){m=m+1;common_factor1[m]<-value} # if verif=only TRUE, value of experimental factor is added as common factor
            if(length(verif2)==1 && verif2==TRUE){n=n+1;common_factor2[n]<-value} # if verif=only TRUE, value of experimental factor is added as common factor
          }
        }
        if(length(common_factor1)>1){                                             # if there are several common factor for conditions in the 1st context
          common_factor1<-toString(common_factor1)                                # conversion list to string
          contx1<-stringr::str_replace(common_factor1,", ","")}else{contx1<-common_factor1}# all common factor are concatenated to become the name of context
        contx1<-stringr::str_replace_all(contx1, "NULL", "")
        if(length(common_factor2)>1){                                             # if there are several common factor for conditions in the 2nd context
          common_factor2<-toString(common_factor2)                                # conversion list to string
          contx2<-stringr::str_replace(common_factor2,", ","")}else{contx2<-common_factor2}# all common factor are concatenated to become the name of context
        contx2<-stringr::str_replace_all(contx2, "NULL", "")
        contrast_asko[i,"context1"]<-contx1                                       # filling contrast data frame with the name of the 1st context
        contrast_asko[i,"context2"]<-contx2                                       # filling contrast data frame with the name of the 2nd context
        contrast_asko[i,"Contrast"]<-paste(contx1,contx2, sep = "vs")             # filling contrast data frame with the name of the contrast
        for(j in seq_along(set_cond1)){                                            # for each condition contained in the complex context (1st):
          list_context<-append(list_context, contx1)                              # verification of the presence of values in each condition
          list_condition<-append(list_condition, set_cond1[j])                    # contained in the 1st context
        }
        for(j in seq_along(set_cond2)){                                            # for each condition contained in the complex context (2nd):
          list_context<-append(list_context, contx2)                              # verification of the presence of values in each condition
          list_condition<-append(list_condition, set_cond2[j])                    # contained in the 1st context
        }
      }
    }
  }
  else{
    for (contrast in colnames(data_list$contrast)){
      i=i+1
      contexts=limma::strsplit2(contrast,"vs")
      contrast_asko[i,"Contrast"]<-contrast
      contrast_asko[i,"context1"]=contexts[1]
      contrast_asko[i,"context2"]=contexts[2]
      set_cond1<-row.names(data_list$contrast)[data_list$contrast[,contrast]>0]
      set_cond2<-row.names(data_list$contrast)[data_list$contrast[,contrast]<0]
      for (cond1 in set_cond1){
        list_context<-append(list_context, contexts[1])
        list_condition<-append(list_condition, cond1)
      }
      for (cond2 in set_cond2){
        list_context<-append(list_context, contexts[2])
        list_condition<-append(list_condition, cond2)
      }
    }
  }

  list_context<-unlist(list_context)
  list_condition<-unlist(list_condition)                                                                    # conversion list to vector
  context_asko<-data.frame(list_context,list_condition)                                                     # creation of the context data frame
  context_asko<-unique(context_asko)
  colnames(context_asko)[colnames(context_asko)=="list_context"]<-"context"                                 # header formatting for askomics
  colnames(context_asko)[colnames(context_asko)=="list_condition"]<-"condition"                             # header formatting for askomics
  asko<-list("condition"=condition_asko, "contrast"=contrast_asko, "context"=context_asko)                  # adding context data frame to asko object
  colnames(context_asko)[colnames(context_asko)=="context"]<-"Context"                                      # header formatting for askomics
  colnames(context_asko)[colnames(context_asko)=="condition"]<-"has@Condition"                              # header formatting for askomics
  colnames(contrast_asko)[colnames(contrast_asko)=="context1"]<-paste("context1_of", "Context", sep="@")    # header formatting for askomics
  colnames(contrast_asko)[colnames(contrast_asko)=="context2"]<-paste("context2_of", "Context", sep="@")    # header formatting for askomics

  # Files creation
  #-------------------
  # creation of condition file for asko
  utils::write.table(data.frame("Condition"=row.names(condition_asko),condition_asko),
                     paste0(asko_dir,"condition.asko.txt"),
                     sep = parameters$sep,
                     row.names = FALSE,
                     quote=FALSE)
  # creation of context file for asko
  utils::write.table(context_asko,
                     paste0(asko_dir, "context.asko.txt"),
                     sep=parameters$sep,
                     col.names = TRUE,
                     row.names = FALSE,
                     quote=FALSE)
  # creation of contrast file for asko
  utils::write.table(contrast_asko,
                     paste0(asko_dir, "contrast.asko.txt"),
                     sep=parameters$sep,
                     col.names = TRUE,
                     row.names = FALSE,
                     quote=FALSE)
  return(asko)
}

#' @title GEfilt
#'
#' @description
#' \itemize{
#'    \item Filter genes according to cpm threshold value and replicated cpm option value.
#'    \item Plot different graphes to explore data before and after filtering.
#' }
#'
#' @param data_list, list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations)
#' @param parameters, list that contains all arguments charged in Asko_start
#' @return filtered_counts, large DGEList with filtered counts and data descriptions.
#'
#' @examples
#' \dontrun{
#'     filtered_counts<-GEfilt(data, parameters)
#' }
#'
#' @export
GEfilt <- function(data_list, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "DataExplore/")

  # plot density before filtering
  #---------------------------------
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

  grDevices::png(paste0(image_dir, parameters$analysis_name, "_raw_data.png"), width=1024, height=1024)
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
  #---------------------------------
  keep.exprs <- rowSums(cpm>parameters$threshold_cpm)>=parameters$replicate_cpm
  filtered_counts <- data_list$dge[keep.exprs,,keep.lib.sizes=FALSE]
  filtered_cpm<-edgeR::cpm(filtered_counts$counts, log=TRUE)

  maxi<-c()
  for (i in seq(nsamples)){
    m=max(stats::density(filtered_cpm[,i])$y)
    maxi<-c(maxi,m)
  }
  ymax<-round(max(maxi),1) + 0.02

  grDevices::png(paste0(image_dir,parameters$analysis_name,"_filtered_data.png"), width=1024, height=1024)
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
  #------------------------------------------------------
  grDevices::png(paste0(image_dir,parameters$analysis_name,"_barplot_logcpm_before_filtering.png"), width=sizeImg, height=sizeImg)
  graphics::hist(logcpm,
                 main= "A. Log2(cpm) distribution before filtering",
                 xlab = "log2(cpm)",
                 col = "grey")
  grDevices::dev.off()

  # histogram cpm values distribution after filtering
  #------------------------------------------------------
  grDevices::png(paste0(image_dir,parameters$analysis_name,"_barplot_logcpm_after_filtering.png"), width=sizeImg, height=sizeImg)
  graphics::hist(filtered_cpm,
                 main= "B. Log2(cpm) distribution after filtering",
                 xlab = "log2(cpm)",
                 col = "grey")
  grDevices::dev.off()

  # boxplot cpm values distribution before filtering
  #------------------------------------------------------
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
  #------------------------------------------------------
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
  #------------------------------------------------------
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
  #------------------------------------------------------
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

#' @title GEnorm
#'
#' @description Normalize counts
#' \itemize{
#'    \item Calculate normalization factors to scale the filtered library sizes.
#'    \item Plot different graphes to explore data before and after normalization.
#'    \item Optionally, write file with mean counts and normalized mean counts in Askomics format.
#' }
#'
#' @param filtered_GE, large DGEList with filtered counts by GEfilt function.
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @param asko_list, list of data.frame contain condition, contrast and context informations made by asko3c.
#' @param data_list, list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations)
#' @return norm_GE, large DGEList with normalized counts and data descriptions.
#'
#' @examples
#' \dontrun{
#'     norm_GE<-GEnorm(filtered_counts, asko, parameters)
#' }
#'
#' @export
GEnorm <- function(filtered_GE, asko_list, data_list, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "DataExplore/")

  # for image size
  nsamples <- ncol(filtered_GE$counts)
  sizeImg=15*nsamples
  if(sizeImg < 480){ sizeImg=480 }

  # Normalization counts
  norm_GE<-edgeR::calcNormFactors(filtered_GE, method = parameters$normal_method)
  if(parameters$norm_factor == TRUE){
    utils::write.table(norm_GE$samples, file=paste0(study_dir, parameters$analysis_name, "_normalization_factors.txt"), col.names=NA, row.names=TRUE, quote=FALSE, sep="\t", dec=".")
  }

  # boxplot log2(cpm) values after normalization
  #----------------------------------------------------
  logcpm_norm <- edgeR::cpm(norm_GE, log=TRUE)
  colnames(logcpm_norm)<-rownames(filtered_GE$samples)

  grDevices::png(paste0(image_dir,parameters$analysis_name,"_boxplot_logcpm_after_norm.png"), width=sizeImg, height=sizeImg)
  graphics::par(oma=c(1,1,1,1))
  graphics::boxplot(logcpm_norm,
                    col=filtered_GE$samples$color,
                    main="C. Log2(cpm) distribution after normalization",
                    cex.axis=0.8,
                    las=2,
                    ylab="Log2(cpm)")
  grDevices::dev.off()

  dScaleFactors <- norm_GE$samples$lib.size * norm_GE$samples$norm.factors
  normCounts <- t(t(filtered_GE$counts)/dScaleFactors)*mean(dScaleFactors)

  # File with normalized counts
  if (parameters$norm_counts == TRUE){
    utils::write.table(normCounts, file=paste0(study_dir, parameters$analysis_name, "_NormCounts.txt"), col.names=NA, row.names=TRUE, quote=FALSE, sep="\t", dec=".")
  }


  # barplot counts after normalization
  #------------------------------------------------------
  grDevices::png(paste0(image_dir,parameters$analysis_name,"_barplot_SumCounts_after_norm.png"), width=sizeImg, height=sizeImg)
  graphics::par(oma=c(1,1,1,1),mar=c(7, 7, 4, 2), mgp=c(5,0.7,0))
  graphics::barplot(colSums(normCounts),
                    col=data_list$dge$samples$color,
                    main=paste0("C. Sum of normalized counts from ",nrow(normCounts)," transcripts"),
                    cex.axis=0.8,
                    las=2,
                    ylab="Counts sum",
                    xlab="Samples")
  grDevices::dev.off()

  # heatmap visualisation
  #----------------------------------------------------
  if(parameters$CompleteHeatmap==TRUE)
  {
    # heatmap cpm value per sample
    #----------------------------------------------------
    cpm_norm  <- edgeR::cpm(norm_GE, log=FALSE)
    cpmscale  <- scale(t(cpm_norm))
    tcpmscale <- t(cpmscale)

    d1 <- stats::dist(cpmscale,  method = parameters$distcluts, diag = FALSE, upper = FALSE)
    d2 <- stats::dist(tcpmscale, method = parameters$distcluts, diag = FALSE, upper = TRUE)
    hc <- stats::hclust(d1, method = parameters$hclust, members = NULL)
    hr <- stats::hclust(d2, method = parameters$hclust, members = NULL)
    my_palette <- grDevices::colorRampPalette(c("green","black","red"), interpolate = "linear")

    grDevices::png(paste0(image_dir,parameters$analysis_name,"_heatmap_CPMcounts_per_sample.png"), width=sizeImg*1.5, height=sizeImg*1.25)
    graphics::par(oma=c(2,1,2,2))
    gplots::heatmap.2(tcpmscale, Colv = stats::as.dendrogram(hc), Rowv = stats::as.dendrogram(hr), density.info="histogram",
                      trace = "none", dendrogram = "column", xlab = "samples", col = my_palette, labRow = FALSE,
                      cexRow = 0.1, cexCol = 1.25, ColSideColors = norm_GE$samples$color, margins = c(10,1),
                      main = paste0("CPM counts per sample\nGenes 1 to ",nrow(norm_GE)))
    grDevices::dev.off()

    # Normalized mean by conditions
    #-------------------------------
    # heatmap mean counts per condition
    n_count <- NormCountsMean(norm_GE, ASKOlist = asko_list)
    countscale  <- scale(t(n_count))
    tcountscale <- t(countscale)

    d1 <- stats::dist(countscale,  method = parameters$distcluts, diag = FALSE, upper = FALSE)
    d2 <- stats::dist(tcountscale, method = parameters$distcluts, diag = FALSE, upper = TRUE)
    hc <- stats::hclust(d1, method = parameters$hclust, members = NULL)
    hr <- stats::hclust(d2, method = parameters$hclust, members = NULL)
    my_palette <- grDevices::colorRampPalette(c("green","black","red"), interpolate = "linear")

    grDevices::png(paste0(image_dir,parameters$analysis_name,"_heatmap_CPMmean_per_condi.png"), width=sizeImg*1.5, height=sizeImg*1.25)
    graphics::par(oma=c(2,1,2,2))
    gplots::heatmap.2(tcountscale, Colv = stats::as.dendrogram(hc), Rowv = stats::as.dendrogram(hr), density.info="histogram",
                      trace = "none", dendrogram = "column", xlab = "Condition", col = my_palette, labRow = FALSE,
                      cexRow = 0.1, cexCol = 1.5, ColSideColors = unique(norm_GE$samples$color), margins = c(10,1),
                      main = paste0("CPM counts per condition (mean)\nGenes 1 to ",nrow(norm_GE)))
    grDevices::dev.off()
  }

  # File with normalized counts in CPM by sample
  cpm_norm <- edgeR::cpm(norm_GE, log=FALSE)
  utils::write.table(cpm_norm, file=paste0(study_dir, parameters$analysis_name, "_CPM_NormCounts.txt"), col.names=NA, row.names=TRUE, quote=FALSE, sep="\t", dec=".")
  # File with normalized mean counts in CPM, grouped by condition
  tempo<-as.data.frame(t(stats::aggregate(t(cpm_norm),list(data_list$samples$condition), mean)))
  colnames(tempo)<-tempo["Group.1",]
  meancpmDEGnorm<-tempo[-1,]
  utils::write.table(meancpmDEGnorm,paste0(study_dir, parameters$analysis_name,"_CPM_NormMeanCounts.txt"), sep="\t", dec=".", row.names=TRUE, col.names=NA, quote=FALSE)

  return(norm_GE)
}

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
#' @param asko_norm, large DGEList with normalized counts by GEnorm function.
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @return none
#'
#' @import ggfortify
#'
#' @examples
#' \dontrun{
#'     GEcorr(asko_norm,parameters)
#' }
#'
#' @export
GEcorr <- function(asko_norm, parameters){
  #library(ggfortify)
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
  #-----------------------------
  cormat<-stats::cor(lcpm)
  color<-grDevices::colorRampPalette(c("black","red","yellow","white"),space="rgb")(35)
  grDevices::png(paste0(image_dir, parameters$analysis_name, "_heatmap_of_sample_correlation.png"), width=sizeImg, height=sizeImg)
  graphics::par(oma=c(4,2,4,1))
  stats::heatmap(cormat, col=color, symm=TRUE, RowSideColors=as.character(asko_norm$samples$color),
                 ColSideColors=as.character(asko_norm$samples$color), main="")
  graphics::title("Sample Correlation Matrix", adj=0.5, outer=TRUE)
  grDevices::dev.off()

  corr<-stats::cor(edgeR::cpm(asko_norm, log=FALSE))
  minCorr=min(corr)
  maxCorr=max(corr)
  grDevices::png(paste0(image_dir, parameters$analysis_name, "_correlogram.png"), width=sizeImg, height=sizeImg)
  corrplot(corr, method="ellipse", type = "lower", tl.col = "black", tl.srt = 45, is.corr = FALSE, cl.lim=c(minCorr,maxCorr))
  #corrplot.mixed(corr, lower = "number", upper = "ellipse", tl.col = "black",is.corr = FALSE, cl.lim=c(minCorr,maxCorr))
  graphics::title("Sample Correlogram", adj=0.5)
  grDevices::dev.off()

  # MDS Plot
  #-----------------------------
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
  #-----------------------------
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
  #-----------------------------
  mat.dist <- stats::dist(t(asko_norm$counts), method = parameters$distcluts)
  clustering <- stats::hclust(mat.dist, method=parameters$hclust)
  grDevices::png(paste0(image_dir, parameters$analysis_name, "_hclust.png"), width=sizeImg, height=sizeImg)
  graphics::par(oma=c(1,1,1,1))
  plot(clustering,
       main = 'Distances Correlation\nHierarchical clustering', sub="",
       xlab=paste0("Samples\nmethod hclust: ",parameters$hclust),
       hang = -1)
  j<-grDevices::dev.off()
}

#' @title plot_glimma
#'
#' @description
#' Use Glimma package for interactive visualization of results from differential
#' expression analyses. Two types of graphs can be created:
#' \itemize{
#'    \item Mean-Difference Plot of Expression Data (aka MA plot).
#'    \item Volcano plot for a specified coefficient/contrast of a linear model.
#' }
#'
#' @param fit, fitted linear model object.
#' @param normGE, large DGEList with normalized counts and data description.
#' @param resDE, vector containing integer values of -1 to represent down-regulated
#' genes, 0 for no differential expression, and 1 for up-regulated genes.
#' @param contrast, coefficient/contrast tested.
#' @param tplot, type of plot selected for display.
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @return none.
#'
#' @examples
#' \dontrun{
#'    plot_glimma(fit, normGE, resDE, contrast, "MD", parameters)  # smear plot
#'    plot_glimma(fit, normGE, resDE, contrast, "VO", parameters)  # volcano plot
#' }
#'
#' @export
plot_glimma <- function(fit, normGE, resDE, contrast, tplot, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "DEanalysis/Images/")

  # Mean-Difference Plot
  #--------------------------
  if(tplot=="MD"){
    if (is.null(normGE$samples$color)==TRUE){
      suppressWarnings(Glimma::glMDPlot(fit, status=resDE[,contrast], counts=normGE, group=normGE$samples$condition,
                                        transform=TRUE, anno=NULL, launch=FALSE, main=contrast,
                                        folder=paste0(image_dir, "Glimma_Plots"), html=paste0("MDPlot_",contrast)))
    }
    else{
      suppressWarnings(Glimma::glMDPlot(fit, status=resDE[,contrast], counts=normGE, group=normGE$samples$condition, transform=TRUE,
                                        sample.cols=normGE$samples$color, anno=NULL, launch=FALSE, main=contrast,
                                        folder=paste0(image_dir, "Glimma_Plots"), html=paste0("MDPlot_",contrast)))
    }
  }

  # Volcano plot
  #--------------------------
  if (tplot=="VO"){
    if (is.null(normGE$samples$color)==TRUE){
      Glimma::glXYPlot(x=fit$table$logFC, y=-log10(fit$table$PValue), status=resDE[,contrast], counts=normGE,
                       group=normGE$samples$condition, xlab="Log2FoldChange", ylab="-log10(pvalue)", main=contrast,
                       launch=FALSE, folder=paste0(image_dir, "Glimma_Plots"), html=paste0("Volcano_",contrast))
    }
    else{
      Glimma::glXYPlot(x=fit$table$logFC, y=-log10(fit$table$PValue), status=resDE[,contrast], counts=normGE, main=contrast,
                       group=normGE$samples$condition, xlab="Log2FoldChange", ylab="-log10(pvalue)", launch=FALSE,
                       sample.cols=normGE$samples$color, folder=paste0(image_dir, "Glimma_Plots"), html=paste0("Volcano_",contrast))
    }
  }
}

#' @title plot_expr
#'
#' @description
#' Function to generate plots showing different aspects of differential expression results.
#' Two types of graphs can be created:
#' \itemize{
#'    \item Mean-Difference Plot of Expression Data (aka MA plot).
#'    \item Volcano plot for a specified coefficient/contrast of a linear model.
#' }
#'
#' @param fit, fitted linear model object.
#' @param normGE, large DGEList with normalized counts and data description.
#' @param resDE, vector containing integer values of -1 to represent down-regulated
#' genes, 0 for no differential expression, and 1 for up-regulated genes.
#' @param contrast, coefficient/contrast tested.
#' @param tplot, type of plot selected for display.
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @return none.
#'
#' @examples
#' \dontrun{
#'    plot_expr(fit, normGE, resDE, contrast, "MD", parameters)  # smear plot
#'    plot_expr(fit, normGE, resDE, contrast, "VO", parameters)  # volcano plot
#' }
#'
#' @export
plot_expr <- function(fit, normGE, resDE, contrast, tplot, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "DEanalysis/Images/")



  # Mean-Difference Plot
  if(tplot=="MD"){
    grDevices::png(paste0(image_dir, contrast, "_MeanDifference_of_ExpressionData.png"))
    edgeR::plotMD.DGELRT(fit, xlab="Average log CPM", ylab="log-fold-change", main=paste0("MD plot - ", contrast),
                         cex=0.5, status=resDE[,contrast], values=c("-1","1"), col=c("blue","red"))
    grDevices::dev.off()
  }

  # Volcano plot
  if(tplot=="VO"){
    tglm<-fit$table
    tglm$FDR<-stats::p.adjust(tglm$PValue, method=parameters$p_adj_method)
    grDevices::png(paste0(image_dir, contrast, "_VolcanoPlot.png"))
    with(tglm, plot(tglm$logFC, -log10(tglm$PValue), pch=16, cex=0.5, xlim=c(min(tglm$logFC)-0.5, max(tglm$logFC)+0.5),
                    ylim=c(min(-log10(tglm$PValue))-0.5, max(-log10(tglm$PValue))+0.5),
                    main=paste0("Volcano plot - ", contrast), xlab="Log2FoldChange", ylab="-log10(pvalue)"))
    with(subset(tglm, tglm$FDR <= parameters$threshold_FDR & tglm$logFC >  parameters$threshold_logFC), graphics::points(logFC, -log10(PValue), pch=16, cex=0.5, col="red"))
    with(subset(tglm, tglm$FDR <= parameters$threshold_FDR & tglm$logFC < -parameters$threshold_logFC), graphics::points(logFC, -log10(PValue), pch=16, cex=0.5, col="blue"))
    grDevices::dev.off()
  }
}

#' @title NormCountsMean
#'
#' @description Calculation mean counts for two contrast or all matrix.
#'
#' @param glmfit, fitted linear model object.
#' @param ASKOlist, list of data.frame contain condition, contrast and context information made by asko3c.
#' @param context, coefficient/contrast tested.
#' @return one of this:
#' \itemize{
#'    \item matrixMean, matrix with mean counts,
#'    \item meanValue for one context/Condition.
#' }
#'
#' @examples
#' \dontrun{
#'     # calculate mean counts in contrast contx1_vs_contx2
#'     mean1<-NormCountsMean(glmfit, ASKOlist, contx1)   # in the 1st context
#'     mean2<-NormCountsMean(glmfit, ASKOlist, contx2)   # in the 2nd context
#'
#'     # for all conditions
#'     n_count<-NormCountsMean(glmfit, ASKOlist, context=NULL)
#' }
#'
#' @export
NormCountsMean <- function(glmfit, ASKOlist, context=NULL){
  lib_size_norm<-glmfit$samples$lib.size*glmfit$samples$norm.factors                          # normalization computation of all library sizes
  if(is.null(context)==TRUE){
    set_condi<-row.names(ASKOlist$condition)
  }else{
    set_condi<-ASKOlist$context$condition[ASKOlist$context$context==context]                  # retrieval of condition names associated to context
  }
  table_c_norm <- data.frame(row.names = row.names(glmfit$counts))

  for (condition in set_condi){
    sample_name<-rownames(glmfit$samples[glmfit$samples$condition==condition,])               # retrieval of the replicate names associated to conditions
    subset_counts<-data.frame(row.names = row.names(glmfit$counts))                           # initialization of data frame as subset of counts table
    for(name in sample_name){
      lib_sample_norm<-glmfit$samples[name,"lib.size"]*glmfit$samples[name,"norm.factors"]    # normalization computation of sample library size
      subset_counts$c<-glmfit$counts[,name]                                                   # addition in subset of sample counts column
      subset_counts$c<-subset_counts$c*mean(lib_size_norm)/lib_sample_norm                    # normalization computation of sample counts
      colnames(subset_counts)[colnames(subset_counts)=="c"]<-name                             # to rename the column with the condition name
    }

    mean_counts<-rowSums(subset_counts)/ncol(subset_counts)                                   # computation of the mean
    table_c_norm$m <- mean_counts

    if(is.null(context)==TRUE){
      colnames(table_c_norm)[colnames(table_c_norm)=="m"]<-condition
    }else{
      colnames(table_c_norm)[colnames(table_c_norm)=="m"]<-paste(context,condition,sep = "/")
    }
  }

  return(table_c_norm)
}

#' @title AskoStats
#'
#' @description Based on result contained in "glm_test":
#' \itemize{
#'    \item print summary result of differential expression analysis
#'    \item format all results in tabulate out file followed parameters given
#'    \item plot heatmap of top differential expressed genes
#' }
#' Create one file by contrast, each file contains for each genes: fold-change and
#' log fold-change values, PValue, Expression, Significance, logCPM, LR, FDR and
#' significance value for each condition/context.
#' By default, LR and logCPM were not displayed, you can switch this parametres
#' to TRUE for display.
#'
#' @param glm_test, tests for one or more coefficients in the linear model (likelihood ratio tests or empirical Bayes quasi-likelihood F-tests).
#' @param fit, fitted linear model object.
#' @param contrast, coefficient/contrast names tested.
#' @param ASKOlist, list of data.frame contain condition, contrast and context informations made by asko3c.
#' @param dge, large DGEList with normalized counts by GEnorm function.
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @return none
#'
#' @examples
#' \dontrun{
#'     AskoStats(glm_test, fit, contrast, ASKOlist, dge, parameters)
#' }
#'
#' @export
AskoStats <- function (glm_test, fit, contrast, ASKOlist, dge, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  asko_dir = paste0(study_dir, "DEanalysis/AskoTables/")
  image_dir = paste0(study_dir, "DEanalysis/Images/")

  contrasko<-ASKOlist$contrast$Contrast[row.names(ASKOlist$contrast)==contrast]   # to retrieve the name of contrast from Asko object
  contx1<-ASKOlist$contrast$context1[row.names(ASKOlist$contrast)==contrast]      # to retrieve the name of 1st context from Asko object
  contx2<-ASKOlist$contrast$context2[row.names(ASKOlist$contrast)==contrast]      # to retrieve the name of 2nd context from Asko object

  ASKO_stat<-glm_test$table
  ASKO_stat$Test_id<-paste(contrasko, rownames(ASKO_stat), sep = "_")             # addition of Test_id column = unique ID
  ASKO_stat$contrast<-contrasko                                                   # addition of the contrast of the test
  ASKO_stat$gene <- row.names(ASKO_stat)                                          # addition of gene column = gene ID
  ASKO_stat$FDR<-stats::p.adjust(ASKO_stat$PValue, method=parameters$p_adj_method)       # computation of False Discovery Rate

  # Between context1 and context2 :
  ASKO_stat$Significance=0
  ASKO_stat$Significance[ASKO_stat$logFC <= -parameters$threshold_logFC & ASKO_stat$FDR <= parameters$threshold_FDR] = -1  # Significance values = -1 for down regulated genes
  ASKO_stat$Significance[ASKO_stat$logFC >= parameters$threshold_logFC  & ASKO_stat$FDR <= parameters$threshold_FDR] = 1   # Significance values =  1 for up regulated genes

  # addition of column "expression"
  ASKO_stat$Expression=NA
  ASKO_stat$Expression[ASKO_stat$Significance==-1]<-paste(contx1, contx2, sep="<")  # the value of attribute "Expression" is a string
  ASKO_stat$Expression[ASKO_stat$Significance==1]<-paste(contx1, contx2, sep=">")   # this attribute is easier to read the Significance
  ASKO_stat$Expression[ASKO_stat$Significance==0]<-paste(contx1, contx2, sep="=")   # of expression between two contexts

  if(parameters$Expression==TRUE){colg="Expression"}else{colg=NULL}
  if(parameters$logFC==TRUE){cola="logFC"}else{cola=NULL}
  # computation of Fold Change from log2FC
  if(parameters$FC==TRUE){colb="FC";ASKO_stat$FC <- 2^abs(ASKO_stat$logFC)}else{colb=NULL}
  if(parameters$Sign==TRUE){colc="Significance"}
  if(parameters$logCPM==TRUE){cold="logCPM"}else{cold=NULL}
  if(parameters$LR==TRUE){cole="LR"}else{cole=NULL}
  if(parameters$FDR==TRUE){colf="FDR"}else{colf=NULL}

  # adding table "stat.table" to the ASKOlist
  ASKOlist$stat.table<-ASKO_stat[,c("Test_id","contrast","gene",cola,colb,"PValue",colg,colc,cold,cole,colf)]

  if(parameters$mean_counts==TRUE){                         # computation of the mean of normalized counts for conditions
    mean1<-NormCountsMean(fit, ASKOlist, contx1)            # in the 1st context
    mean2<-NormCountsMean(fit, ASKOlist, contx2)            # in the 2nd context
    ASKOlist$stat.table<-cbind(ASKOlist$stat.table, mean1)
    ASKOlist$stat.table<-cbind(ASKOlist$stat.table, mean2)
  }
  print(table(ASKO_stat$Expression))
  colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="gene"] <- paste("is", "gene", sep="@")                  # header formatting for askomics
  colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="contrast"] <- paste("measured_in", "Contrast", sep="@") # header formatting for askomics
  o <- order(ASKOlist$stat.table$FDR)                                                                                   # ordering genes by FDR value
  ASKOlist$stat.table<-ASKOlist$stat.table[o,]

  utils::write.table(ASKOlist$stat.table,paste0(asko_dir, parameters$organism, "_", contrasko, ".txt"), sep=parameters$sep, col.names = TRUE, row.names = FALSE, quote=FALSE)

  # for image size
  nsamples<-ncol(dge$counts)
  sizeImg=15*nsamples
  if(sizeImg < 480) {sizeImg=480}

  # heatmap of Most Differential Genes Expression
  if(parameters$heatmap==TRUE){
    cpm_gstats<-edgeR::cpm(dge, log=TRUE)[o,][seq(parameters$numhigh),]
    grDevices::png(paste0(image_dir, contrast, "_topDGE_heatmap.png"), width=sizeImg*1.5, height=sizeImg*1.5)
    graphics::par(oma=c(2,2,2,2))
    gplots::heatmap.2(cpm_gstats,
                      trace="none",
                      scale="row",
                      labCol=dge$samples$Name,
                      main = contrasko,
                      xlab = "samples",
                      ColSideColors = dge$samples$color,
                      margins = c(12,12),
                      Rowv = FALSE,
                      dendrogram="col")
    grDevices::dev.off()
  }
}

#' @title DEanalysis
#'
#' @description Genewise statistical tests for a given coefficient or contrast, with edgeR method.
#'
#' @param norm_GE, large DGEList with normalized counts and data description.
#' @param data_list, list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations).
#' @param asko_list, list of data.frame contain condition, contrast and context informations made by asko3c.
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @return SumMat, list (TestResults format class limma) contains for each contrast the significance expression (1/0/-1) for all gene.
#'
#' @import edgeR
#' @import limma
#'
#' @examples
#' \dontrun{
#'     sum_table<-DEanalysis(norm_GE, data_list, asko_list, parameters)
#' }
#'
#' @export
DEanalysis <- function(norm_GE, data_list, asko_list, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "DEanalysis/Images/")

  # for image size
  nsamples <- ncol(data_list$dge$counts)
  sizeImg=15*nsamples
  if(sizeImg < 480){ sizeImg=480 }

  # Checks Contrasts
  if(is.null(parameters$select_sample) & is.null(parameters$rm_sample)){
    c1<-unique(data_list$samples$condition)
    len1<-length(c1)
    c2<-rownames(data_list$contrast)
    len2<-length(c2)
    if(len1 > len2){
      cat("\n\n")
      stop("Contrast files must be contain all conditions (in rows).\n\n")
    }
    if(len1 < len2){
      cat("\n\n")
      stop("Too much condtions in contrast file!\n\n")
    }
    if(length(setdiff(c1,c2)) > 0 && length(setdiff(c2,c1)) > 0){
      cat("\n\n")
      stop("Erronate or unknown conditions names in contrast file!\n\n")
    }
  }

  # Estimate Common, Trended and Tagwise Dispersion for Negative Binomial GLMs
  if(parameters$glm_disp==TRUE)
  {
    normGEdisp <- estimateGLMCommonDisp(norm_GE, data_list$design)
    normGEdisp <- estimateGLMTrendedDisp(normGEdisp, data_list$design)
    normGEdisp <- estimateGLMTagwiseDisp(normGEdisp, data_list$design)
  }
  # Estimate Common, Trended and Tagwise Negative Binomial dispersions by weighted likelihood empirical Bayes
  else
  {
    normGEdisp <- estimateDisp(norm_GE, data_list$design)
  }
  grDevices::png(paste0(image_dir, "Biological_coefficient_of_variation.png"), width=sizeImg, height=sizeImg)
  plotBCV(normGEdisp)
  grDevices::dev.off()

  # Genewise Negative Binomial Generalized Linear Models
  if(parameters$glm=="lrt"){
    fit <- glmFit(normGEdisp, data_list$design, robust = TRUE)
  }
  # Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests
  else if(parameters$glm=="qlf"){
    fit <- glmQLFit(normGEdisp, data_list$design, robust = TRUE)
    grDevices::png(paste0(image_dir, parameters$analysis_name, "_quasi-likelihood_dispersion.png"), width=sizeImg, height=sizeImg)
    plotQLDisp(fit)
    grDevices::dev.off()
  }

  # data frame combine all status genes results for summary file
  sum<-data.frame(row.names = rownames(fit))
  sum2=list()
  # if only one contrast ask
  if(length(data_list$contrast)==1){
    contrast<-makeContrasts(contrasts = data_list$contrast, levels = data_list$design)
    colnames(contrast)<-colnames(data_list$contrast)
    # likelihood ratio tests for one or more coefficients in the linear model.
    if(parameters$glm=="lrt"){
      glm_test<-glmLRT(fit, contrast=contrast)
    }
    # similar to glmLRT except that it replaces likelihood ratio tests with empirical Bayes quasi-likelihood F-tests
    if(parameters$glm=="qlf"){
      glm_test<-glmQLFTest(fit, contrast=contrast)
    }
    sum[,colnames(contrast)]<-decideTestsDGE(glm_test, adjust.method = parameters$p_adj_method, lfc=parameters$threshold_logFC, p.value=parameters$threshold_FDR)
    AskoStats(glm_test, fit, colnames(contrast), asko_list, normGEdisp, parameters)

    # display grahes (volcano or/and MD)
    if(parameters$plotMD==TRUE) { plot_expr(glm_test, normGEdisp, sum, colnames(contrast), "MD", parameters) }
    if(parameters$plotVO==TRUE) { plot_expr(glm_test, normGEdisp, sum, colnames(contrast), "VO", parameters) }
    if(parameters$glimMD==TRUE) { plot_glimma(glm_test, normGEdisp, sum, colnames(contrast), "MD", parameters) }
    if(parameters$glimVO==TRUE) { plot_glimma(glm_test, normGEdisp, sum, colnames(contrast), "VO", parameters) }
  }
  # for more than one contrast
  else{
    for (contrast in colnames(data_list$contrast)){
      # likelihood ratio tests for one or more coefficients in the linear model.
      if(parameters$glm=="lrt"){
        glm_test<-glmLRT(fit, contrast=data_list$contrast[,contrast])
      }
      # similar to glmLRT except that it replaces likelihood ratio tests with empirical Bayes quasi-likelihood F-tests
      if(parameters$glm=="qlf"){
        glm_test<-glmQLFTest(fit, contrast=data_list$contrast[,contrast])
      }
      sum[,contrast]<-decideTestsDGE(glm_test, adjust.method = parameters$p_adj_method, lfc=parameters$threshold_logFC, p.value=parameters$threshold_FDR)
      AskoStats(glm_test, fit, contrast, asko_list, normGEdisp, parameters)

      # display grahes (volcano or/and MD)
      if(parameters$plotMD==TRUE) { plot_expr(glm_test, normGEdisp, sum, contrast, "MD", parameters) }
      if(parameters$plotVO==TRUE) { plot_expr(glm_test, normGEdisp, sum, contrast, "VO", parameters) }
      if(parameters$glimMD==TRUE) { plot_glimma(glm_test, normGEdisp, sum, contrast, "MD", parameters) }
      if(parameters$glimVO==TRUE) { plot_glimma(glm_test, normGEdisp, sum, contrast, "VO", parameters) }
    }
  }

  # Create summary file with annotations (if available) and contrast value for each gene
  #---------------------------------------------------------------------------------------
  cat("\nCreate Summary file\n\n")
  sumFile<-paste0(study_dir,"DEanalysis/",parameters$analysis_name,"_summary_DE.txt")
  if(is.null(data_list$annot)==FALSE)
  {
    rnames<-row.names(sum)                        # get Genes DE names
    annDE<-as.matrix(data_list$annot[rnames,])    # get annotations for each genes DE
    rownames(annDE)<-rnames
    colnames(annDE)<-colnames(data_list$annot)
    SumMat<-cbind(sum,annDE)                      # merge the two matrix

    utils::write.table(SumMat, file=sumFile, col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')
  }
  else
  {
    utils::write.table(sum, file=sumFile, col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')
  }

  # reformate summary result table
  newMat <- as.data.frame(matrix(unlist(sum), nrow=nrow(sum)))
  rownames(newMat)<-rownames(sum)
  colnames(newMat)<-colnames(sum)

  return(newMat)
}

#' @title UpSetGraph
#'
#' @description Generate upsetR graphs.
#'
#' @param resDEG, data frame contains for each contrast the significance expression (1/0/-1) for all gene.
#' @param data_list, list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations).
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @return none
#'
#' @examples
#' \dontrun{
#'    UpSetGraph(sumDEG, data_list, parameters)
#' }
#'
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
  #---------------------------------------------------------------------------------------
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
  #---------------------------------------------------------------------------------------
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
          grDevices::png(paste0(subset_dir, parameters$analysis_name,"_UpSetR_",comparaison,"_allDEG.png"), width=1280, height=1024)
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
          grDevices::png(paste0(subset_dir, parameters$analysis_name,"_UpSetR_",comparaison,"_upDEG.png"), width=1280, height=1024)
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
          grDevices::png(paste0(subset_dir, parameters$analysis_name,"_UpSetR_",comparaison,"_downDEG.png"), width=1280, height=1024)
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
          grDevices::png(paste0(subset_dir, parameters$analysis_name,"_UpSetR_",comparaison,"_mixedDEG.png"), width=1280, height=1024)
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

#' @title VD
#'
#' @description Plot Venn Diagram to compare different contrast
#'
#' @param resDEG, data frame contains for each contrast the significance expression (1/0/-1) for all gene.
#' @param asko_list, list of data.frame contain condition, contrast and context informations made by asko3c.
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @return none.
#'
#' @examples
#' \dontrun{
#'    VD(resDEG, parameters, asko_list)
#' }
#'
#' @export
VD <- function(resDEG, parameters, asko_list){
  options(warn = -1)
  # check parameters
  if(is.null(parameters$VD)==TRUE){ return(NULL) }
  if(is.null(parameters$compaVD)==TRUE || parameters$compaVD==""){
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

  cat("\nCreated VennDiagrams ")
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

#' @title GOenrichment
#'
#' @description Perform GO enrichment analysis with topGO package.
#' This package provides tools for testing GO terms while accounting for
#' the topology of the GO graph. Different test statistics and different
#' methods for eliminating local similarities and dependencies between GO
#' terms can be implemented and applied.
#'
#' @param resDEG, data frame contains for each contrast the significance expression (1/0/-1) for all gene.
#' @param data, list contain all data and metadata (DGEList, samples descriptions, contrast, design and annotations).
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @param list, gene list of interest if you want to apply GOenrichment function on a specific gene list
#' @param title, name of the gene list if you want to apply GOenrichment function on a specific gene list
#' @return none.
#'
#' @import topGO
#' @import goSTAG
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#'    GOenrichment(resDEG, data, parameters, list, title)
#' }
#'
#' @export
GOenrichment<-function(resDEG, data_list, parameters, list=NULL, title=NULL){
  study_dir = paste0(parameters$dir_path, parameters$analysis_name, "/")
  input_path = "/import/"
  GO_dir  = paste0(study_dir, "GOenrichment/")
  if(dir.exists(GO_dir)==FALSE){
    dir.create(GO_dir)
    cat("\n\nDirectory: ",GO_dir," created\n")
  }

  # Get GO annotations
  geneID2GO <- readMappings(file = paste0(input_path,parameters$geneID2GO_file))
  geneNames <- names(geneID2GO)

  if (is.null(list) == FALSE){
    img_go_dir = paste0(GO_dir, "OnSpecificList_",title,"/")
    if(dir.exists(img_go_dir)==FALSE){
      dir.create(img_go_dir)
      cat("Directory: ",img_go_dir," created\n")
    }
    GeneListName = title
    geneSelected = list
  }
  else {
    img_go_dir = paste0(GO_dir, "OnContrasts/")
    if(dir.exists(img_go_dir)==FALSE){
      dir.create(img_go_dir)
      cat("\n\nDirectory: ",img_go_dir," created\n")
    }
    GeneListName = colnames(data_list$contrast)
  }


  for(contrast in GeneListName){


    if (is.null(list) == TRUE){
      if(is.null(parameters$GO)==TRUE){ return(NULL) }

      if(parameters$GO == "both"){
        #print(contrast)
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

    img_GOtoGene_dir = paste0(img_go_dir, contrast,"_SignificantGO_to_Genes/")
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
      #print(GOdata)

      resultTest <- runTest(GOdata, algorithm = parameters$GO_algo, statistic = parameters$GO_stats)
      #print(resultTest)

      resGenTab <- GenTable(GOdata, numChar = 1000000, statisticTest = resultTest, orderBy = "statisticTest", topNodes=length(graph::nodes(graph(GOdata))) )
      resGenTab$Ratio = as.numeric(as.numeric(resGenTab$Significant)/as.numeric(resGenTab$Expected))
      resGenTab$GO_cat <- ontology

      if (is.null(parameters$annotation)==FALSE){
        annot<-read.csv(paste0(input_path, parameters$annotation), header = T, row.names = 1, sep = '\t', quote = "")
      }

      # import normalized MEAN counts in CPM
      moys<-read.csv(paste0(study_dir, parameters$analysis_name,"_CPM_NormMeanCounts.txt"), header=TRUE, sep="\t", row.names=1)
      moys = as.matrix(moys)

      # Create files of genes for each enrichied GO
      myterms = as.character(resGenTab$GO.ID[as.numeric(resGenTab$statisticTest)<=parameters$GO_threshold])

      if (length(myterms) != "0"){
        cat("\nAskoR is saving one file per enriched GO-term (category ", ontology, ").\n")
        mygenes <- genesInTerm(GOdata, myterms)
        noms=names(mygenes)
        nomss=str_replace(noms,":","_")
        for (z in 1:length(mygenes)){
          listes=mygenes[[z]][mygenes[[z]] %in% geneSelected == TRUE]
          GOtab <- data.frame(Gene=listes)
          #GOtab$Gene_cluster = clustered
          rownames(GOtab)=GOtab$Gene
          if (is.null(parameters$annotation)==FALSE){
            GOtab = merge(GOtab, annot, by="row.names")
            GOtab = GOtab[,-1]
            GOtab = GOtab[,1:2]
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
          TabSigCompl<-TabSigCompl[1:maxi,]
        }else{
          cat("\n\n->",contrast," - ontology: ",ontology," - No enrichment can pe performed - there are no feasible GO terms!\n\n")
          TabSigCompl<-as.data.frame(setNames(replicate(8,numeric(0), simplify = F),c("GO.ID","Term","Annotated","Significant","Expected","statisticTest","Ratio","GO_cat") ))
        }

      }else{
        TabCompl=rbind(TabCompl,resGenTab)
        resGenTab[resGenTab=="< 1e-30"]<-"1.0e-30"

        if(nrow(resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold & resGenTab$Significant >= parameters$GO_min_sig_genes,])!=0){
          maxi<-parameters$GO_max_top_terms
          tempSig<-resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold & resGenTab$Significant >= parameters$GO_min_sig_genes,]
          if(maxi > nrow(tempSig)){ maxi<-nrow(tempSig) }
          TabSigCompl=rbind(TabSigCompl,tempSig[1:maxi,])
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
          ggplot(TabSigCompl[TabSigCompl$GO_cat==ontology,], aes(x=stringr::str_wrap(Term, 55), y=Ratio,fill=-1*log10(as.numeric(statisticTest)))) +
            coord_flip()+
            geom_col()+
            theme_classic()+
            geom_text(aes(label=Significant), position=position_stack(0.5),color="white")+
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
          ggsave(filename=paste0(img_go_dir,contrast,"_",ontology,"_GOgraph.png"),width=10, height = 8)
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

        if (is.null(list) == FALSE){
          GraphTitle = paste0("GO Enrichment for list\n",contrast, "\n (",length(which(geneList==1)), " annotated genes among ",length(geneSelected)," genes)")
        }
        else{
          GraphTitle = paste0("GO Enrichment for contrast\n",contrast, "\n (",length(which(geneList==1)), " annotated genes among ",length(geneSelected)," genes)")
        }

        # Ratio Graph
        ggplot(TabSigCompl, aes(x=Ratio, y=Term, size=Significant, color=GO_cat)) +
          geom_point(alpha=1) +
          labs(title = GraphTitle, x="Ratio Significant/Expected", y="GOterm")+
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
        ggplot(TabSigCompl, aes(x=as.numeric(statisticTest), y=Term, size=Significant, color=GO_cat)) +
          geom_point(alpha=1) + labs(title = GraphTitle,x="Pvalue",y="GOterm")+
          scale_color_manual(values=coul,labels=comp_names,name="GO categories")+
          facet_grid(GO_cat~., scales="free", space = "free",labeller = as_labeller(comp_names2))+
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
#' @param asko_norm, large DGEList with normalized counts by GEnorm function.
#' @param resDEG, data frame contains for each contrast the significance expression (1/0/-1) for all genes coming from DEanalysis function.
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @param data, list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations)
#' @param list, gene list of interest if you want to apply ClustAndGO function on a specific gene list
#' @param title, name of the gene list if you want to apply ClustAndGO function on a specific gene list
#' @return clust, data frame with clusters of each gene
#'
#' @example
#' \dontrun{
#'    clust<-ClustAndGO(asko_norm, resDEG, parameters, data, list, title)
#' }
#'
#' @export
ClustAndGO <- function(asko_norm, resDEG, parameters, data, list=NULL, title=NULL){
  study_dir = paste0(parameters$dir_path, parameters$analysis_name, "/")

  CLUST_dir  = paste0(study_dir, "Clustering/")
  if(dir.exists(CLUST_dir)==FALSE){
    dir.create(CLUST_dir)
    cat("\n\nDirectory: ",CLUST_dir," created\n")
  }

  input_path = "/import/"

  if (is.null(list) == TRUE){
    img_Clustering_dir = paste0(CLUST_dir, "OnDEgenes/")
    if(dir.exists(img_Clustering_dir)==FALSE){
      dir.create(img_Clustering_dir)
      cat("Directory: ",img_Clustering_dir," created\n")
    }
  }
  else {
    img_Clustering_dir = paste0(CLUST_dir, "OnSpecificList_",title,"/")
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
  moys<-read.csv(paste0(study_dir, parameters$analysis_name,"_CPM_NormMeanCounts.txt"), header=TRUE, sep="\t", row.names=1)
  moys = as.matrix(moys)

  # import normalized counts (all samples) in CPM
  object=read.csv(paste0(study_dir, parameters$analysis_name,"_CPM_NormCounts.txt"), header=TRUE, sep="\t", row.names=1)



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

  ###############
  ## run coseq ##
  ###############
  library("coseq")

  if (parameters$coseq_ClustersNb > 25){
    detach("package:coseq")
    stop("TOO MANY CLUSTERS : Please set parameters$coseq_ClustersNb to default or under 25")
  }

  if (parameters$coseq_data == 'LogScaledData' & parameters$coseq_transformation != 'none' ){
    detach("package:coseq")
    stop("WRONG TRANSFORMATION CHOSEN  : You try to cluster data already transformed in log2+1. So you don't want to cluster transformed expression profiles (as recommended by coseq creators). Please set parameters$coseq_transformation to 'none' or parameters$coseq_data to 'ExpressionProfiles' ")
  }

  if (length(parameters$coseq_ClustersNb)==1 & parameters$coseq_model=="kmeans"){
    coexpr=coseq(object, K=2:25, model = parameters$coseq_model, transformation = parameters$coseq_transformation,normFactors = "none", seed = 12345)
    clust=as.data.frame(clusters(coexpr, K=parameters$coseq_ClustersNb))
    names(clust)=c("clusters(coexpr)")
  }
  else{
    coexpr=coseq(object, K=parameters$coseq_ClustersNb, model = parameters$coseq_model, transformation = parameters$coseq_transformation,normFactors = "none", seed = 12345)
    clust=as.data.frame(clusters(coexpr))
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

    write.table(GeneToClustersSummary,paste0(img_transfo_dir, parameters$analysis_name,"_ClusteringSUMMARY_",parameters$coseq_model,"_",parameters$coseq_transformation,".txt"),sep="\t",dec=".",row.names = TRUE,col.names = NA)
  }
  else
  {
    write.table(GeneToClustersSummary,paste0(img_transfo_dir, parameters$analysis_name,"_ClusteringSUMMARY_",parameters$coseq_model,"_",parameters$coseq_transformation,".txt"),sep="\t",dec=".",row.names = TRUE,col.names = NA)
  }

  ###################
  ## Global graphs ##
  ###################
  # Boxplots (scaled expression)
  GeneToClustersScaled=GeneToClusters
  GeneToClustersScaled=GeneToClustersScaled[,-2]
  rownames(GeneToClustersScaled)=GeneToClustersScaled$Row.names
  GeneToClustersScaled=GeneToClustersScaled[,-1]
  GeneToClustersScaled=t(apply(as.matrix(GeneToClustersScaled), 1, scale))
  colnames(GeneToClustersScaled)=colnames(GeneToClusters[,-c(1:2)])

  final=data.frame()
  n=as.numeric(ncol(GeneToClustersScaled))
  for (i in 1:n) {
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
    ggsave(filename=paste0(img_transfo_dir, parameters$analysis_name, "_Boxplots_ScaledCPM_",parameters$coseq_model,"_",parameters$coseq_transformation,".png"),width=12,height=8)
  }
  else if (length(unique(final$cluster)) <= 3) {
    ggsave(filename=paste0(img_transfo_dir, parameters$analysis_name, "_Boxplots_ScaledCPM_",parameters$coseq_model,"_",parameters$coseq_transformation,".png"),width=12,height=4)
  }
  else {
    ggsave(filename=paste0(img_transfo_dir, parameters$analysis_name, "_Boxplots_ScaledCPM_",parameters$coseq_model,"_",parameters$coseq_transformation,".png"),width=12,height=12)
  }

  #Global expression profiles with probability (red genes are under proba 0.8)
  proba = 0.80
  for (thresh in proba){
    p <- plot(coexpr, graphs="profiles", K=length(unique(clust$`clusters(coexpr)`)), threshold=thresh)
    p2 <- p$profiles + ggtitle(paste0("Red genes are affiliated to the cluster with a probability lower than ",thresh))

    b <- plot(coexpr, graphs="probapost_barplots", K=length(unique(clust$`clusters(coexpr)`)), threshold=thresh)
    b2 <- b$probapost_barplots + ggtitle(paste0("Threshold probability = ",thresh)) +scale_fill_manual(values = c("black", "red"),name="Max\nconditional\nprobability")

    plot_grid(p2, b2, ncol = 2, nrow = 1,rel_widths = c(1,1))
    ggsave(filename=paste0(img_transfo_dir, parameters$analysis_name, "_ClusterProbabilities_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Threshold",thresh,".png"),width=16,height=8)
  }


  # Heatmap on ScaledCPM
  m=n+2
  mat = as.matrix(GeneToClusters[, 3:m])
  mat_scaled = t(apply(mat, 1, scale))
  colnames(mat_scaled)=colnames(mat)
  rownames(mat_scaled)=GeneToClusters[,1]

  cluster=GeneToClusters[,2]
  min=0
  max=0
  for (i in ncol(mat_scaled)) {
    min2 = min(mat_scaled[,i])
    if (min2<min){
      min =min2
    }
  }
  for (i in ncol(mat_scaled)) {
    max2 = max(mat_scaled[,i])
    if (max2>max){
      max =max2
    }
  }

  if (parameters$coseq_HeatmapOrderSample==TRUE){
    ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")
    ht_list = Heatmap(t(mat_scaled),cluster_rows = FALSE,column_order=order(cluster), name = "Scaled CPM expression",column_split = cluster,
                      heatmap_legend_param = list(title_position = "topcenter",legend_direction = "horizontal"),
                      col = colorRamp2(c(min, 0, max), c("steelblue", "white", "red")),
                      show_column_names = F,
                      column_title_gp = gpar(fill = grey.colors(0.5), col="white", font = 2, fontsize=15),
                      row_gap = unit(2, "mm"), column_gap = unit(2, "mm")
    )
    png(paste0(img_transfo_dir, parameters$analysis_name, "_Heatmap_ScaledCPM_",parameters$coseq_model,"_",parameters$coseq_transformation,"_MySampleOrder.png"), width=sizeImg*1.75, height=sizeImg/4*1.25)
    draw(ht_list,column_title_gp = gpar(font=2, fontsize=20), heatmap_legend_side = "bottom",column_title = paste0("Heatmap on Clusters (parameters : ",parameters$coseq_model," and ",parameters$coseq_transformation," transformation)"))
    dev.off()
  }
  else{
    ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")
    ht_list = Heatmap(t(mat_scaled),column_order=order(cluster), name = "Scaled CPM expression",column_split = cluster,
                      heatmap_legend_param = list(title_position = "topcenter",legend_direction = "horizontal"),
                      col = colorRamp2(c(min, 0, max), c("steelblue", "white", "red")),
                      show_column_names = F,
                      column_title_gp = gpar(fill = grey.colors(0.5), col="white",font=2, fontsize=15),
                      row_gap = unit(2, "mm"), column_gap = unit(2, "mm")
    )
    png(paste0(img_transfo_dir, parameters$analysis_name, "_Heatmap_ScaledCPM_",parameters$coseq_model,"_",parameters$coseq_transformation,".png"), width=sizeImg*1.75, height=sizeImg/4*1.25)
    draw(ht_list,column_title_gp = gpar(font=2, fontsize=20), heatmap_legend_side = "bottom",column_title = paste0("Heatmap on Clusters (parameters : ",parameters$coseq_model," and ",parameters$coseq_transformation," transformation)"))
    dev.off()
  }

  detach("package:coseq")

  #############################
  ## Graphs for each cluster ##
  #############################

  # import data and create vectors for color and cluster
  uniqClust=unique(GeneToClusters$`clusters(coexpr)`)

  GoCoul=c("palegreen", "skyblue", "lightsalmon", "thistle", "tan", "pink", "aquamarine", "violetred", "darkorange", "yellow", "mediumpurple", "wheat","palegreen", "skyblue", "lightsalmon", "thistle", "tan", "pink", "aquamarine", "violetred", "darkorange", "yellow", "mediumpurple", "wheat","palegreen")

  if (is.null(list) == TRUE){
    # create file with matrix of DE genes and cluster for each gene
    resDEG3=resDEG2[which(rowSums(resDEG2)>=1),]
    # delete contrasts with no DE genes
    resDEG3=resDEG3[,(apply(resDEG3,2,sum)!=0)]  }
  else {
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
    ForUpset = ForUpset[,(apply(ForUpset,2,sum)!=0)]
    ForUpset = data.frame(ForUpset)

    if (ncol(ForUpset)>1) {
      png(paste0(img_CLUST_dir,parameters$analysis_name,"_UpSet_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_",clustered,".png"), width=1600, height=1024, units = "px")
      print(upset(data=ForUpset, sets=rev(colnames(ForUpset)), nsets=ncol(ForUpset), keep.order=TRUE, att.color ="black" ,sets.bar.color=GoCoul[clustered],point.size = 5, line.size = 1.5, nintersects=NA, text.scale = 2))
      grid.text(paste0("All differentially expressed genes (up+down) in cluster ",clustered), x=0.65, y=0.95, gp=gpar(fontsize=20))
      dev.off()
    }

    # Scaled expression of each condition in the cluster
    if (nrow(GeneToClusters[GeneToClusters$`clusters(coexpr)`==clustered,])<=750) {alph=1} else {alph=0.2}

    ggplot(final[which(final$cluster==clustered),],aes(x=sample, y=expression))+
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
    ggsave(filename=paste0(img_CLUST_dir,parameters$analysis_name,"_ScaledExpression_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_",clustered,".png"),width=10, height=10)

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
        if (chisq.test(obs1,p=proba)$p.value<0.001 & obs2>=proportion) {
          FileForContrast$ChiTest[a]<-"***"
          FileForContrast$ObservedProportion[a]<-paste0(FileForContrast$ObservedProportion[a],FileForContrast$ChiTest[a])
        }
        else if (chisq.test(obs1,p=proba)$p.value>=0.001 & chisq.test(obs1,p=proba)$p.value<0.01  & obs2>=proportion){
          FileForContrast$ChiTest[a]<-"**"
          FileForContrast$ObservedProportion[a]<-paste0(FileForContrast$ObservedProportion[a],FileForContrast$ChiTest[a])
        }
        else if (chisq.test(obs1,p=proba)$p.value>=0.01 & chisq.test(obs1,p=proba)$p.value<0.05  & obs2>=proportion){
          FileForContrast$ChiTest[a]<-"*"
          FileForContrast$ObservedProportion[a]<-paste0(FileForContrast$ObservedProportion[a],FileForContrast$ChiTest[a])
        }
      }



      ggplot(FileForContrast[FileForContrast$cluster==clustered,], aes(x=contrast, y=GenesOfContrastInCluster)) +
        coord_flip()+
        geom_col(fill=GoCoul[clustered])+
        theme_classic()+
        geom_text(aes(label=ObservedProportion), position=position_stack(0.5),color="black")+
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
      ggsave(filename=paste0(img_CLUST_dir,parameters$analysis_name,"_GenesInContrasts_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_",clustered,".png"),width=10, height = 8)

    }


    # GO enrichment in the cluster for MF, CC, and BP category (if annotation file is provided)
    if(is.null(parameters$geneID2GO_file)==FALSE){

      img_GOtoGene_dir = paste0(img_CLUST_dir,"SignificantGO_to_Genes/")
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
        GOdata <- new("topGOdata",
                      nodeSize = parameters$GO_min_num_genes,
                      ontology = ontology,
                      allGenes = geneList,
                      annot = annFUN.gene2GO,
                      gene2GO = geneID2GO)

        resultTest <- runTest(GOdata, algorithm = parameters$GO_algo, statistic = parameters$GO_stats)
        resGenTab <- GenTable(GOdata, numChar = 1000,statisticTest = resultTest, orderBy = "statisticTest", topNodes=length(nodes(graph(GOdata))) )
        resGenTab$Ratio = as.numeric(as.numeric(resGenTab$Significant)/as.numeric(resGenTab$Expected))
        resGenTab$GO_cat <- ontology

        if (is.null(parameters$annotation)==FALSE){
          annot<-read.csv(paste0(input_path, parameters$annotation), header = T, row.names = 1, sep = '\t', quote = "")
        }

        myterms = as.character(resGenTab$GO.ID[as.numeric(resGenTab$statisticTest)<=parameters$GO_threshold])

        if (length(myterms) != "0"){
          cat("\nAskoR is saving one file per enriched GO-term in cluster ", clustered, " (category ", ontology, ").\n")
          mygenes <- genesInTerm(GOdata, myterms)
          noms=names(mygenes)
          nomss=str_replace(noms,":","_")
          for (z in 1:length(mygenes)){
            listes=mygenes[[z]][mygenes[[z]] %in% GeneToClusters[which(GeneToClusters$`clusters(coexpr)`==clustered),1] == TRUE]
            GOtab <- data.frame(Gene=listes)
            GOtab$Gene_cluster = clustered
            rownames(GOtab)=GOtab$Gene
            if (is.null(parameters$annotation)==FALSE){
              GOtab = merge(GOtab, annot, by="row.names")
              GOtab = GOtab[,-1]
              GOtab = GOtab[,1:3]
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
            write.table(GOtab,paste0(img_GOtoGene_dir, ontology, "_", nomss[z],".txt"), sep="\t", dec=".", row.names = FALSE, col.names = TRUE, quote=FALSE)
          }
        }

        if(ontology == "MF"){
          TabCompl<-resGenTab
          resGenTab[resGenTab=="< 1e-30"]<-"1.0e-30"


          if(nrow(resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold & resGenTab$Significant >= parameters$GO_min_sig_genes,])!=0){
            maxi<-parameters$GO_max_top_terms
            TabSigCompl<-resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold & resGenTab$Significant >= parameters$GO_min_sig_genes,]
            if(maxi > nrow(TabSigCompl)){ maxi<-nrow(TabSigCompl) }
            TabSigCompl<-TabSigCompl[1:maxi,]
          }else{
            cat("\n\n->Cluster ",clustered," - ontology: ",ontology," - No enrichment can pe performed - there are no feasible GO terms!\n\n")
            TabSigCompl<-as.data.frame(setNames(replicate(8,numeric(0), simplify = F),c("GO.ID","Term","Annotated","Significant","Expected","statisticTest","Ratio","GO_cat") ))
          }
        }else{
          TabCompl=rbind(TabCompl,resGenTab)
          resGenTab[resGenTab=="< 1e-30"]<-"1.0e-30"

          if(nrow(resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold & resGenTab$Significant >= parameters$GO_min_sig_genes,])!=0){
            maxi<-parameters$GO_max_top_terms
            tempSig<-resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold & resGenTab$Significant >= parameters$GO_min_sig_genes,]
            if(maxi > nrow(tempSig)){ maxi<-nrow(tempSig) }
            TabSigCompl=rbind(TabSigCompl,tempSig[1:maxi,])
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
            ggplot(TabSigCompl[TabSigCompl$GO_cat==ontology,], aes(x=stringr::str_wrap(Term, 55), y=Ratio,fill=-1*log10(as.numeric(statisticTest)))) +
              coord_flip()+
              geom_col()+
              theme_classic()+
              geom_text(aes(label=Significant), position=position_stack(0.5),color="white")+
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
            ggsave(filename=paste0(img_CLUST_dir,parameters$analysis_name,"_GOEnrichment_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_",clustered,"_", ontology,".png"),width=10, height = 8)
          }
        }
      }

      TabCompl<-TabCompl[TabCompl$Significant > 0,]
      write.table(TabCompl, file=paste0(img_CLUST_dir,parameters$analysis_name,"_GOEnrichmentTable_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_",clustered, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

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
          ggplot(TabSigCompl, aes(x=Ratio, y=Term, size=Significant, color=GO_cat)) +
            geom_point(alpha=1) +
            labs(title = paste0("GO Enrichment for Cluster ",clustered, "\n (",length(which(geneList==1)), " annotated genes among the ",length(which(GeneToClusters[,2]==clustered))," in the cluster)"), x="Ratio Significant / Expected", y="GOterm") +
            scale_color_manual(values=coul,labels=comp_names,name="GO categories") +
            facet_grid(GO_cat~., scales="free", space = "free",labeller = as_labeller(comp_names2)) +
            scale_size_continuous(name="Number of genes") + scale_x_continuous(expand = expansion(add = minR)) +
            scale_y_discrete(labels = function(x) str_wrap(x, 70)) +
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
          ggsave(filename=paste0(img_CLUST_dir,parameters$analysis_name,"_Ratio_BUBBLESgraph_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_",clustered, ".png"),width=10, height=10)
        }
      }else{
        cat("\n\nToo few results to display the graph.\n\n")
      }

    }

  }
  return(clust)
}

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
#' @param data, list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations)
#' @param asko_norm, large DGEList with normalized counts by GEnorm function.
#' @param resDEG, data frame contains for each contrast the significance expression (1/0/-1) for all genes coming from DEanalysis function.
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @param clustering, data frame with clusters of each gene produced by ClustAndGO function
#' @return none
#'
#' @example
#' \dontrun{
#'    IncludeNonDEgenes_InClustering(data, asko_norm, resDEG, parameters, clustering)
#' }
#'
#' @export
IncludeNonDEgenes_InClustering <- function(data, asko_norm, resDEG, parameters, clustering){
  study_dir = paste0(parameters$dir_path, parameters$analysis_name, "/")
  input_path = "/import/"
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
  moys<-read.csv(paste0(study_dir, parameters$analysis_name,"_CPM_NormMeanCounts.txt"), header=TRUE, sep="\t", row.names=1)

  # import GeneToClusters
  GeneToClusters<-read.csv(paste0(img_transfo_dir, parameters$analysis_name,"_ClusteringSUMMARY_",parameters$coseq_model,"_",parameters$coseq_transformation,".txt"), header=TRUE, sep="\t", row.names=1)
  nbCond = length(unique(asko_norm$samples$condition))
  GeneToClusters = GeneToClusters[,1:(1+nbCond)]

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

    write.table(GeneToClustersSummary, paste0(img_transfo_dir, parameters$analysis_name,"_ClusteringSUMMARY_WithNonDEgenes_",parameters$coseq_model,"_",parameters$coseq_transformation,".txt"),sep="\t",dec=".",row.names = TRUE,col.names = NA)
  }
  else
  {
    write.table(GeneToClustersSummary, paste0(img_transfo_dir, parameters$analysis_name,"_ClusteringSUMMARY_WithNonDEgenes_",parameters$coseq_model,"_",parameters$coseq_transformation,".txt"),sep="\t",dec=".",row.names = TRUE,col.names = NA)
  }


  # Boxplots (scaled expression)
  GeneToClustersScaled=GeneToClusters
  GeneToClustersScaled=GeneToClustersScaled[,-1]
  GeneToClustersScaled=t(apply(as.matrix(GeneToClustersScaled), 1, scale))
  colnames(GeneToClustersScaled)=colnames(GeneToClusters[,-1])

  final=data.frame()
  n=as.numeric(ncol(GeneToClustersScaled))
  for (i in 1:n) {
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
  min=0
  max=0
  for (i in ncol(mat_scaled)) {
    min2 = min(mat_scaled[,i])
    if (min2<min){
      min =min2
    }
  }
  for (i in ncol(mat_scaled)) {
    max2 = max(mat_scaled[,i])
    if (max2>max){
      max =max2
    }
  }


  if (parameters$coseq_HeatmapOrderSample==TRUE){
    ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")
    ht_list = Heatmap(t(mat_scaled),cluster_rows = FALSE,column_order=order(cluster), name = "Scaled CPM expression",column_split = cluster,
                      heatmap_legend_param = list(title_position = "topcenter",legend_direction = "horizontal"),
                      col = colorRamp2(c(min, 0, max), c("steelblue", "white", "red")),
                      show_column_names = F,
                      column_title = c(c(1:(length(lab)-1)),"NOT DE"),
                      column_title_gp = gpar(fill = grey.colors(0.5), col="white", font = 2, fontsize=15),
                      row_gap = unit(2, "mm"), column_gap = unit(2, "mm")
    )
    png(paste0(img_transfo_dir, parameters$analysis_name, "_Heatmap_ScaledCPM_WithNonDEgenes_",parameters$coseq_model,"_",parameters$coseq_transformation,"_MySampleOrder.png"), width=sizeImg*1.75, height=sizeImg/4*1.25)
    draw(ht_list,column_title_gp = gpar(font=2, fontsize=20), heatmap_legend_side = "bottom",column_title = paste0("Heatmap on Clusters (parameters : ",parameters$coseq_model," and ",parameters$coseq_transformation," transformation)"))
    dev.off()
  }
  else{
    ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")
    ht_list = Heatmap(t(mat_scaled),column_order=order(cluster), name = "Scaled CPM expression",column_split = cluster,
                      heatmap_legend_param = list(title_position = "topcenter",legend_direction = "horizontal"),
                      col = colorRamp2(c(min, 0, max), c("steelblue", "white", "red")),
                      show_column_names = F,
                      column_title = c(c(1:(length(lab)-1)),"NOT DE"),
                      column_title_gp = gpar(fill = grey.colors(0.5), col="white",font=2, fontsize=15),
                      row_gap = unit(2, "mm"), column_gap = unit(2, "mm")
    )
    png(paste0(img_transfo_dir, parameters$analysis_name, "_Heatmap_ScaledCPM_WithNonDEgenes_",parameters$coseq_model,"_",parameters$coseq_transformation,".png"), width=sizeImg*1.75, height=sizeImg/4*1.25)
    draw(ht_list, column_title_gp = gpar(font=2, fontsize=20),heatmap_legend_side = "bottom",column_title = paste0("Heatmap on Clusters (parameters : ",parameters$coseq_model," and ",parameters$coseq_transformation," transformation)"))
    dev.off()
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

    img_GOtoGene_dir = paste0(img_CLUST_dir,"SignificantGO_to_Genes/")
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
      next
    }

    if(sum(levels(geneList)==1)==0){
      cat("\n -> No DE genes with GO annotation!\n")
      next
    }

    GO=NULL

    listOnto <- c("MF","BP","CC")
    for(ontology in listOnto){
      GOdata <- new("topGOdata",
                    nodeSize = parameters$GO_min_num_genes,
                    ontology = ontology,
                    allGenes = geneList,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)

      resultTest <- runTest(GOdata, algorithm = parameters$GO_algo, statistic = parameters$GO_stats)
      resGenTab <- GenTable(GOdata, numChar = 1000,statisticTest = resultTest, orderBy = "statisticTest", topNodes=length(nodes(graph(GOdata))) )
      resGenTab$Ratio = as.numeric(as.numeric(resGenTab$Significant)/as.numeric(resGenTab$Expected))
      resGenTab$GO_cat <- ontology

      if (is.null(parameters$annotation)==FALSE){
        annot<-read.csv(paste0(input_path, parameters$annotation), header = T, row.names = 1, sep = '\t', quote = "")
      }

      myterms = as.character(resGenTab$GO.ID[as.numeric(resGenTab$statisticTest)<=parameters$GO_threshold])


      if (length(myterms) != "0"){
        cat("\nAskoR is saving one file per enriched GO-term in cluster NOT DE (category ", ontology, ").\n")
        mygenes <- genesInTerm(GOdata, myterms)
        noms=names(mygenes)
        nomss=str_replace(noms,":","_")
        for (z in 1:length(mygenes)){
          listes=mygenes[[z]][mygenes[[z]] %in% rownames(GeneToClusters[which(GeneToClusters$clusters.coexpr.==clustered),]) == TRUE]
          GOtab <- data.frame(Gene=listes)
          GOtab$Gene_cluster = "NOT DE"
          rownames(GOtab)=GOtab$Gene
          if (is.null(parameters$annotation)==FALSE){
            GOtab = merge(GOtab, annot, by="row.names")
            GOtab = GOtab[,-1]
            GOtab = GOtab[,1:3]
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
          write.table(GOtab,paste0(img_GOtoGene_dir, ontology, "_", nomss[z],".txt"), sep="\t", dec=".", row.names = FALSE, col.names = TRUE, quote=FALSE)
        }
      }

      if(ontology == "MF"){
        TabCompl<-resGenTab
        resGenTab[resGenTab=="< 1e-30"]<-"1.0e-30"


        if(nrow(resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold & resGenTab$Significant >= parameters$GO_min_sig_genes,])!=0){
          maxi<-parameters$GO_max_top_terms
          TabSigCompl<-resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold & resGenTab$Significant >= parameters$GO_min_sig_genes,]
          if(maxi > nrow(TabSigCompl)){ maxi<-nrow(TabSigCompl) }
          TabSigCompl<-TabSigCompl[1:maxi,]
        }else{
          cat("\n\n->Cluster NOT DE - ontology: ",ontology," - No enrichment can pe performed - there are no feasible GO terms!\n\n")
          TabSigCompl<-as.data.frame(setNames(replicate(8,numeric(0), simplify = F),c("GO.ID","Term","Annotated","Significant","Expected","statisticTest","Ratio","GO_cat") ))
        }
      }else{
        TabCompl=rbind(TabCompl,resGenTab)
        resGenTab[resGenTab=="< 1e-30"]<-"1.0e-30"

        if(nrow(resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold & resGenTab$Significant >= parameters$GO_min_sig_genes,])!=0){
          maxi<-parameters$GO_max_top_terms
          tempSig<-resGenTab[as.numeric(resGenTab$statisticTest) <= parameters$GO_threshold & resGenTab$Ratio >= parameters$Ratio_threshold & resGenTab$Significant >= parameters$GO_min_sig_genes,]
          if(maxi > nrow(tempSig)){ maxi<-nrow(tempSig) }
          TabSigCompl=rbind(TabSigCompl,tempSig[1:maxi,])
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
          ggplot(TabSigCompl[TabSigCompl$GO_cat==ontology,], aes(x=stringr::str_wrap(Term, 55), y=Ratio,fill=-1*log10(as.numeric(statisticTest)))) +
            coord_flip()+
            geom_col()+
            theme_classic()+
            geom_text(aes(label=Significant), position=position_stack(0.5),color="white")+
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
    write.table(TabCompl, file=paste0(img_CLUST_dir,parameters$analysis_name,"_GOEnrichmentTable_",parameters$coseq_model,"_",parameters$coseq_transformation,"_Cluster_NOT_DE.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

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
        ggplot(TabSigCompl, aes(x=Ratio, y=Term, size=Significant, color=GO_cat)) +
          geom_point(alpha=1) +
          labs(title = paste0("GO Enrichment for Cluster NOT DE \n(",length(which(geneList==1)), " annotated genes among the ",length(which(GeneToClusters[,1]==clustered))," in the cluster)"), x="Ratio Significant / Expected", y="GOterm") +
          scale_color_manual(values=coul,labels=comp_names,name="GO categories") +
          facet_grid(GO_cat~., scales="free", space = "free",labeller = as_labeller(comp_names2)) +
          scale_size_continuous(name="Number of genes") + scale_x_continuous(expand = expansion(add = minR)) +
          scale_y_discrete(labels = function(x) str_wrap(x, 70)) +
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
      }
    }else{
      cat("\n\nToo few results to display the graph.\n\n")
    }


  }
}

#' @title GeneInfo_OnList
#'
#' @description
#' Produce a heatmap of the expression of the genes included in the list with information on DE status of each gene and resume all informations on the genes in a table
#' \itemize{
#'    \item Table with gene expression, DE status, clustering, and gene description
#'    \item Heatmap of gene expression and DE status
#' }
#'
#' @param list, list contain all the genes you want to get information on
#' @param resDEG, data frame contains for each contrast the significance expression (1/0/-1) for all genes coming from DEanalysis function.
#' @param data, list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations)
#' @param title, name of the gene list
#' @param clustering, data frame with clusters of each gene produced by ClustAndGO function
#' @param conditions, list of the conditions you want to see in graph and table
#' @param contrasts, list of the conditions you want to see in graph and table
#' @return none
#'
#' @example
#' \dontrun{
#'    GeneInfo_OnList(list, resDEG, data, title, clustering=NULL, conditions=NULL, contrasts=NULL)
#' }
#'
#' @export
GeneInfo_OnList<-function(list, resDEG, data, title, clustering=NULL, conditions=NULL, contrasts=NULL){
  study_dir  = paste0(parameters$dir_path, parameters$analysis_name, "/")
  input_path = "/import/"

  list_dir = paste0(study_dir, "GeneListExplore/")
  if(dir.exists(list_dir)==FALSE){
    dir.create(list_dir)
    cat("\n\nDirectory: ",list_dir," created\n")
  }

  img_InfosOnGenes_dir = paste0(list_dir, title, "/")
  if(dir.exists(img_InfosOnGenes_dir)==FALSE){
    dir.create(img_InfosOnGenes_dir)
    cat("\n\nDirectory: ",img_InfosOnGenes_dir," created\n")
  }

  # import normalized MEAN counts in CPM
  moys<-read.csv(paste0(study_dir, parameters$analysis_name,"_CPM_NormMeanCounts.txt"), header=TRUE, sep="\t", row.names=1)
  if(is.null(conditions) == FALSE) {
    moys = moys[,colnames(moys) %in% conditions]
  }

  if(is.null(contrasts) == FALSE) {
    resDEG = resDEG[,colnames(resDEG) %in% contrasts]
  }

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
  write.table(totalData,paste0(img_InfosOnGenes_dir, title, "_SummaryGeneList.txt"), sep="\t", dec=".", row.names = TRUE, col.names = NA)

  if(is.null(conditions) == FALSE | is.null(contrasts) == FALSE) {
    suff = "_Subset"
  }
  else{
    suff = "_Complete"
  }


  n = ncol(moys)
  totalDataLOG = log2(totalData[, 1:n]+1)
  mat = as.matrix(totalDataLOG[, 1:n])
  mat_scaled = t(apply(mat, 1, scale))
  colnames(mat_scaled)=colnames(mat)

  min=0
  max=0
  for (i in ncol(mat_scaled)) {
    min2 = min(mat_scaled[,i])
    if (min2<min){
      min =min2
    }
  }
  for (i in ncol(mat_scaled)) {
    max2 = max(mat_scaled[,i])
    if (max2>max){
      max =max2
    }
  }

  if (is.null(clustering)==FALSE){
    graphTitle = "Clusters"
  }
  else{
    graphTitle=""
  }

  hc = rowAnnotation("DE Status in contrasts" = as.matrix(totalData[,(n+1):(n+ncol(resDEG))]),simple_anno_size = unit(0.5, "cm"),gp=gpar(pch=1,col="white",lwd = 4),col = list("DE Status in contrasts" = c("-1" = "green", "0" = "lightgrey", "1" = "red")),
                     annotation_legend_param = list(
                       at = c(-1, 0, 1),
                       legend_height = unit(4, "cm"),
                       title_position = "topleft",
                       legend_side = "bottom", direction ="horizontal")
  )
  ht_opt$TITLE_PADDING = unit(c(7, 7), "points")
  if (is.null(clustering)==FALSE){
    ht_list = Heatmap((mat_scaled), name = "Expression \nLog2(cpm+1)",
                      heatmap_legend_param = list(
                        at = c(-2, 0, 2),
                        legend_height = unit(4, "cm"),
                        title_position = "topleft", direction = "horizontal"
                      ),
                      col = colorRamp2(c(min, 0, max), c("green", "white", "red")),
                      row_split = totalData$CoExpression_Cluster,
                      row_title_gp = gpar(fill = grey.colors(0.5), col="white", font = 2, fontsize=10),
                      row_title_rot = 0,
                      show_row_dend = F,
                      show_column_names = T,
                      show_column_dend = F,
                      row_names_side = "left",
                      column_order = sort(colnames(mat)),
                      row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
                      right_annotation = hc,
                      width=ncol(mat_scaled)*unit(20,"mm"),
                      height=nrow(mat_scaled)*unit(5,"mm")
    )
  }
  else{
    ht_list = Heatmap((mat_scaled), name = "Expression \nLog2(cpm+1)",
                      heatmap_legend_param = list(
                        at = c(-2, 0, 2),
                        legend_height = unit(4, "cm"),
                        title_position = "topleft", direction = "horizontal"
                      ),
                      col = colorRamp2(c(min, 0, max), c("green", "white", "red")),
                      show_row_dend = T,
                      show_column_names = T,
                      show_column_dend = F,
                      row_names_side = "left",
                      column_order = sort(colnames(mat)),
                      row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
                      right_annotation = hc,
                      width=ncol(mat_scaled)*unit(20,"mm"),
                      height=nrow(mat_scaled)*unit(5,"mm")
    )
  }
  if (nrow(mat_scaled)<15){
    png(paste0(img_InfosOnGenes_dir, title, suff, "_heatmap.png"), width=sizeImg*0.8, height=nrow(mat_scaled)*40)
  }
  else if (nrow(mat_scaled)>14 & nrow(mat_scaled)<30){
    png(paste0(img_InfosOnGenes_dir, title, suff, "_heatmap.png"), width=sizeImg*0.8, height=nrow(mat_scaled)*25)
  }
  else if (nrow(mat_scaled)>29 & nrow(mat_scaled)<100){
    png(paste0(img_InfosOnGenes_dir, title, suff, "_heatmap.png"), width=sizeImg*0.8, height=nrow(mat_scaled)*20)
  }
  else{
    png(paste0(img_InfosOnGenes_dir, title, suff, "_heatmap.png"), width=sizeImg*0.8, height=nrow(mat_scaled)*17)
  }
  draw(ht_list, row_title = graphTitle, column_title_gp = gpar(font=2, fontsize=15), heatmap_legend_side = "bottom",column_title = paste0("Expression and DE status \n on genes from list '", title, "'"))
  dev.off()

}


