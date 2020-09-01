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
  pkgs<-c("limma","statmod","edgeR","VennDiagram","RColorBrewer", "UpSetR", "grid", "topGO","ggfortify",
          "Rgraphviz","ggplot2","ggrepel","gplots","stringr","optparse","goSTAG","Glimma","ComplexHeatmap","circlize")
  for(p in pkgs) suppressWarnings(suppressMessages(library(p, quietly=TRUE, character.only=TRUE, warn.conflicts=FALSE)))

  # Specify desired options in a list
  option_list = list(
    make_option(c("-o", "--out"), type="character", default="DE_analysis",dest="analysis_name",
                help="output directory name [default= %default]", metavar="character"),
    make_option(c("-d", "--dir"), type="character", default=".",dest="dir_path",
                help="data directory path [default= %default]", metavar="character"),
    make_option(c("-O", "--org"), type="character", default="Asko", dest="organism",
                help="Organism name [default= %default]", metavar="character"),
    make_option(c("-f", "--fileofcount"), type="character", default=NULL, dest="fileofcount",
                help="file of counts [default= %default]", metavar="character"),
    make_option(c("-G", "--col_genes"), type="integer", default=1, dest="col_genes",
                help="col of genes ids in count files [default= %default]", metavar="integer"),
    make_option(c("-C", "--col_counts"), type="integer", default=7,dest="col_counts",
                help="col of counts in count files [default= %default (featureCounts) ]", metavar="integer"),
    make_option(c("-t", "--sep"), type="character", default="\t", dest="sep",
                help="field separator for count files or count matrix [default= %default]", metavar="character"),
    make_option(c("-a", "--annotation"), type="character", default=NULL, dest="annotation",
                help="annotation file [default= %default]", metavar="character"),
    make_option(c("-s", "--sample"), type="character", default="Samples.txt", dest="sample_file",
                help="Samples file [default= %default]", metavar="character"),
    make_option(c("-c", "--contrasts"), type="character", default="Contrasts.txt",dest="contrast_file",
                help="Contrasts file [default= %default]", metavar="character"),
    make_option(c("-k", "--mk_context"), type="logical", default=FALSE,dest="mk_context",
                help="generate automatically the context names [default= %default]", metavar="logical"),
    make_option(c("--palette"), type="character", default="Set2", dest="palette",
                help="Color palette (ggplot)[default= %default]", metavar="character"),
    make_option(c("-R", "--regex"), type="logical", default=FALSE, dest="regex",
                help="use regex when selecting/removing samples [default= %default]", metavar="logical"),
    make_option(c("-S", "--select"), type="character", default=NULL, dest="select_sample",
                help="selected samples [default= %default]", metavar="character"),
    make_option(c("-r", "--remove"), type="character", default=NULL, dest="rm_sample",
                help="removed samples [default= %default]", metavar="character"),
    make_option(c("--th_cpm"), type="double", default=0.5, dest="threshold_cpm",
                help="CPM's threshold [default= %default]", metavar="double"),
    make_option(c("--rep"), type="integer", default=3, dest="replicate_cpm",
                help="Minimum number of replicates [default= %default]", metavar="integer"),
    make_option(c("--th_FDR"), type="double", default=0.05, dest="threshold_FDR",
                help="FDR threshold [default= %default]", metavar="double"),
    make_option(c("-n", "--normalization"), type="character", default="TMM", dest="normal_method",
                help="normalization method (TMM/RLE/upperquartile/none) [default= %default]", metavar="character"),
    make_option(c("--adj"), type="character", default="fdr", dest="p_adj_method",
                help="p-value adjust method (holm/hochberg/hommel/bonferroni/BH/BY/fdr/none) [default= %default]", metavar="character"),
    make_option("--glm", type="character", default="qlf", dest="glm",
                help="GLM method (lrt/qlf) [default= %default]", metavar="character"),
    make_option("--glmDisp", type="logical", default=FALSE, dest="glm_disp",
                help="Estimate Common, Trended and Tagwise Negative Binomial dispersions GLMs (TRUE/FALSE) [default= %default]", metavar="logical"),
    make_option(c("--lfc"), type="logical", default=TRUE, dest="logFC",
                help="logFC in the summary table [default= %default]", metavar="logical"),
    make_option(c("--th_lfc"), type="double", default=1, dest="threshold_logFC",
                help="logFC threshold [default= %default]", metavar="double"),
    make_option("--fc", type="logical", default=TRUE, dest="FC",
                help="FC in the summary table [default= %default]", metavar="logical"),
    make_option(c("--lcpm"), type="logical", default=FALSE, dest="logCPM",
                help="logCPm in the summary table [default= %default]", metavar="logical"),
    make_option("--fdr", type="logical", default=TRUE, dest="FDR",
                help="FDR in the summary table [default= %default]", metavar="logical"),
    make_option("--lr", type="logical", default=FALSE, dest="LR",
                help="LR in the summary table [default= %default]", metavar="logical"),
    make_option(c("--sign"), type="logical", default=TRUE, dest="Sign",
                help="Significance (1/0/-1) in the summary table [default= %default]", metavar="logical"),
    make_option(c("--expr"), type="logical", default=TRUE, dest="Expression",
                help="Significance expression in the summary table [default= %default]", metavar="logical"),
    make_option(c("--mc"), type="logical", default=TRUE, dest="mean_counts",
                help="Mean counts in the summary table [default= %default]", metavar="logical"),
    make_option(c("--dclust"), type="character", default="euclidean", dest="distcluts",
                help="The distance measure to be used : euclidean, maximum, manhattan, canberra, binary or minkowski [default= %default]", metavar="character"),
    make_option(c("--hclust"), type="character", default="complete", dest="hclust",
                help="The agglomeration method to be used : ward.D, ward.D2, single, complete, average, mcquitty, median or centroid [default= %default]", metavar="character"),
    make_option(c("--hm"), type="logical", default=TRUE, dest="heatmap",
                help="generation of the expression heatmap [default= %default]", metavar="logical"),
    make_option(c("--nh"), type="integer", default="50", dest="numhigh",
                help="number of genes in the heatmap [default= %default]", metavar="integer"),
    make_option(c("--norm_factor"), type="logical", default=FALSE, dest="norm_factor",
                help="generate file with normalize factor for each condition/sample [default= %default]", metavar="logical"),
    make_option(c("--norm_counts"), type="logical", default=FALSE, dest="norm_counts",
                help="Generate files with mormalized counts [default= %default]", metavar="logical"),
    make_option(c("--VD"), type="character", default=NULL, dest="VD",
                help="Plot VennDiagram, precise type of comparison: all, down, up or both [default=%default]", metavar = "character"),
    make_option(c("--compaVD"), type="character", default=NULL, dest="compaVD",
                help="Contrast comparison list to display in VennDiagram", metavar="character"),
    make_option(c("--GO"), type="character", default=NULL, dest="GO",
                help="GO enrichment analysis for gene expressed 'up', 'down' or 'both', NULL for no GO enrichment.", metavar="character"),
    make_option(c("--ID2GO"), type="character", default=NULL, dest="geneID2GO_file",
                help="GO annotation file [default= %default]", metavar="character"),
    make_option(c("--GO_threshold"), type="numeric", default="0.05", dest="GO_threshold",
                help="the significant threshold used to filter p-values", metavar="double"),
    make_option(c("--GO_min_num_genes"), type="integer", default="10", dest="GO_min_num_genes",
                help="the minimum number of genes for each GO terms", metavar="integer"),
    make_option(c("--GO_max_top_terms"), type="integer", default="10", dest="GO_max_top_terms",
                help="the maximum number of GO terms plot", metavar="integer"),
    make_option(c("--GO_algo"), type="character", default="weight01", dest="GO_algo",
                help="algorithms which are accessible via the runTest function: shown by the whichAlgorithms() function, [default=%default]", metavar="character"),
    make_option(c("--GO_stats"), type="character", default="fisher", dest="GO_stats",
                help = "statistical tests which are accessible via the runTest function: shown by the whichTests() function, [default=%default]", metavar = "character"),
    make_option(c("--Ratio_threshold"), type="numeric", default="0", dest="Ratio_threshold",
                help="the minimum ratio value to display GO in graph", metavar="double"),
    make_option(c("--plotMD"),type="logical", default=FALSE, dest="plotMD", metavar="logical",
                help="Mean-Difference Plot of Expression Data (aka MA plot) [default= %default]"),
    make_option(c("--plotVO"),type="logical", default=FALSE, dest="plotVO", metavar="logical",
                help="Volcano plot for a specified coefficient/contrast of a linear model [default= %default]"),
    make_option(c("--glimMD"),type="logical", default=FALSE, dest="glimMD", metavar="logical",
                help="Glimma - Interactif Mean-Difference Plot of Expression Data (aka MA plot) [default= %default]"),
    make_option(c("--glimVO"),type="logical", default=FALSE, dest="glimVO", metavar="logical",
                help="Glimma - Interactif Volcano plot for a specified coefficient/contrast of a linear model [default= %default]"),
    make_option(c("--dens_bottom_mar"), type="integer", default="20", dest="densbotmar", metavar="integer",
                help="Set bottom margin of density plot to help position the legend [default= %default]"),
    make_option(c("--dens_inset"), type="double", default="0.45", dest="densinset", metavar="double",
                help="Set position the legend in bottom density graphe [default= %default]"),
    make_option(c("--legend_col"), type="integer", default="6", dest="legendcol", metavar="integer",
                help="Set numbers of column for density plot legends [default= %default]"),
    make_option(c("--upset_basic"),type="character", default=NULL, dest="upset_basic",
                help="Display UpSetR charts for all contrasts, precise type of comparison: all, down, up, mixed [default=%default].", metavar = "character"),
    make_option(c("--upset_type"),type="character", default=NULL, dest="upset_type",
                help="Display UpSetR charts for list of contrasts, precise type of comparison: all, down, up, mixed [default=%default].", metavar = "character"),
    make_option(c("--upset_list"), type="character", default=NULL, dest="upset_list",
                help="Contrast comparison list to display in UpSetR chart. See documentation. [default=%default]", metavar="character"),
    make_option("--coseq_model", type="character", default="kmeans", dest="coseq_model",
                help="Coseq model (kmeans, Normal) [default= %default]", metavar="character"),
    make_option("--coseq_transformation", type="character", default="clr", dest="coseq_transformation",
                help="Coseq tranformation (arcsin, logit, logclr, clr, alr, ilr, none) [default= %default]", metavar="character"),
    make_option("--coseq_ClustersNb", type="double", default=2:12, dest="coseq_ClustersNb",
                help="Coseq : number of clusters desired (2:12 (auto), number from 2 to 12) [default= %default]", metavar="double"),
    make_option("--coseq_ContrastsThreshold", type="integer", default="1", dest="coseq_ContrastsThreshold",
                help="Coseq : number of contrasts in which DE genes are found for clustering [default= %default]", metavar="integer"),
    make_option("--coseq_normFactors", type="character", default="none", dest="coseq_normFactors",
                help="Coseq normalization factor (TC , UQ, Med , DESeq, TMM , none) [default= %default]", metavar="character")
  )
  # Get command line options
  opt_parser = OptionParser(option_list=option_list)
  parameters = parse_args(opt_parser)

  if(is.null(parameters$rm_sample) == FALSE ) {
    str_replace_all(parameters$rm_sample, " ", "")
    parameters$rm_sample<-strsplit2(parameters$rm_sample, ",")
  }

  if(is.null(parameters$select_sample) == FALSE ) {
    str_replace_all(parameters$select_sample, " ", "")
    parameters$select_sample<-strsplit2(parameters$select_sample, ",")
  }

  return(parameters)
}

#' @title loadData
#'
#' @description
#' Function to load :
#'   - Count data :
#'      ... count matrix : 1 file with all counts for each samples/conditions or mutiple
#'    OR
#'      ... list of files : 1 file of count per conditions, files names contained in sample file
#'   - Metatdata :
#'      ... sample file : file describing the samples and the experimental design
#'      ... contrast file : matrix which specifies which comparisons you would like to make between the samples
#'      ... (optional) annotation file : functional/genomic annotation for each genes
#'      ... (optional) GO terms annotations files : GO annotations for each genes
#'
#' Three output directory will be create :
#' \itemize{
#'   \item images : contains all the images created when running this pipeline with the exception of the upsetR and Venn graphs.
#'   \item vennDiagram : contain all venn diagrams created by VD function
#'   \item UpSetR_graphs : contain all upset graphs created by UpSetGraph function
#'   \item Askomics : files compatible with Askomics Software
#' }
#'
#' @param parameters, list that contains all arguments charged in Asko_start
#' @return data, list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations)
#'
#' @example
#'    data<-loadData(parameters)
#'
#' @export
loadData <- function(parameters){

  # Folders for output files
  #---------------------------------------------------------
  cat("\nCreate directories:\n")
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  if(dir.exists(study_dir)==FALSE){ dir.create(study_dir) }
  cat("\t",study_dir,"\n")

  image_dir = paste0(study_dir, "images/")
  if(dir.exists(image_dir)==FALSE){ dir.create(image_dir) }
  cat("\t",image_dir,"\n")

  if(is.null(parameters$VD)==FALSE){
    venn_dir = paste0(study_dir, "vennDiagram/")
    if(dir.exists(venn_dir)==FALSE){ dir.create(venn_dir) }
    cat("\t",venn_dir,"\n")
  }

  if((is.null(parameters$upset_basic)==FALSE) || (is.null(parameters$upset_list)==FALSE && is.null(parameters$upset_type)==FALSE)){
    upset_dir = paste0(study_dir, "UpSetR_graphs/")
    if(dir.exists(upset_dir)==FALSE){ dir.create(upset_dir) }
    cat("\t",upset_dir,"\n")
  }

  asko_dir = paste0(study_dir, "Askomics/")
  if(dir.exists(asko_dir)==FALSE){ dir.create(asko_dir) }
  cat("\t",asko_dir,"\n")

  # Management of input files
  #---------------------------------------------------------
  input_path = paste0(parameters$dir_path, "/input/")

  # Sample file
  sample_path<-paste0(input_path, parameters$sample_file)
  samples<-read.csv(sample_path, header=TRUE, sep="\t", row.names=1)

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
    else{ color<-colorRampPalette(brewer.pal(11,"Spectral"))(length(condition)) }
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
    dge<-readDGE(paste0(input_path,samples$file), labels=rownames(samples), columns=c(parameters$col_genes,parameters$col_counts), header=TRUE, comment.char="#")
    dge<-DGEList(counts=dge$counts, samples=samples)
  }
  # Matrix file with all counts for all conditions
  else{
    cat("\nSamples:\n")
    print(rownames(samples))
    count_path<-paste0(input_path, parameters$fileofcount)
    if(grepl(".csv", parameters$fileofcount)==TRUE){
      count<-read.csv(count_path, header=TRUE, sep = "\t", row.names = parameters$col_genes, comment.char="#")
    }
    else{
      count<-read.table(count_path, header=TRUE, sep = "\t", row.names = parameters$col_genes, comment.char="#")
    }

    # If you ask for some samples were removed for analysis
    select_counts<-row.names(samples)
    countT<-count[,select_counts]

    # Creates a DGEList object from a table of counts
    dge<-DGEList(counts=countT, samples=samples)
  }

  # Experimental design
  #---------------------------------------------------------
  Group<-factor(samples$condition)
  cat("\nConditions :\n")
  print(Group)
  designExp<-model.matrix(~0+Group)
  rownames(designExp) <- row.names(samples)
  colnames(designExp) <- levels(Group)

  # Contrast for DE analysis
  #---------------------------------------------------------
  contrast_path<-paste0(input_path, parameters$contrast_file)
  contrastab<-read.table(contrast_path, sep="\t", header=TRUE, row.names = 1, comment.char="#", stringsAsFactors = FALSE)

  # Verify if some colunms will be not use for analysis
  rmcol<-list()
  for(condition_name in row.names(contrastab)){
    test<-match(condition_name, colnames(designExp),nomatch = 0)
    if(test==0){
      rm<-grep("0", contrastab[condition_name,], invert = T)
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
    annot<-read.csv(paste0(input_path, parameters$annotation), header = T, row.names = 1, sep = '\t', quote = "")
    data[["annot"]]=annot
  }
  return(data)
}

#' @title asko3c
#'
#' @description Create contrast/condition/context file in format readable by Askomics Software.
#'
#' @param data_list, list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations)
#' @param parameters, list that contains all arguments charged in Asko_start
#' @return asko, list of data.frame contain condition, contrast and context informations
#'
#' @example
#'    asko<-asko3c(data_list, parameters)
#'
#' @export
asko3c <- function(data_list, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  asko_dir = paste0(study_dir, "Askomics/")
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

      if(length(set_cond1)==1){complex1=F}else{complex1=T}                        # to determine if we have complex contrast (multiple conditions
      if(length(set_cond2)==1){complex2=F}else{complex2=T}                        # compared to multiple conditions) or not
      if(complex1==F && complex2==F){                                             # Case 1: one condition against one condition
        contrast_asko[i,"context1"]<-set_cond1                                    # filling contrast data frame with the name of the 1st context
        contrast_asko[i,"context2"]<-set_cond2                                    # filling contrast data frame with the name of the 2nd context
        contrast_name<-paste(set_cond1,set_cond2, sep = "vs")                     # creation of contrast name by associating the names of contexts
        contrast_asko[i,"Contrast"]<-contrast_name                                # filling contrast data frame with contrast name
        list_context<-append(list_context, set_cond1)                             #
        list_condition<-append(list_condition, set_cond1)                         # adding respectively to the lists "context" and "condition" the context name
        list_context<-append(list_context, set_cond2)                             # and the condition name associated
        list_condition<-append(list_condition, set_cond2)                         #
      }
      if(complex1==F && complex2==T){                                             # Case 2: one condition against multiple condition
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
            verif<-unique(str_detect(set_cond2, value))                           # verification of the presence of values in each condition contained in the set
            if(length(verif)==1 && verif==TRUE){common_factor[l]<-value}          # if verif contains only TRUE, value of experimental factor
          }                                                                       # is added as common factor
        }
        if(length(common_factor)>1){                                              # if there are several common factor
          common_factor<-toString(common_factor)                                  # the list is converted to string
          contx<-str_replace(common_factor,", ","")
          contx<-str_replace_all(contx, "NULL", "")}else{contx<-common_factor}    # and all common factor are concatenated to become the name of context
        contrast_asko[i,"context2"]<-contx                                        # filling contrast data frame with the name of the 2nd context
        contrast_name<-paste(set_cond1,contx, sep = "vs")                         # concatenation of context names to make the contrast name
        contrast_asko[i,"Contrast"]<-contrast_name                                # filling contrast data frame with the contrast name
        for(j in seq_along(set_cond2)){                                            # for each condition contained in the complex context (2nd):
          list_context<-append(list_context, contx)                               # adding condition name with the context name associated
          list_condition<-append(list_condition, set_cond2[j])                    # to their respective list
        }
      }
      if(complex1==T && complex2==F){                                             # Case 3: multiple conditions against one condition
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
            verif<-unique(str_detect(set_cond1, value))                           # verification of the presence of values in each condition contained in the set
            if(length(verif)==1 && verif==TRUE){common_factor[l]<-value}          # if verif contains only TRUE, value of experimental factor
          }                                                                       # is added as common factor
        }
        if(length(common_factor)>1){                                              # if there are several common factor
          common_factor<-toString(common_factor)                                  # the list is converted to string
          contx<-str_replace(common_factor,", ","")
          contx<-str_replace_all(contx, "NULL", "")}else{contx<-common_factor}    # and all common factor are concatenated to become the name of context
        contrast_asko[i,"context1"]<-contx                                        # filling contrast data frame with the name of the 1st context
        contrast_name<-paste(contx,set_cond2, sep = "vs")                         # concatenation of context names to make the contrast name
        contrast_asko[i,"Contrast"]<-contrast_name                                # filling contrast data frame with the contrast name
        for(j in seq_along(set_cond1)){                                            # for each condition contained in the complex context (1st):
          list_context<-append(list_context, contx)                               # adding condition name with the context name associated
          list_condition<-append(list_condition, set_cond1[j])                    # to their respective list
        }
      }
      if(complex1==T && complex2==T){                                             # Case 4: multiple conditions against multiple conditions
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
            verif1<-unique(str_detect(set_cond1, value))                          # verification of the presence of values in each condition contained in the 1st context
            verif2<-unique(str_detect(set_cond2, value))                          # verification of the presence of values in each condition contained in the 2nd context

            if(length(verif1)==1 && verif1==TRUE){m=m+1;common_factor1[m]<-value} # if verif=only TRUE, value of experimental factor is added as common factor
            if(length(verif2)==1 && verif2==TRUE){n=n+1;common_factor2[n]<-value} # if verif=only TRUE, value of experimental factor is added as common factor
          }
        }
        if(length(common_factor1)>1){                                             # if there are several common factor for conditions in the 1st context
          common_factor1<-toString(common_factor1)                                # conversion list to string
          contx1<-str_replace(common_factor1,", ","")}else{contx1<-common_factor1}# all common factor are concatenated to become the name of context
        contx1<-str_replace_all(contx1, "NULL", "")
        if(length(common_factor2)>1){                                             # if there are several common factor for conditions in the 2nd context
          common_factor2<-toString(common_factor2)                                # conversion list to string
          contx2<-str_replace(common_factor2,", ","")}else{contx2<-common_factor2}# all common factor are concatenated to become the name of context
        contx2<-str_replace_all(contx2, "NULL", "")
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
      contexts=strsplit2(contrast,"vs")
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
  ctime<-format(Sys.time(), "%d-%m-%Y_%Hh%Mm%Ss")
  # creation of condition file for asko
  write.table(data.frame("Condition"=row.names(condition_asko),condition_asko),
              paste0(asko_dir,"condition.asko",ctime,".txt"),
              sep = parameters$sep,
              row.names = F,
              quote=F)
  # creation of context file for asko
  write.table(context_asko,
              paste0(asko_dir, "context.asko",ctime,".txt"),
              sep=parameters$sep,
              col.names = T,
              row.names = F,
              quote=F)
  # creation of contrast file for asko
  write.table(contrast_asko,
              paste0(asko_dir, "contrast.asko",ctime,".txt"),
              sep=parameters$sep,
              col.names = T,
              row.names = F,
              quote=F)
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
#' @example
#'    filtered_counts<-GEfilt(data_list, parameters)
#'
#' @export
GEfilt <- function(data_list, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "images/")

  # plot density before filtering
  #---------------------------------
  cpm<-cpm(data_list$dge)
  logcpm<-cpm(data_list$dge, log=TRUE)
  colnames(logcpm)<-rownames(data_list$dge$samples)
  nsamples <- ncol(data_list$dge$counts)

  maxi<-c()
  for (i in seq(nsamples)){
    m=max(density(logcpm[,i])$y)
    maxi<-c(maxi,m)
  }
  ymax<-round(max(maxi),1) + 0.02

  sizeImg=15*nsamples
  if(sizeImg < 480){ sizeImg=480 }
  btm=round((nsamples/6),0)+0.5

  png(paste0(image_dir, parameters$analysis_name, "_raw_data.png"), width=1024, height=1024)
  par(oma=c(2,2,2,0), mar=c(parameters$densbotmar,5,5,5))
  plot(density(logcpm[,1]),
       col=as.character(data_list$dge$samples$color[1]),
       lwd=1,
       las=2,
       ylim=c(0,ymax),
       main="A. Raw data",
       xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den<-density(logcpm[,i])
    lines(den$x, col=as.character(data_list$dge$samples$color[i]), den$y, lwd=1)
  }

  legend("bottom", fill=data_list$dge$samples$color, bty="n", ncol=parameters$legendcol,
         legend=rownames(data_list$dge$samples), xpd=TRUE, inset=-parameters$densinset)
  dev.off()

  # plot density after filtering
  #---------------------------------
  keep.exprs <- rowSums(cpm>parameters$threshold_cpm)>=parameters$replicate_cpm
  filtered_counts <- data_list$dge[keep.exprs,,keep.lib.sizes=F]
  filtered_cpm<-cpm(filtered_counts$counts, log=TRUE)

  maxi<-c()
  for (i in seq(nsamples)){
    m=max(density(filtered_cpm[,i])$y)
    maxi<-c(maxi,m)
  }
  ymax<-round(max(maxi),1) + 0.02

  png(paste0(image_dir,parameters$analysis_name,"_filtered_data.png"), width=1024, height=1024)
  par(oma=c(2,2,2,0), mar=c(parameters$densbotmar,5,5,5))
  plot(density(filtered_cpm[,1]),
       col=as.character(data_list$dge$samples$color[1]),
       lwd=1,
       ylim=c(0,ymax),
       las=2,
       main="B. Filtered data",
       xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- density(filtered_cpm[,i])
    lines(den$x,col=as.character(data_list$dge$samples$color[i]), den$y, lwd=1)
  }
  legend("bottom", fill=data_list$dge$samples$color, bty="n", ncol=parameters$legendcol,
         legend=rownames(data_list$dge$samples), xpd=TRUE, inset=-parameters$densinset)
  dev.off()

  # histogram cpm values distribution before filtering
  #------------------------------------------------------
  png(paste0(image_dir,parameters$analysis_name,"_barplot_logcpm_before_filtering.png"), width=sizeImg, height=sizeImg)
  hist(logcpm,
       main= "A. Log2(cpm) distribution before filtering",
       xlab = "log2(cpm)",
       col = "grey")
  dev.off()

  # histogram cpm values distribution after filtering
  #------------------------------------------------------
  png(paste0(image_dir,parameters$analysis_name,"_barplot_logcpm_after_filtering.png"), width=sizeImg, height=sizeImg)
  hist(filtered_cpm,
       main= "B. Log2(cpm) distribution after filtering",
       xlab = "log2(cpm)",
       col = "grey")
  dev.off()

  # boxplot cpm values distribution before filtering
  #------------------------------------------------------
  png(paste0(image_dir,parameters$analysis_name,"_boxplot_logcpm_before_filtering.png"), width=sizeImg, height=sizeImg)
  par(oma=c(1,1,1,1))
  boxplot(logcpm,
          col=data_list$dge$samples$color,
          main="A. Log2(cpm) distribution before filtering",
          cex.axis=0.8,
          las=2,
          ylab="log2(cpm)")
  dev.off()

  # boxplot cpm values distribution after filtering
  #------------------------------------------------------
  png(paste0(image_dir,parameters$analysis_name,"_boxplot_logcpm_after_filtering.png"), width=sizeImg, height=sizeImg)
  par(oma=c(1,1,1,1))
  boxplot(filtered_cpm,
          col=data_list$dge$samples$color,
          main="B. Log2(cpm) distribution after filtering",
          cex.axis=0.8,
          las=2,
          ylab="log2(cpm)")
  dev.off()

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
#' @example
#'    norm_GE<-GEnorm(filtered_GE, asko_list, parameters)
#'
#' @export
GEnorm <- function(filtered_GE, asko_list, data, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "images/")

  # for image size
  nsamples <- ncol(filtered_GE$counts)
  sizeImg=15*nsamples
  if(sizeImg < 480){ sizeImg=480 }

  # Normalization counts
  norm_GE<-calcNormFactors(filtered_GE, method = parameters$normal_method)
  if(parameters$norm_factor == TRUE){
    write.table(norm_GE$samples, file=paste0(study_dir, parameters$analysis_name, "_normalization_factors.txt"), col.names=NA, row.names=TRUE, quote=FALSE, sep="\t", dec=".")
  }

  # boxplot log2(cpm) values after normalization
  #----------------------------------------------------
  logcpm_norm <- cpm(norm_GE, log=TRUE)
  colnames(logcpm_norm)<-rownames(filtered_GE$samples)

  png(paste0(image_dir,parameters$analysis_name,"_boxplot_logcpm_after_norm.png"), width=sizeImg, height=sizeImg)
  par(oma=c(1,1,1,1))
  boxplot(logcpm_norm,
          col=filtered_GE$samples$color,
          main="B. Log2(cpm) distribution after normalization",
          cex.axis=0.8,
          las=2,
          ylab="Log2(cpm)")
  dev.off()

  # heatmap visualisation
  #----------------------------------------------------
  if(nrow(filtered_GE$counts) <= 30000)
  {
    # heatmap cpm value per sample
    #----------------------------------------------------
    cpm_norm  <- cpm(norm_GE, log=FALSE)
    cpmscale  <- scale(t(cpm_norm))
    tcpmscale <- t(cpmscale)

    d1 <- dist(cpmscale,  method = parameters$distcluts, diag = FALSE, upper = FALSE)
    d2 <- dist(tcpmscale, method = parameters$distcluts, diag = FALSE, upper = TRUE)
    hc <- hclust(d1, method = parameters$hclust, members = NULL)
    hr <- hclust(d2, method = parameters$hclust, members = NULL)
    my_palette <- colorRampPalette(c("green","black","red"), interpolate = "linear")

    png(paste0(image_dir,parameters$analysis_name,"_heatmap_CPMcounts_per_sample.png"), width=sizeImg*1.5, height=sizeImg*1.25)
    par(oma=c(2,1,2,2))
    heatmap.2(tcpmscale, Colv = as.dendrogram(hc), Rowv = as.dendrogram(hr), density.info="histogram",
              trace = "none", dendrogram = "column", xlab = "samples", col = my_palette, labRow = FALSE,
              cexRow = 0.1, cexCol = 1.25, ColSideColors = norm_GE$samples$color, margins = c(10,1),
              main = paste0("CPM counts per sample\nGenes 1 to ",nrow(norm_GE)))
    dev.off()
    write.table(cpm_norm, file=paste0(study_dir, parameters$analysis_name, "_CPM_NormCounts.txt"), col.names=NA, row.names=TRUE, quote=FALSE, sep="\t", dec=".")

    # Normalized mean by conditions
    #-------------------------------
    # heatmap mean counts per condition
    n_count <- NormCountsMean(norm_GE, ASKOlist = asko_list)
    countscale  <- scale(t(n_count))
    tcountscale <- t(countscale)

    d1 <- dist(countscale,  method = parameters$distcluts, diag = FALSE, upper = FALSE)
    d2 <- dist(tcountscale, method = parameters$distcluts, diag = FALSE, upper = TRUE)
    hc <- hclust(d1, method = parameters$hclust, members = NULL)
    hr <- hclust(d2, method = parameters$hclust, members = NULL)
    my_palette <- colorRampPalette(c("green","black","red"), interpolate = "linear")

    png(paste0(image_dir,parameters$analysis_name,"_heatmap_meanCounts_per_condi.png"), width=sizeImg*1.5, height=sizeImg*1.25)
    par(oma=c(2,1,2,2))
    heatmap.2(tcountscale, Colv = as.dendrogram(hc), Rowv = as.dendrogram(hr), density.info="histogram",
              trace = "none", dendrogram = "column", xlab = "Condition", col = my_palette, labRow = FALSE,
              cexRow = 0.1, cexCol = 1.5, ColSideColors = unique(norm_GE$samples$color), margins = c(10,1),
              main = paste0("Mean count per condition\nGenes 1 to ",nrow(norm_GE)))
    dev.off()
  }

  # File with normalized counts
  if (parameters$norm_counts == TRUE){
    dScaleFactors <- norm_GE$samples$lib.size * norm_GE$samples$norm.factors
    normCounts <- t(t(filtered_GE$counts)/dScaleFactors)*mean(dScaleFactors)
    write.table(normCounts, file=paste0(study_dir, parameters$analysis_name, "_NormCounts.txt"), col.names=NA, row.names=TRUE, quote=FALSE, sep="\t", dec=".")
  }

  # File with normalized mean counts in CPM, grouped by condition
  meancpmDEGnorm<-as.data.frame(t(aggregate(t(cpm_norm),list(data$samples$condition), mean)))[-1,]
  colnames(meancpmDEGnorm)<-unique(data$samples$condition)
  write.table(meancpmDEGnorm,paste0(study_dir, parameters$analysis_name,"_CPM_NormMeanCounts.txt"), sep="\t", dec=".", row.names=TRUE, col.names=NA, quote=FALSE)

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
#' @example
#'    GEcorr(asko_norm,parameters)
#'
#' @export
GEcorr <- function(asko_norm, parameters){
  options(warn = -1)
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "images/")

  # for image size
  nsamples <- ncol(asko_norm$counts)
  sizeImg=15*nsamples
  if(sizeImg < 1024){ sizeImg=1024 }

  lcpm<-cpm(asko_norm, log=TRUE)
  colnames(lcpm)<-rownames(asko_norm$samples)

  # Heatmap sample correlation
  #-----------------------------
  cormat<-cor(lcpm)
  color<-colorRampPalette(c("black","red","yellow","white"),space="rgb")(35)
  png(paste0(image_dir, parameters$analysis_name, "_heatmap_of_sample_correlation.png"), width=sizeImg, height=sizeImg)
  par(oma=c(4,2,4,1))
  heatmap(cormat, col=color, symm=TRUE, RowSideColors=as.character(asko_norm$samples$color),
          ColSideColors=as.character(asko_norm$samples$color), main="")
  title("Sample Correlation Matrix", adj=0.5, outer=TRUE)
  dev.off()

  # MDS Plot
  #-----------------------------
  mds <- cmdscale(dist(t(lcpm)),k=3, eig=TRUE)
  eigs<-round((mds$eig)*100/sum(mds$eig[mds$eig>0]),2)
  dfmds<-as.data.frame(mds$points)
  # Axe 1 and 2
  png(paste0(image_dir, parameters$analysis_name, "_MDS_corr_axe1_2.png"), width=sizeImg*1.25, height=sizeImg*1.25)
  mds1<-ggplot(dfmds, aes(dfmds$V1, dfmds$V2, label=rownames(mds$points))) + labs(title="MDS Axes 1 and 2") +
    theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin=margin(20,30,20,15)) +
    geom_point(aes(color=as.character(asko_norm$samples$condition)),shape=17,size=15 ) +
    xlab(paste('dim 1 [', eigs[1], '%]')) + ylab(paste('dim 2 [', eigs[2], "%]")) +
    geom_label_repel(aes(label = rownames(mds$points),fill = factor(as.character(asko_norm$samples$condition))), color = 'white',size = 7) +
    theme(
      axis.text.y = element_text(face="bold",size=30),
      axis.text.x = element_text(face="bold",size=30),
      legend.text = element_text(size = 30),
      legend.title = element_text(size=30,face="bold"),
      strip.text.y = element_text(size=30, face="bold"),
      axis.title = element_text(size=30,face="bold"),
      plot.title = element_text(size=35)) +
    guides(color=FALSE, fill = guide_legend(title = "Condition", override.aes = aes(label = ""), ncol=parameters$legendcol)) +
    theme(legend.position="bottom",legend.direction="vertical",legend.margin=margin(5,5,5,5),legend.box.margin=margin(10,10,10,10))
  print(mds1)
  dev.off()
  # Axe 2 and 3
  png(paste0(image_dir, parameters$analysis_name, "_MDS_corr_axe2_3.png"), width=sizeImg*1.25, height=sizeImg*1.25)
  mds2<-ggplot(dfmds, aes(dfmds$V2, dfmds$V3, label = rownames(mds$points))) + labs(title="MDS Axes 2 and 3") +
    theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin=margin(20,30,20,15)) +
    geom_point(aes(color=as.character(asko_norm$samples$condition)),shape=17,size=15 ) +
    xlab(paste('dim 2 [', eigs[2], '%]')) + ylab(paste('dim 3 [', eigs[3], "%]")) +
    geom_label_repel(aes(label = rownames(mds$points),fill = factor(as.character(asko_norm$samples$condition))), color = 'white',size = 7) +
    theme(
      axis.text.y = element_text(face="bold",size=30),
      axis.text.x = element_text(face="bold",size=30),
      legend.text = element_text(size = 30),
      legend.title = element_text(size=30,face="bold"),
      strip.text.y = element_text(size=30, face="bold"),
      axis.title=element_text(size=30,face="bold"),
      plot.title = element_text(size=35)) +
    guides(color=FALSE, fill = guide_legend(title = "Condition", override.aes = aes(label = ""), ncol=parameters$legendcol)) +
    theme(legend.position="bottom",legend.direction="vertical",legend.margin=margin(5,5,5,5),legend.box.margin=margin(10,10,10,10))
  print(mds2)
  dev.off()
  # Axe 1 and 3
  png(paste0(image_dir, parameters$analysis_name, "_MDS_corr_axe1_3.png"), width=sizeImg*1.25, height=sizeImg*1.25)
  mds3<-ggplot(dfmds, aes(dfmds$V1, dfmds$V3, label = rownames(mds$points))) + labs(title="MDS Axes 1 and 3") +
    theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin=margin(20,30,20,15)) +
    geom_point(aes(color=as.character(asko_norm$samples$condition)),shape=17,size=15 ) +
    xlab(paste('dim 1 [', eigs[1], '%]')) + ylab(paste('dim 3 [', eigs[3], "%]")) +
    geom_label_repel(aes(label = rownames(mds$points),fill = factor(as.character(asko_norm$samples$condition))), color = 'white',size = 7) +
    theme(
      axis.text.y = element_text(face="bold",size=30),
      axis.text.x = element_text(face="bold",size=30),
      legend.text = element_text(size = 30),
      legend.title = element_text(size=30,face="bold"),
      strip.text.y = element_text(size=30, face="bold"),
      axis.title=element_text(size=30,face="bold"),
      plot.title = element_text(size=35)) +
    guides(color=FALSE, fill = guide_legend(title = "Condition",override.aes = aes(label = ""), ncol=parameters$legendcol)) +
    theme(legend.position="bottom",legend.direction="vertical",legend.margin=margin(5,5,5,5),legend.box.margin=margin(10,10,10,10))
  print(mds3)
  dev.off()

  # PCA
  #-----------------------------
  lcpm2=as.data.frame(t(lcpm))
  lcpm3=lcpm2
  lcpm3$cond=asko_norm$samples$condition

  # Axe 1 and 2
  png(paste0(image_dir, parameters$analysis_name, "_PCA_axe1_2.png"), width=sizeImg*1.25, height=sizeImg*1.25)
  p=autoplot(prcomp(lcpm2), data = lcpm3, size=15,colour="cond",show.legend = FALSE, x=1, y=2,shape=17)
  test=p + geom_label_repel(aes(label = rownames(lcpm3),fill = factor(cond)), color = 'white',size = 7) +
  labs(title="PCA Axes 1 and 2") + theme(
    plot.margin=margin(20,30,20,15),
    axis.text.y = element_text(face="bold",size=30),
    axis.text.x = element_text(face="bold",size=30),
    legend.text = element_text(size = 30),
    legend.title = element_text(size=30,face="bold"),
    strip.text.y = element_text(size=30, face="bold"),
    axis.title=element_text(size=30,face="bold"),
    plot.title = element_text(size=35, hjust=0.5)) +
    guides(colour=FALSE, fill = guide_legend(title = "Condition",override.aes = aes(label = ""), ncol=parameters$legendcol)) +
    theme(legend.position="bottom",legend.direction="vertical",legend.margin=margin(5,5,5,5),legend.box.margin=margin(10,10,10,10))
  print(test)
  dev.off()

  # Axe 2 and 3
  png(paste0(image_dir, parameters$analysis_name, "_PCA_axe2_3.png"), width=sizeImg*1.25, height=sizeImg*1.25)
  p=autoplot(prcomp(lcpm2), data = lcpm3, size=15,colour="cond",show.legend = FALSE, x=2, y=3,shape=17)
  test=p + geom_label_repel(aes(label = rownames(lcpm3),fill = factor(cond)), color = 'white',size = 7) +
  labs(title="PCA Axes 2 and 3") + theme(
    plot.margin=margin(20,30,20,15),
    axis.text.y = element_text(face="bold",size=30),
    axis.text.x = element_text(face="bold",size=30),
    legend.text = element_text(size = 30),
    legend.title = element_text(size=30,face="bold"),
    strip.text.y = element_text(size=30, face="bold"),
    axis.title=element_text(size=30,face="bold"),
    plot.title = element_text(size=35, hjust=0.5)) +
    guides(colour=FALSE, fill = guide_legend(title = "Condition",override.aes = aes(label = ""), ncol=parameters$legendcol)) +
    theme(legend.position="bottom",legend.direction="vertical",legend.margin=margin(5,5,5,5),legend.box.margin=margin(10,10,10,10))
  print(test)
  dev.off()

  # Axe 1 and 3
  png(paste0(image_dir, parameters$analysis_name, "_PCA_axe1_3.png"), width=sizeImg*1.25, height=sizeImg*1.25)
  p=autoplot(prcomp(lcpm2), data = lcpm3, size=15,colour="cond",show.legend = FALSE, x=1, y=3,shape=17)
  test=p + geom_label_repel(aes(label = rownames(lcpm3),fill = factor(cond)), color = 'white',size = 7) +
  labs(title="PCA Axes 1 and 3") + theme(
    plot.margin=margin(20,30,20,15),
    axis.text.y = element_text(face="bold",size=30),
    axis.text.x = element_text(face="bold",size=30),
    legend.text = element_text(size = 30),
    legend.title = element_text(size=30,face="bold"),
    strip.text.y = element_text(size=30, face="bold"),
    axis.title=element_text(size=30,face="bold"),
    plot.title = element_text(size=35, hjust=0.5)) +
    guides(colour=FALSE, fill = guide_legend(title = "Condition",override.aes = aes(label = ""), ncol=parameters$legendcol)) +
    theme(legend.position="bottom",legend.direction="vertical",legend.margin=margin(5,5,5,5),legend.box.margin=margin(10,10,10,10))
  print(test)
  dev.off()

  # hierarchical clustering
  #-----------------------------
  mat.dist <- dist(t(asko_norm$counts), method = parameters$distcluts)
  clustering <- hclust(mat.dist, method=parameters$hclust)
  png(paste0(image_dir, parameters$analysis_name, "_hclust.png"), width=sizeImg, height=sizeImg)
  par(oma=c(1,1,1,1))
  plot(clustering,
       main = 'Distances Correlation\nHierarchical clustering', sub="",
       xlab=paste0("Samples\nmethod hclust: ",parameters$hclust),
       hang = -1)
  j<-dev.off()
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
#'    plot_glimma(fit, normGE, resDE, contrast, "MD", parameters)  # smear plot
#'    plot_glimma(fit, normGE, resDE, contrast, "VO", parameters)  # volcano plot
#'
plot_glimma <- function(fit, normGE, resDE, contrast, tplot, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "images/")

  # Mean-Difference Plot
  #--------------------------
  if(tplot=="MD"){
    if (is.null(normGE$samples$color)==TRUE){
      suppressWarnings(glMDPlot(fit, status=resDE[,contrast], counts=normGE, group=normGE$samples$condition,
                                transform=TRUE, anno=NULL, launch=FALSE, main=contrast,
                                folder=paste0(image_dir, "Glimma_Plots"), html=paste0("MDPlot_",contrast)))
    }
    else{
      suppressWarnings(glMDPlot(fit, status=resDE[,contrast], counts=normGE, group=normGE$samples$condition, transform=TRUE,
                                sample.cols=normGE$samples$color, anno=NULL, launch=FALSE, main=contrast,
                                folder=paste0(image_dir, "Glimma_Plots"), html=paste0("MDPlot_",contrast)))
    }
  }

  # Volcano plot
  #--------------------------
  if (tplot=="VO"){
    if (is.null(normGE$samples$color)==TRUE){
      glXYPlot(x=fit$table$logFC, y=-log10(fit$table$PValue), status=resDE[,contrast], counts=normGE,
               group=normGE$samples$condition, xlab="Log2FoldChange", ylab="-log10(pvalue)", main=contrast,
               launch=FALSE, folder=paste0(image_dir, "Glimma_Plots"), html=paste0("Volcano_",contrast))
    }
    else{
      glXYPlot(x=fit$table$logFC, y=-log10(fit$table$PValue), status=resDE[,contrast], counts=normGE, main=contrast,
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
#'    plot_expr(fit, normGE, resDE, contrast, "MD", parameters)  # smear plot
#'    plot_expr(fit, normGE, resDE, contrast, "VO", parameters)  # volcano plot
#'
plot_expr <- function(fit, normGE, resDE, contrast, tplot, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "images/")



  # Mean-Difference Plot
  if(tplot=="MD"){
    png(paste0(image_dir, contrast, "_MeanDifference_of_ExpressionData.png"))
    plotMD.DGELRT(fit, xlab="Average log CPM", ylab="log-fold-change", main=paste0("MD plot - ", contrast),
                  cex=0.5, status=resDE[,contrast], values=c("-1","1"), col=c("blue","red"))
    dev.off()
  }

  # Volcano plot
  if(tplot=="VO"){
    tglm<-fit$table
    tglm$FDR<-p.adjust(tglm$PValue, method=parameters$p_adj_method)
    png(paste0(image_dir, contrast, "_VolcanoPlot.png"))
    with(tglm, plot(tglm$logFC, -log10(tglm$PValue), pch=16, cex=0.5, xlim=c(min(tglm$logFC)-0.5, max(tglm$logFC)+0.5),
                    ylim=c(min(-log10(tglm$PValue))-0.5, max(-log10(tglm$PValue))+0.5),
                    main=paste0("Volcano plot - ", contrast), xlab="Log2FoldChange", ylab="-log10(pvalue)"))
    with(subset(tglm, tglm$FDR <= parameters$threshold_FDR & tglm$logFC >  parameters$threshold_logFC), points(logFC, -log10(PValue), pch=16, cex=0.5, col="red"))
    with(subset(tglm, tglm$FDR <= parameters$threshold_FDR & tglm$logFC < -parameters$threshold_logFC), points(logFC, -log10(PValue), pch=16, cex=0.5, col="blue"))
    # abline(h=-log10(parameters$threshold_FDR), v=c(-parameters$threshold_logFC,parameters$threshold_logFC), col="darkgreen")
    dev.off()
  }
}

#' @title NormCountsMean
#'
#' @description Calcul mean counts for two contrast or all matrix.
#'
#'  @param glmfit, fitted linear model object.
#'  @param ASKOlist, list of data.frame contain condition, contrast and context informations made by asko3c.
#'  @param context, coefficient/contrast tested.
#'  @return
#'      matrixMean, matrix with mean counts,
#'      OR
#'      meanValue for one context/Condition.
#'
#'  @examples
#'     # calculate mean counts in contrast contx1_vs_contx2
#'     mean1<-NormCountsMean(glmfit, ASKOlist, contx1)   # in the 1st context
#'     mean2<-NormCountsMean(glmfit, ASKOlist, contx2)   # in the 2nd context
#'
#'     # for all conditions
#'     n_count<-NormCountsMean(glmfit, ASKOlist, context=NULL)
#'
#'  @export
NormCountsMean <- function(glmfit, ASKOlist, context=NULL){
  lib_size_norm<-glmfit$samples$lib.size*glmfit$samples$norm.factors                          # normalization computation of all library sizes
  if(is.null(context)==T){
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

    if(is.null(context)==T){
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
#' @example
#'    AskoStats(glm_test, fit, contrast, ASKOlist, dge, parameters)
#'
#'  @export
AskoStats <- function (glm_test, fit, contrast, ASKOlist, dge, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  asko_dir = paste0(study_dir, "Askomics/")
  image_dir = paste0(study_dir, "images/")

  contrasko<-ASKOlist$contrast$Contrast[row.names(ASKOlist$contrast)==contrast]   # to retrieve the name of contrast from Asko object
  contx1<-ASKOlist$contrast$context1[row.names(ASKOlist$contrast)==contrast]      # to retrieve the name of 1st context from Asko object
  contx2<-ASKOlist$contrast$context2[row.names(ASKOlist$contrast)==contrast]      # to retrieve the name of 2nd context from Asko object

  ASKO_stat<-glm_test$table
  ASKO_stat$Test_id<-paste(contrasko, rownames(ASKO_stat), sep = "_")             # addition of Test_id column = unique ID
  ASKO_stat$contrast<-contrasko                                                   # addition of the contrast of the test
  ASKO_stat$gene <- row.names(ASKO_stat)                                          # addition of gene column = gene ID
  ASKO_stat$FDR<-p.adjust(ASKO_stat$PValue, method=parameters$p_adj_method)       # computation of False Discovery Rate

  # Between context1 and context2 :
  ASKO_stat$Significance=0
  ASKO_stat$Significance[ASKO_stat$logFC < -parameters$threshold_logFC & ASKO_stat$FDR <= parameters$threshold_FDR] = -1  # Significance values = -1 for down regulated genes
  ASKO_stat$Significance[ASKO_stat$logFC > parameters$threshold_logFC  & ASKO_stat$FDR <= parameters$threshold_FDR] = 1   # Significance values =  1 for up regulated genes

  # addition of column "expression"
  ASKO_stat$Expression=NA
  ASKO_stat$Expression[ASKO_stat$Significance==-1]<-paste(contx1, contx2, sep="<")  # the value of attribute "Expression" is a string
  ASKO_stat$Expression[ASKO_stat$Significance==1]<-paste(contx1, contx2, sep=">")   # this attribute is easier to read the Significance
  ASKO_stat$Expression[ASKO_stat$Significance==0]<-paste(contx1, contx2, sep="=")   # of expression between two contexts

  if(parameters$Expression==TRUE){colg="Expression"}else{colg=NULL}
  if(parameters$logFC==T){cola="logFC"}else{cola=NULL}
  # computation of Fold Change from log2FC
  if(parameters$FC==T){colb="FC";ASKO_stat$FC <- 2^abs(ASKO_stat$logFC)}else{colb=NULL}
  if(parameters$Sign==T){colc="Significance"}
  if(parameters$logCPM==T){cold="logCPM"}else{cold=NULL}
  if(parameters$LR==T){cole="LR"}else{cole=NULL}
  if(parameters$FDR==T){colf="FDR"}else{colf=NULL}

  # adding table "stat.table" to the ASKOlist
  ASKOlist$stat.table<-ASKO_stat[,c("Test_id","contrast","gene",cola,colb,"PValue",colg,colc,cold,cole,colf)]

  if(parameters$mean_counts==T){                            # computation of the mean of normalized counts for conditions
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

  write.table(ASKOlist$stat.table,paste0(asko_dir, parameters$organism, "_", contrasko, ".txt"), sep=parameters$sep, col.names = T, row.names = F, quote=FALSE)

  # for image size
  nsamples<-ncol(dge$counts)
  sizeImg=15*nsamples
  if(sizeImg < 480) {sizeImg=480}

  # heatmap of Most Differential Genes Expression
  if(parameters$heatmap==TRUE){
    cpm_gstats<-cpm(dge, log=TRUE)[o,][seq(parameters$numhigh),]
    png(paste0(image_dir, contrast, "_topDGE_heatmap.png"), width=sizeImg*1.5, height=sizeImg*1.5)
    par(oma=c(2,2,2,2))
    heatmap.2(cpm_gstats,
              trace="none",
              scale="row",
              labCol=dge$samples$Name,
              main = contrasko,
              xlab = "samples",
              ColSideColors = dge$samples$color,
              margins = c(12,12),
              Rowv = FALSE,
              dendrogram="col")
    dev.off()
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
#' @example
#'    sum_table<-DEanalysis(norm_GE, data_list, asko_list, parameters)
#'
#' @export
DEanalysis <- function(norm_GE, data_list, asko_list, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "images/")

  # for image size
  nsamples <- ncol(data_list$dge$counts)
  sizeImg=15*nsamples
  if(sizeImg < 480){ sizeImg=480 }

  # Checks Contrasts
  if(is.null(parameters$select_sample) & is.null(parameters$rm_sample)){
    c1<-levels(data_list$samples$condition)
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
  png(paste0(image_dir, parameters$analysis_name, "_biological_coefficient_of_variation.png"), width=sizeImg, height=sizeImg)
  plotBCV(normGEdisp)
  dev.off()

  # Genewise Negative Binomial Generalized Linear Models
  if(parameters$glm=="lrt"){
    fit <- glmFit(normGEdisp, data_list$design, robust = T)
  }
  # Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests
  else if(parameters$glm=="qlf"){
    fit <- glmQLFit(normGEdisp, data_list$design, robust = T)
    png(paste0(image_dir, parameters$analysis_name, "_quasi-likelihood_dispersion.png"), width=sizeImg, height=sizeImg)
    plotQLDisp(fit)
    dev.off()
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
    sum[,colnames(contrast)]<-decideTestsDGE(glm_test, adjust.method = parameters$p_adj_method, lfc=parameters$threshold_logFC)
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
      sum[,contrast]<-decideTestsDGE(glm_test, adjust.method = parameters$p_adj_method, lfc=parameters$threshold_logFC)
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
  ctime<-format(Sys.time(), "%d-%m-%Y_%Hh%Mm%Ss")
  sumFile<-paste0(study_dir,parameters$analysis_name,"_summary_DE_",ctime,".csv")
  if(is.null(data_list$annot)==FALSE)
  {
    rnames<-row.names(sum)                        # get Genes DE names
    annDE<-as.matrix(data_list$annot[rnames,])    # get annotations for each genes DE
    rownames(annDE)<-rnames
    colnames(annDE)<-colnames(data_list$annot)
    SumMat<-cbind(sum,annDE)                      # merge the two matrix

    write.table(SumMat, file=sumFile, col.names=NA, row.names = T, quote=F, sep='\t')
  }
  else
  {
    write.table(sum, file=sumFile, col.names=NA, row.names = T, quote=F, sep='\t')
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
#' @param sumDEG, data frame contains for each contrast the significance expression (1/0/-1) for all gene.
#' @param data_list, list contain all data and metadata (DGEList, samples descritions, contrast, design and annotations).
#' @param parameters, list that contains all arguments charged in Asko_start.
#'
#' @example
#'    UpSetGraph(sumDEG, data_list, parameters)
#'
#' @export
UpSetGraph <- function(resDEG, data_list, parameters){
  options(warn = -1)
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  image_dir = paste0(study_dir, "UpSetR_graphs/")
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
        png(paste0(global_dir, parameters$analysis_name,"_UpSetR_allDEG.png"), width=1280, height=1024)
        print(upset(data=ordDEG, sets=rev(sets), nsets=ncol(ordDEG), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA, text.scale = tsc))
        grid.text("All differentially expressed genes (up+down)", x=0.65, y=0.95, gp=gpar(fontsize=26))
        dev.off()
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
        png(paste0(global_dir, parameters$analysis_name,"_UpSetR_upDEG.png"), width=1280, height=1024)
        print(upset(data=ordDEG, sets=rev(sets), nsets=ncol(ordDEG), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA, text.scale = tsc))
        grid.text("Genes expressed \"UP\"", x=0.65, y=0.95, gp=gpar(fontsize=26))
        dev.off()
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
        png(paste0(global_dir, parameters$analysis_name,"_UpSetR_downDEG.png"), width=1280, height=1024)
        print(upset(data=downDEG, sets=rev(colnames(downDEG)), nsets=ncol(downDEG), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA, text.scale = tsc))
        grid.text("Genes expressed \"DOWN\"", x=0.65, y=0.95, gp=gpar(fontsize=26))
        dev.off()
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
        png(paste0(global_dir, parameters$analysis_name,"_UpSetR_mixedDEG.png"), width=1280, height=1024)
        print(upset(data=ordDEG, sets=rev(sets), nsets=ncol(ordDEG), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA,
              text.scale = tsc, set.metadata = list(data = metadata, plots = list(list(type = "matrix_rows",
              column = "SENS", colors = c(UP = "#FF9999", DOWN = "#99FF99"), alpha = 0.5)))))
        grid.text("Genes expressed \"UP\" and \"DOWN\"", x=0.65, y=0.95, gp=gpar(fontsize=26))
        dev.off()
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
      compa<-as.vector(strsplit2(comparaison, "-"))
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
          png(paste0(subset_dir, parameters$analysis_name,"_UpSetR_",comparaison,"_allDEG.png"), width=1280, height=1024)
          print(upset(data=ordDEG, sets=rev(compa), nsets=length(compa), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA, text.scale = tsc))
          grid.text("All differentially expressed genes (up+down)", x=0.65, y=0.95, gp=gpar(fontsize=26))
          dev.off()
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
          png(paste0(subset_dir, parameters$analysis_name,"_UpSetR_",comparaison,"_upDEG.png"), width=1280, height=1024)
          print(upset(data=ordDEG, sets=rev(compa), nsets=length(compa), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA, text.scale = tsc))
          grid.text("Genes expressed \"UP\"", x=0.65, y=0.95, gp=gpar(fontsize=26))
          dev.off()
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
          png(paste0(subset_dir, parameters$analysis_name,"_UpSetR_",comparaison,"_downDEG.png"), width=1280, height=1024)
          print(upset(data=ordDEG, sets=rev(newcompa), nsets=length(newcompa), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA, text.scale = tsc))
          grid.text("Genes expressed \"DOWN\"", x=0.65, y=0.95, gp=gpar(fontsize=26))
          dev.off()
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
          png(paste0(subset_dir, parameters$analysis_name,"_UpSetR_",comparaison,"_mixedDEG.png"), width=1280, height=1024)
          print(upset(data=ordDEG, sets=rev(sets), nsets=length(sets), keep.order=TRUE, sets.bar.color="#56B4E9", nintersects=NA,
                      text.scale = tsc, set.metadata = list(data = metadata,
                      plots = list(list(type = "matrix_rows",column = "SENS", colors = c(UP = "#FF9999", DOWN = "#99FF99"), alpha = 0.5)))))
          grid.text("Genes expressed \"UP\" and \"DOWN\"", x=0.65, y=0.95, gp=gpar(fontsize=26))
          dev.off()
        }
      }
    }
  }
}

#' @title VD
#'
#' @description Plot Venn Diagram to compare different contrast
#'
#' @param sumDEG, data frame contains for each contrast the significance expression (1/0/-1) for all gene.
#' @param asko_list, list of data.frame contain condition, contrast and context informations made by asko3c.
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @return none.
#'
#' @example
#'    VD(sumDEG, parameters, asko_list)
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
  venn_dir = paste0(study_dir, "vennDiagram/")
  if(dir.exists(venn_dir)==FALSE){
    dir.create(venn_dir)
    cat("\n\nDirectory: ",venn_dir," created\n")
  }

  # don't write log file
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

  cat("\nCreate VennDiagrams ")
  if(parameters$VD == "all"){
    cat("for all differentially expressed genes.\n")
    for(comparaison in parameters$compaVD){
      compa<-strsplit2(comparaison, "-")
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

      color <- brewer.pal(nbCompa, parameters$palette)
      venn.diagram(input,
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
      compa<-strsplit2(comparaison, "-")
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
      venn<-venn.diagram(input, main="Genes expressed \"UP\" and \"DOWN\"",
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
      compa<-strsplit2(comparaison, "-")
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

      color <- brewer.pal(nbCompa, parameters$palette)
      venn.diagram(input,
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
      compa<-strsplit2(comparaison, "-")
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

      color <- brewer.pal(nbCompa, parameters$palette)
      venn.diagram(input,
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
#' @param sumDEG, data frame contains for each contrast the significance expression (1/0/-1) for all gene.
#' @param parameters, list that contains all arguments charged in Asko_start.
#' @return none.
#'
#' @examples
#'    GOenrichment(sumDEG, parameters)
#'
#' @export
GOenrichment<-function(resDEG, data, parameters){
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

  for(contrast in colnames(data$contrast)){
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
      GOdata <- new("topGOdata",
                    nodeSize = parameters$GO_min_num_genes,
                    ontology = ontology,
                    allGenes = geneList,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)
      print(GOdata)

      resultTest <- runTest(GOdata, algorithm = parameters$GO_algo, statistic = parameters$GO_stats)
      print(resultTest)

      resGenTab <- GenTable(GOdata, numChar = 1000000, statisticTest = resultTest, orderBy = "statisticTest", topNodes=length(nodes(graph(GOdata))) )
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
    write.table(TabCompl, file=paste0(img_go_dir, parameters$analysis_name, "_", contrast, "_Complet_GOenrichment.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

    if(exists("TabSigCompl")==TRUE){
      if(nrow(TabSigCompl)>=1){
        if (parameters$GO_max_top_terms > 10) {
          TabSigCompl$Term = str_trunc(TabSigCompl$Term, 67)
        }else{
          TabSigCompl$Term = str_trunc(TabSigCompl$Term, 137)
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
          scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) + theme_linedraw() +
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
          scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) + theme_linedraw() +
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

