
Asko_start <- function(){
  # Loading libraries in silent mode (only error messages will be displayed) 
  pkgs<-c("limma","statmod","edgeR","VennDiagram","RColorBrewer",
          "ggplot2","ggrepel","gplots","stringr","optparse","goSTAG")
  for(p in pkgs) suppressPackageStartupMessages(library(p, quietly=TRUE, character.only=TRUE))
  
  # Options List
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
                help="col of ids in count files [default= %default]", metavar="integer"),
    make_option(c("-C", "--col_counts"), type="integer", default=7,dest="col_counts",
                help="col of counts in count files [default= %default (featureCounts) ]", metavar="integer"),
    make_option(c("-t", "--sep"), type="character", default="\t", dest="sep",
                help="col separator [default= %default]", metavar="character"),
    make_option(c("--GO_MF"), type="character", default=NULL, dest="GO_MF",
                help="GO Molecular Fonction annotation file [default= %default]", metavar="character"),
    make_option(c("--GO_BP"), type="character", default=NULL, dest="GO_BP",
                help="GO Biological Process annotation file [default= %default]", metavar="character"),
    make_option(c("--GO_CC"), type="character", default=NULL, dest="GO_CC",
                help="GO Cellular Component annotation file [default= %default]", metavar="character"),
    make_option(c("-a", "--annotation"), type="character", default=NULL, dest="annotation",
                help="annotation file [default= %default]", metavar="character"),
    make_option(c("-s", "--sample"), type="character", default="Samples.txt", dest="sample_file",
                help="Samples file [default= %default]", metavar="character"),
    make_option(c("-c", "--contrasts"), type="character", default="Contrasts.txt",dest="contrast_file",
                help="Contrasts file [default= %default]", metavar="character"),
    make_option(c("-k", "--mk_context"), type="logical", default=FALSE,dest="mk_context",
                help="generate automatically the context names [default= %default]", metavar="logical"),
    make_option(c("-p", "--palette"), type="character", default="Set2", dest="palette",
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
                help=" GLM method (lrt/qlf) [default= %default]", metavar="character"),
    make_option(c("--lfc"), type="logical", default="TRUE", dest="logFC",
                help="logFC in the summary table [default= %default]", metavar="logical"),
    make_option(c("--th_lfc"), type="double", default=1, dest="threshold_logFC",
                help="logFC threshold [default= %default]", metavar="double"),
    make_option("--fc", type="logical", default="TRUE", dest="FC",
                help="FC in the summary table [default= %default]", metavar="logical"),
    make_option(c("--lcpm"), type="logical", default="FALSE", dest="logCPM",
                help="logCPm in the summary table [default= %default]", metavar="logical"),
    make_option("--fdr", type="logical", default="TRUE", dest="FDR",
                help="FDR in the summary table [default= %default]", metavar="logical"),
    make_option("--lr", type="logical", default="FALSE", dest="LR",
                help="LR in the summary table [default= %default]", metavar="logical"),
    make_option(c("--sign"), type="logical", default="TRUE", dest="Sign",
                help="Significance (1/0/-1) in the summary table [default= %default]", metavar="logical"),
    make_option(c("--expr"), type="logical", default="TRUE", dest="Expression",
                help="Significance expression in the summary table [default= %default]", metavar="logical"),
    make_option(c("--mc"), type="logical", default="TRUE", dest="mean_counts",
                help="Mean counts in the summary table [default= %default]", metavar="logical"),
    make_option(c("--hm"), type="logical", default="TRUE", dest="heatmap",
                help="generation of the expression heatmap [default= %default]", metavar="logical"),
    make_option(c("--nh"), type="integer", default="50", dest="numhigh",
                 help="number of genes in the heatmap [default= %default]", metavar="integer"),
    make_option(c("--VD"), type = "character", default = NULL, dest = "VD",
                help = "", metavar = ""),
    make_option(c("--compaVD"), type = "character", default = NULL, dest = "compaVD",
                help = "", metavar = ""),
    make_option(c("--GO"), type="character", default=NULL, dest="GO",
                help = "gene set chosen for GO enrichment analysis 'up', 'down', 'both', or NULL", metavar="character"),
    make_option(c("--GO_filt_meth"), type="character", default = "p.adjust", dest="GO_filt_meth",
                help = "Use 'pval' to filter on nominal p-value or 'p.adjust' to filter on adjusted p-value", metavar="character"),
    make_option(c("--GO_padj_meth"), type="character", default = "BH", dest="GO_padj_meth",
                help = "correction method used to adjust p-values; available option : 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'", metavar = "character"),
    make_option(c("--GO_threshold"), type="numeric", default = "0.05", dest="GO_threshold",
                help = "the significant threshold used to filter p-values", metavar = "integer"),
    make_option(c("--GO_min_num_terms"), type="integer", default = "10", dest="GO_min_num_terms",
                help = "the minimum number of GO terms required to plot a cluster label", metavar = "integer")
  )
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

loadData <- function(parameters){

  # Dossiers pour les fichiers de sortie
  #---------------------------------------------------------
  # initialise le répertoire de travail
  cat("\n\nCreate directories:\n")
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/") 
  if(dir.exists(study_dir)==FALSE){ dir.create(study_dir) }
  cat("\t",study_dir,"\n")
  
  # initialise les noms des dossiers et les crée si inexistant
  image_dir = paste0(study_dir, "images/") 
  if(dir.exists(image_dir)==FALSE){ dir.create(image_dir) }
  cat("\t",image_dir,"\n")
  
  venn_dir = paste0(study_dir, "vennDiagram/") 
  if(dir.exists(venn_dir)==FALSE){ dir.create(venn_dir) }
  cat("\t",venn_dir,"\n")
  
  asko_dir = paste0(study_dir, "Askomics/") 
  if(dir.exists(asko_dir)==FALSE){ dir.create(asko_dir) }
  cat("\t",asko_dir,"\n")
  
  # Gestions des fichiers d'entrée 
  #---------------------------------------------------------
  # répertoire contenant les fichiers d'entrée
  input_path = paste0(parameters$dir_path, "/input/") 
  
  # fichier décrivant les échantillons
  sample_path<-paste0(input_path, parameters$sample_file)
  samples<-read.csv(sample_path, header=TRUE, sep="\t", row.names=1, comment.char = "#")

  # Selection de certains sample (parametre select sample)
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

  # Suppression de samples (parametres rm_sample)
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
  
  # décrit les conditions et les couleurs pour chaque échantillons
  condition<-unique(samples$condition)
  if(length(condition)<3){ color=c("#FF9999","#99CCFF") }
  else{ color<-brewer.pal(length(condition), parameters$palette) }
  samples$color<-NA
  j=0
  for(name in condition){
    j=j+1
    samples$color[samples$condition==name]<-color[j]
  }

  # Fichier(s) de comptages.
  # Deux possibilités : 1 fichier de comptage par condition 
  # ou une matrice de comptage pour toutes les conditions
  #---------------------------------------------------------
  # fichiers de comptages (1 par conditions)
  if(is.null(parameters$fileofcount)){
    cat("\nFiles of counts:\n")
    print(samples$file)
    cat("\nSamples:\n")
    print(rownames(samples))
    
    # creates a DGEList object from a table of counts
    dge<-readDGE(paste0(input_path,samples$file), labels=rownames(samples), columns=c(parameters$col_genes,parameters$col_counts), header=TRUE, comment.char="#")
    dge<-DGEList(counts=dge$counts, samples=samples)
  }
  # matrix file with all counts for all conditions
  else{
    cat("\nSamples:\n")
    print(rownames(samples))
    count_path<-paste0(input_path, parameters$fileofcount)
    if(grepl(".csv", parameters$fileofcount)==TRUE){
      count<-read.csv(count_path, header=TRUE, sep = "\t", row.names = parameters$col_genes)
    }
    else{
      count<-read.table(count_path, header=TRUE, sep = "\t", row.names = parameters$col_genes)
    }
    
    # if you ask for some samples were removed for analysis
    select_counts<-row.names(samples)
    countT<-count[,select_counts]
    
    # creates a DGEList object from a table of counts
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

  # verify if some colunms will be not use for analysis
  rmcol<-list()
  for(condition_name in row.names(contrastab)){
    test<-match(condition_name, colnames(designExp),nomatch = 0)
    if(test==0){
      rm<-grep("0", contrastab[condition_name,], invert = T)
      if(is.null(rmcol)){rmcol=rm}else{rmcol<-append(rmcol, rm)}
    }
  }
  # if it's the case then it delete them
  if (length(rmcol)> 0){
    rmcol<-unlist(rmcol)
    rmcol<-unique(rmcol)
    rmcol=-rmcol
    contrastab<-contrastab[,rmcol]
  }
  # sort contrast table if more than one contrast in contrastab
  if(length(contrastab)>1){
    ord<-match(colnames(designExp),row.names(contrastab), nomatch = 0)
    contrast_table<-contrastab[ord,]
  }
  else{
    contrast_table<-contrastab
  }
  cat("\nContrasts:\n")
  print(contrast_table)

  # format contrast : convert "+" to "1" and - to "-1"
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

  # Annotations and GOterms files
  #---------------------------------------------------------
  # annnotation file
  if(is.null(parameters$annotation)==FALSE){
    annot<-read.csv(paste0(input_path, parameters$annotation), header = T, row.names = 1, sep = '\t', quote = "")
    data[["annot"]]=annot
  }
  # GO terms for molecular function
  if(is.null(parameters$GO_MF)==FALSE) {
    goMF<-read.csv(paste0(input_path, parameters$GO_MF), header = F, sep = '\t', quote = "")
    data[["GO_MF"]]=goMF
  }
  # GO terms for biological process
  if(is.null(parameters$GO_BP)==FALSE) {
    goBP<-read.csv(paste0(input_path, parameters$GO_BP), header = F, sep = '\t', quote = "")
    data[["GO_BP"]]=goBP
  }
  # GO terms for cellular component
  if(is.null(parameters$GO_CC)==FALSE) {
    goCC<-read.csv(paste0(input_path, parameters$GO_CC), header = F, sep = '\t', quote = "")
    data[["GO_CC"]]=goCC
  }

  return(data)
}

asko3c <- function(data_list){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")
  asko_dir = paste0(study_dir, "Askomics/")
  asko<-list()

  ##### Condition #####
  
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

  ##### Contrast + Context #####  
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
      for(j in 1:length(set_cond2)){                                            # for each condition contained in the complex context (2nd):
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
      for(j in 1:length(set_cond1)){                                            # for each condition contained in the complex context (1st):
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
          #print(value)
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
      for(j in 1:length(set_cond1)){                                            # for each condition contained in the complex context (1st):
        list_context<-append(list_context, contx1)                              # verification of the presence of values in each condition
        list_condition<-append(list_condition, set_cond1[j])                    # contained in the 1st context
      }
      for(j in 1:length(set_cond2)){                                            # for each condition contained in the complex context (2nd):
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
  
  ##### Files creation #####
  ctime<-format(Sys.time(), "%d-%m-%Y_%Hh%Mm%Ss")
  write.table(data.frame("Condition"=row.names(condition_asko),condition_asko), 
              paste0(asko_dir,"condition.asko",ctime,".txt"), 
              sep = parameters$sep, 
              row.names = F, 
              quote=F)            # creation of condition file for asko 
  write.table(context_asko,  
              paste0(asko_dir, "context.asko",ctime,".txt"), 
              sep=parameters$sep, 
              col.names = T, 
              row.names = F,
              quote=F)            # creation of context file for asko
  write.table(contrast_asko,  
              paste0(asko_dir, "contrast.asko",ctime,".txt"), 
              sep=parameters$sep, 
              col.names = T, 
              row.names = F, 
              quote=F)            # creation of contrast file for asko
  return(asko)
}

GEfilt <- function(data_list, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/") 
  image_dir = paste0(study_dir, "images/")
  
##### plot density before filtering #####
  cpm<-cpm(data_list$dge)
  logcpm<-cpm(data_list$dge, log=TRUE)
  
  colnames(logcpm)<-rownames(data_list$dge$samples)
  nsamples <- ncol(data_list$dge$counts)                                                                    # creation nouveau plot 
  png(paste0(image_dir, parameters$analysis_name, "_raw_data.png"))
  plot(density(logcpm[,1]),
       col=as.character(data_list$dge$samples$color[1]),      # plot exprimant la densit? de chaque g?ne   
       lwd=1,
       ylim=c(0,0.21),
       las=2,
       main="A. Raw data",
       xlab="Log-cpm")        # en fonction de leurs valeurs d'expression
  abline(v=0, lty=3)
  for (i in 2:nsamples){                                                        # on boucle sur chaque condition restante
    den<-density(logcpm[,i])                                                    # et les courbes sont rajout?es dans le plot
    lines(den$x, col=as.character(data_list$dge$samples$color[i]), den$y, lwd=1)   #
  }
  legend("topright", rownames(data_list$dge$samples), 
         text.col=as.character(data_list$dge$samples$color), 
         bty="n",
         text.width=6,
         cex=0.5)
  dev.off()
#####plot density after filtering#####                         # rowSums compte le nombre de score (cases) pour chaque colonne Sup ? 0.5
  keep.exprs <- rowSums(cpm>parameters$threshold_cpm)>=parameters$replicate_cpm      # en ajoutant >=3 cela donne un test conditionnel
  filtered_counts <- data_list$dge[keep.exprs,,keep.lib.sizes=F] # si le comptage respecte la condition alors renvoie TRUE
  filtered_cpm<-cpm(filtered_counts$counts, log=TRUE)
  png(paste0(image_dir ,parameters$analysis_name,"_filtered_data.png"))
  plot(density(filtered_cpm[,1]),
       col=as.character(data_list$dge$samples$color[1]),
       lwd=2,
       ylim=c(0,0.21),
       las=2,
       main="B. Filtered data", xlab="Log-cpm")
  abline(v=0, lty=3) 
  for (i in 2:nsamples){
    den <- density(filtered_cpm[,i])
    lines(den$x,col=as.character(data_list$dge$samples$color[i]), den$y, lwd=1)
  } 
  legend("topright", rownames(data_list$dge$samples),
         text.col=as.character(data_list$dge$samples$col),
         bty="n",
         text.width=6,
         cex=0.5)
  dev.off()
##### histogram cpm values distribution before filtering #####
  png(paste0(image_dir,parameters$analysis_name,"_read_count_distrib_before_filtering.png"))
  countNotNull <- data_list$dge$counts[-which(rowSums(data_list$dge$counts) == 0),]
  #hist(log2(countNotNull+1),
  hist(logcpm,
       main= "A. Log2(cpm) distribution before filtering",
       xlab = "log2(cpm)",
       col = "grey")
  dev.off()
##### histogram cpm values distribution after filtering #####
  png(paste0(image_dir,parameters$analysis_name,"_read_count_distrib_after_filtering.png"))
  #hist(log2(filtered_counts$counts+1), 
  hist(filtered_cpm,
       main= "B. Log2(cpm) distribution after filtering",
       #xlab = "log(counts +1)",
       xlab = "log2(cpm)",
       col = "grey")
  dev.off()
##### boxplot cpm values distribution before filtering ##### 
  png(paste0(image_dir,parameters$analysis_name,"_boxplot_read_count_distrib1.png"))
  #boxplot(log(countNotNull+1),
  boxplot(logcpm,
          col=data_list$dge$samples$color,        
          main="A. Distribution of log2(cpm) values per sample before filtering",
          cex.axis=0.5,
          las=2,
          ylab="log2(cpm)")
  dev.off()
##### boxplot cpm values distribution after filtering #####
  png(paste0(image_dir,parameters$analysis_name,"_boxplot_read_count_distrib2.png"))
  #boxplot(log(filtered_counts$counts+1),
  boxplot(filtered_cpm,
          col=data_list$dge$samples$color,       
          main="B. Distribution of log2(cpm) values per sample after filtering",
          cex.axis=0.5,
          las=2,
          ylab="log2(cpm)")
  dev.off()
##### histogram libsize before normalization #####
  png(paste0(image_dir,parameters$analysis_name,"_barplot_libsize.png"))
  barplot(filtered_counts$sample$lib.size,
          main = "Distribution of library sizes per samples",
          ylab = "Library size",
          cex.axis = 0.6,
          cex.names = 0.7,
          names.arg = rownames(data_list$dge$samples),
          beside = T,
          las = 2,
          col = data_list$dge$samples$color,
          space = 0.3)
  dev.off()
  return(filtered_counts)
}

GEnorm <- function(filtered_GE, parameters){
  
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/") 
  image_dir = paste0(study_dir, "images/")
  
  filtered_cpm <- cpm(filtered_GE, log=TRUE)   #nouveau calcul Cpm sur donnees filtrees, si log=true alors valeurs cpm en log2 
  colnames(filtered_cpm)<-rownames(filtered_GE$samples)
  ##### boxplot log2(cpm) values before normalization #####
  png(paste0(image_dir, parameters$analysis_name,"_boxplot_logcpm_before_norm.png"))
  boxplot(filtered_cpm,
          col=filtered_GE$samples$color,         #boxplot des scores cpm non normalis?s
          main="A. Before normalization",
          cex.axis=0.5,
          las=2,
          ylab="Log-cpm")
  dev.off()
  
  norm_GE<-calcNormFactors(filtered_GE, method = parameters$normal_method)                      # normalisation de nos comptages par le methode TMM, estimation du taux de production d'un ARN                                                                      # en estimant l'?chelle des facteurs entre echantillons -> but : pouvoir comparer nos ech entre eux
  
  ##### boxplot log2(cpm) values after normalization #####
  logcpm_norm <- cpm(norm_GE, log=TRUE)
  colnames(logcpm_norm)<-rownames(filtered_GE$samples)
  png(paste0(image_dir, parameters$analysis_name,"_boxplot_logcpm_after_norm.png"))
  boxplot(logcpm_norm,
          col=filtered_GE$samples$color, 
          main="B. After normalization",
          cex.axis=0.5,
          las=2,
          ylab="Log-cpm")
  dev.off()
  ##### histogram libsize after normalization PAS BON A REFAIRE ##### 
  library.size <- norm_GE$samples$lib.size * norm_GE$samples$norm.factors
  png(paste0(image_dir,parameters$analysis_name,"_barplot_effective_libsize.png"))
  barplot(library.size,
          main = "Distribution of effective library sizes per samples",
          ylab = "Effective library size",
          cex.axis = 0.6,
          cex.names = 0.7,
          names.arg = rownames(norm_GE$samples),
          beside = T,
          las = 2,
          col = norm_GE$samples$color,
          space = 0.3)
  dev.off()
  return(norm_GE)
}

GEcorr <- function(dge, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/") 
  image_dir = paste0(study_dir, "images/")
  
  lcpm<-cpm(dge, log=TRUE)
  colnames(lcpm)<-rownames(dge$samples)
  cormat<-cor(lcpm)
 # color<- colorRampPalette(c("yellow", "white", "green"))(20)
  color<-colorRampPalette(c("black","red","yellow","white"),space="rgb")(35)
  png(paste0(image_dir, parameters$analysis_name, "_heatmap_of_sample_correlation.png"))
  heatmap(cormat, col=color, symm=TRUE,RowSideColors =as.character(dge$samples$color), ColSideColors = as.character(dge$samples$color))
  dev.off()
  #MDS
  mds <- cmdscale(dist(t(lcpm)),k=3, eig=TRUE)
  eigs<-round((mds$eig)*100/sum(mds$eig[mds$eig>0]),2)
  
  png(paste0(image_dir, parameters$analysis_name, "_MDS_corr_axe1_2.png"))
  mds1<-ggplot(as.data.frame(mds$points), aes(V1, V2, label = rownames(mds$points))) + labs(title="MDS Axes 1 and 2") + geom_point(color =as.character(dge$samples$color) ) + xlab(paste('dim 1 [', eigs[1], '%]')) +ylab(paste('dim 2 [', eigs[2], "%]")) + geom_label_repel(aes(label = rownames(mds$points)), color = 'black',size = 3.5)
  print(mds1)
  dev.off()
  #ggsave("mds_corr1-2.tiff")
  #ggtitle("MDS Axes 2 and 3")
  png(paste0(image_dir, parameters$analysis_name, "_MDS_corr_axe2_3.png"))
  mds2<-ggplot(as.data.frame(mds$points), aes(V2, V3, label = rownames(mds$points))) + labs(title="MDS Axes 2 and 3") + geom_point(color =as.character(dge$samples$color) ) + xlab(paste('dim 2 [', eigs[2], '%]')) +ylab(paste('dim 3 [', eigs[3], "%]")) + geom_label_repel(aes(label = rownames(mds$points)), color = 'black',size = 3.5)
  print(mds2)
  dev.off()
  # ggtitle("MDS Axes 1 and 3")
  #ggsave("mds_corr2-3.tiff")
  png(paste0(image_dir, parameters$analysis_name, "_MDS_corr_axe1_3.png"))
  mds3<-ggplot(as.data.frame(mds$points), aes(V1, V3, label = rownames(mds$points))) + labs(title="MDS Axes 1 and 3") + geom_point(color =as.character(dge$samples$color) ) + xlab(paste('dim 1 [', eigs[1], '%]')) +ylab(paste('dim 3 [', eigs[3], "%]")) + geom_label_repel(aes(label = rownames(mds$points)), color = 'black',size = 3.5)
  print(mds3)
  dev.off()
  #ggsave("mds_corr1-3.tiff")
  
  #####hierarchical clustering##### 
  mat.dist <- dist(t(asko_norm$counts), method = "euclidean")
  clustering <- hclust(mat.dist, method = "complete")
  png(paste0(image_dir, parameters$analysis_name, "_hclust.png"))
  plot(clustering,
       main = 'cluster dendrogram',
       hang = -1)
  dev.off()
  
}

.NormCountsMean <- function(glmfit, ASKOlist, context=NULL){
  
  lib_size_norm<-glmfit$samples$lib.size*glmfit$samples$norm.factors                          # normalization computation of all library sizes 
  if(is.null(context)==T){
    set_condi<-row.names(ASKOlist$condition)
  }else{
    set_condi<-ASKOlist$context$condition[ASKOlist$context$context==context]                    # retrieval of condition names associated to context
  }
  table_c_norm <- data.frame(row.names = row.names(glmfit$counts))

  for (condition in set_condi){
    sample_name<-rownames(glmfit$samples[glmfit$samples$condition==condition,])               # retrieval of the replicate names associated to conditions
    subset_counts<-data.frame(row.names = row.names(glmfit$counts))                            # initialization of data frame as subset of counts table
    for(name in sample_name){
      lib_sample_norm<-glmfit$samples[name,"lib.size"]*glmfit$samples[name,"norm.factors"]    # normalization computation of sample library size 
      subset_counts$c<-glmfit$counts[,name]                                                   # addition in subset of sample counts column
      subset_counts$c<-subset_counts$c*mean(lib_size_norm)/lib_sample_norm                    # normalization computation of sample counts
      colnames(subset_counts)[colnames(subset_counts)=="c"]<-name                             # to rename the column with the condition name
    }
    mean_counts<-rowSums(subset_counts)/ncol(subset_counts)                                   # computation of the mean
    #print(length(mean_counts))
    table_c_norm$m <- mean_counts
    #colnames(matrix_c_nor)<-paste(context,condition,sep = "/")             
    if(is.null(context)==T){
      colnames(table_c_norm)[colnames(table_c_norm)=="m"]<-condition          
    }else{
      colnames(table_c_norm)[colnames(table_c_norm)=="m"]<-paste(context,condition,sep = "/")             
    }
    #print(table_c_norm)
    #ASKOlist$stat.table$mean<-mean_counts                                                     # subset integration in the glm_result table 
  }                                                                                           # to rename the column with the context name
  return(table_c_norm)
  #return(ASKOlist$stat.table)                                                                 # return the glm object
}

VD <- function(decideTestTable, parameters, asko_list){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")  #
  venn_dir = paste0(study_dir, "vennDiagram/") 
  
  cat("\nCreate VD graphs\n")
  if(is.null(parameters$VD)==TRUE || parameters$VD=="" || parameters$VD == "all"){
    for(comparaison in parameters$compaVD){
      name<-c()
      title<-c()
      input<-list()
      compa<-strsplit2(comparaison, "-")
      nbCompa <- length(compa)
      print(nbCompa)
      for(n in 1:nbCompa){
        col_num <- which(colnames(decideTestTable)==compa[n])
        na <- paste0(asko_list$contrast$context1[col_num],"/",asko_list$contrast$context2[col_num])
        name <- append(name, na)
        ti <- asko_list$contrast$Contrast[rownames(asko_list$contrast)==compa[n]]
        title <- append(title, ti)
        all <- rownames(decideTestTable)[decideTestTable[,col_num]!=0]
        print(length(all))
        input[[n]] <- append(input, all)
        print(name)
      }
      if(nbCompa==2){
        input_list = list(c1=input[[1]] ,c2=input[[2]])
        color<-c("palegreen","skyblue")
      }
      if(nbCompa==3){input_list = list(c1=input[[1]] ,c2=input[[2]], c3=input[[3]])}
      if(nbCompa==4){input_list = list(c1=input[[1]] ,c2=input[[2]], c3=input[[3]], c4=input[[4]])}
      if(nbCompa==5){input_list = list(c1=input[[1]] ,c2=input[[2]], c3=input[[3]], c4=input[[4]], c5=input[[5]])}
      color <- brewer.pal(nbCompa, parameters$palette)
      title_file <- paste(title, sep = "-", collapse = "-")
      filename <- paste0(title_file,"all")
      venn.diagram(input_list, main=title_file,
                   filename=paste0(venn_dir, filename, ".png"),
                   height = 3000,
                   width = 3500,
                   imagetype = "png",
                   cat.cex = 1.5,
                   cex = 1.5,
                   # cat.pos = c(340, 20, 180),
                   fill = color,
                   category.names = name,
                   col=0,euler.d = FALSE,scaled=FALSE
      )
    }
  }
  else if(parameters$VD == "both"){
    for(comparaison in parameters$compaVD){
      compa<-strsplit2(comparaison, "-")
      column1<-which(colnames(decideTestTable)==compa[1])
      column2<-which(colnames(decideTestTable)==compa[2])
  
      na1<-paste0(asko_list$contrast$context1[column1],"<",asko_list$contrast$context2[column1])
      na2<-paste0(asko_list$contrast$context1[column1],">",asko_list$contrast$context2[column1])
      na3<-paste0(asko_list$contrast$context1[column2],"<",asko_list$contrast$context2[column2])
      na4<-paste0(asko_list$contrast$context1[column2],">",asko_list$contrast$context2[column2])
      
      name_c1<-asko_list$contrast$Contrast[rownames(asko_list$contrast)==compa[1]]
      name_c2<-asko_list$contrast$Contrast[rownames(asko_list$contrast)==compa[2]]
      
      Gup_c1<-rownames(decideTestTable)[decideTestTable[,column1]==-1]
      Gdown_c1<-rownames(decideTestTable)[decideTestTable[,column1]==1]
      Gup_c2<-rownames(decideTestTable)[decideTestTable[,column2]==-1]
      Gdown_c2<-rownames(decideTestTable)[decideTestTable[,column2]==1]
      input<-list(up_1=Gup_c1,
                  down_1=Gdown_c1,
                  up_2=Gup_c2,
                  down_2=Gdown_c2)
      
      filename = paste0(name_c1,"-",name_c2)
      print(filename)
      venn<-venn.diagram(input, main=paste(name_c1, name_c2, sep = "/"),
                   filename=paste0(venn_dir, filename, ".png"),
                   height = 3000,
                   width = 3000,
                   imagetype = "png",
                   main.cex = 1.5,
                   cat.cex = 1.2,
                   cex=1.5,
                   cat.dist = c(-0.4,-0.4,0.1,0.1),
                   cat.col = c( "red1","royalblue1", "red3", "royalblue4"),
                   category.names = c(na1, na2, na3, na4),
                   col=c( "red1","royalblue1", "red3", "royalblue4"),
                   euler.d = FALSE,
                   scaled=FALSE)
    }
  }
  else if(parameters$VD == "up"){
    for(comparaison in parameters$compaVD){
      name<-c()
      title<-c()
      input<-list()
      compa<-strsplit2(comparaison, "-")
      nbCompa <- length(compa)
      print(nbCompa)
      for(n in 1:nbCompa){
        col_num <- which(colnames(decideTestTable)==compa[n])
        na <- paste0(asko_list$contrast$context1[col_num],"<",asko_list$contrast$context2[col_num])
        name <- append(name, na)
        
        ti <- asko_list$contrast$Contrast[rownames(asko_list$contrast)==compa[n]]
        title <- append(title, ti)
        up <- rownames(decideTestTable)[decideTestTable[,col_num]==-1]
        print(length(up))
        input[[n]] <- append(input, up)
        
      }
      if(nbCompa==2){
        input_list = list(c1=input[[1]] ,c2=input[[2]])
        color<-c("palegreen","skyblue")
      }
      if(nbCompa==3){input_list = list(c1=input[[1]] ,c2=input[[2]], c3=input[[3]])}
      if(nbCompa==4){input_list = list(c1=input[[1]] ,c2=input[[2]], c3=input[[3]], c4=input[[4]])}
      if(nbCompa==5){input_list = list(c1=input[[1]] ,c2=input[[2]], c3=input[[3]], c4=input[[4]], c5=input[[5]])}
      color <- brewer.pal(nbCompa, parameters$palette)
      title_file <- paste(title, sep = "-", collapse = "-")
      filename <- paste0(title_file,"up")
      print(name)
      venn.diagram(input_list, main=title_file,
                   filename=paste0(venn_dir, filename, ".png"),
                   imagetype = "png",
                   fill = color,
                   cex=1.5,
                   cat.cex=1,
                   category.names = name,
                   col=0,euler.d = FALSE,scaled=FALSE
                   )
    }
  }
  else if(parameters$VD == "down"){
    for(comparaison in parameters$compaVD){
      name<-c()
      title<-c()
      input<-list()
      compa<-strsplit2(comparaison, "-")
      nbCompa <- length(compa)
      for(n in 1:nbCompa){
        col_num <- which(colnames(decideTestTable)==compa[n])
        na <- paste0(asko_list$contrast$context1[col_num],">",asko_list$contrast$context2[col_num])
        name <- append(name, na)
        ti <- asko_list$contrast$Contrast[rownames(asko_list$contrast)==compa[n]]
        title <- append(title, ti)
        down <- rownames(decideTestTable)[decideTestTable[,col_num]==1]
        print(length(down))
        input[[n]] <- append(input, down)
        
      }
      if(nbCompa==2){
        input_list = list(c1=input[[1]] ,c2=input[[2]])
        color<-c("palegreen","skyblue")
      }
      if(nbCompa==3){input_list = list(c1=input[[1]] ,c2=input[[2]], c3=input[[3]])}
      if(nbCompa==4){input_list = list(c1=input[[1]] ,c2=input[[2]], c3=input[[3]], c4=input[[4]])}
      if(nbCompa==5){input_list = list(c1=input[[1]] ,c2=input[[2]], c3=input[[3]], c4=input[[4]], c5=input[[5]])}
      color <- brewer.pal(nbCompa, parameters$palette)
      title_file <- paste(title, sep = "-", collapse = "-")
      filename <- paste0(title_file,"down")
      venn.diagram(input_list, main=title_file,
                   filename=paste0(venn_dir, filename, ".png"),
                   imagetype = "png",
                   fill = color,
                   cex=1.5,
                   cat.cex=1,
                   category.names = name,
                   col=0,euler.d = FALSE,scaled=FALSE
      )
    }
  }
  cat("\n")
}

AskoStats <- function (glm_test, fit, contrast, ASKOlist, dge, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/") 
  asko_dir = paste0(study_dir, "Askomics/")
  image_dir = paste0(study_dir, "images/")

  contrasko<-ASKOlist$contrast$Contrast[row.names(ASKOlist$contrast)==contrast]         # to retrieve the name of contrast from Asko object
  contx1<-ASKOlist$contrast$context1[row.names(ASKOlist$contrast)==contrast]            # to retrieve the name of 1st context from Asko object 
  contx2<-ASKOlist$contrast$context2[row.names(ASKOlist$contrast)==contrast]            # to retrieve the name of 2nd context from Asko object

  ASKO_stat<-glm_test$table
  ASKO_stat$Test_id<-paste(contrasko, rownames(ASKO_stat), sep = "_")                   # addition of Test_id column = unique ID
  ASKO_stat$contrast<-contrasko                                                         # addition of the contrast of the test
  ASKO_stat$gene <- row.names(ASKO_stat)                                                # addition of gene column = gene ID
  ASKO_stat$FDR<-p.adjust(ASKO_stat$PValue, method=parameters$p_adj_method)                                # computation of False Discovery Rate
 
  ASKO_stat$Significance=0                                                              # Between context1 and context2 :
  ASKO_stat$Significance[ASKO_stat$logFC < -parameters$threshold_logFC & ASKO_stat$FDR <= parameters$threshold_FDR] = -1       # Significance values = -1 for down regulated genes
  ASKO_stat$Significance[ASKO_stat$logFC > parameters$threshold_logFC  & ASKO_stat$FDR <= parameters$threshold_FDR] = 1        # Significance values =  1 for up regulated genes

  ASKO_stat$Expression=NA                                                             # addition of column "expression" 
  ASKO_stat$Expression[ASKO_stat$Significance==-1]<-paste(contx1, contx2, sep="<")    # the value of attribute "Expression" is a string
  ASKO_stat$Expression[ASKO_stat$Significance==1]<-paste(contx1, contx2, sep=">")     # this attribute is easier to read the Significance
  ASKO_stat$Expression[ASKO_stat$Significance==0]<-paste(contx1, contx2, sep="=")     # of expression between two contexts
  
  if(parameters$Expression==TRUE){colg="Expression"}else{colg=NULL}
  if(parameters$logFC==T){cola="logFC"}else{cola=NULL}                                             #
  if(parameters$FC==T){colb="FC";ASKO_stat$FC <- 2^abs(ASKO_stat$logFC)}else{colb=NULL}            # computation of Fold Change from log2FC
  if(parameters$Sign==T){colc="Significance"}                                                      #
  if(parameters$logCPM==T){cold="logCPM"}else{cold=NULL}                                           #
  if(parameters$LR==T){cole="LR"}else{cole=NULL}                                                   #
  if(parameters$FDR==T){colf="FDR"}else{colf=NULL}

  # adding table "stat.table" to the ASKOlist  
  ASKOlist$stat.table<-ASKO_stat[,c("Test_id","contrast","gene",cola,colb,"PValue",colg,colc,cold,cole,colf)]
  
  if(parameters$mean_counts==T){                             # computation of the mean of normalized counts for conditions
    mean1<-.NormCountsMean(fit, ASKOlist, contx1)            # in the 1st context
    mean2<-.NormCountsMean(fit, ASKOlist, contx2)            # in the 2nd context
    ASKOlist$stat.table<-cbind(ASKOlist$stat.table, mean1)
    ASKOlist$stat.table<-cbind(ASKOlist$stat.table, mean2)
  }
  print(table(ASKO_stat$Expression))
  colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="gene"] <- paste("is", "gene", sep="@")                  # header formatting for askomics
  colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="contrast"] <- paste("measured_in", "Contrast", sep="@") # header formatting for askomics
  o <- order(ASKOlist$stat.table$FDR)                                                                                   # ordering genes by FDR value
  ASKOlist$stat.table<-ASKOlist$stat.table[o,]
  
  write.table(ASKOlist$stat.table,paste0(asko_dir, parameters$organism, "_", contrasko, ".txt"),                                    #
              sep=parameters$sep, col.names = T, row.names = F, quote=FALSE)
  
  if(parameters$heatmap==TRUE){
    cpm_gstats<-cpm(dge, log=TRUE)[o,][1:parameters$numhigh,]
    png(paste0(image_dir, contrast, "_topDEG_heatmap.png"))#, height = 3500, width = 3500)
    heatmap.2(cpm_gstats,
              cexRow=0.8, 
              cexCol=0.8,
              lwid = c(2, 8),
              lhei = c(2, 8),
              scale="row", 
              labCol=dge$samples$Name, 
              main = contrasko,
              xlab = "samples",
              ColSideColors = dge$samples$color,
              sepwidth = c(0.05,0.2),
              #colCol = dge$samples$color,
              srtRow = 0,
              adjRow = c(0,0),
              srtCol = 90,
              adjCol = c(1,1),
              margins = c(6,8),
              Rowv = FALSE, 
              dendrogram="col")
    dev.off()  
    }
}


DEanalysis <- function(norm_GE, data_list, asko_list, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/") 
  image_dir = paste0(study_dir, "images/")
  
  # Vérifie que les contrast sont renseigné correctement
  c1<-levels(data$samples$condition)
  len1<-length(c1)
  c2<-rownames(data$contrast)
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
  
  normGEdisp <- estimateDisp(norm_GE, data_list$design)
  if(parameters$glm=="lrt"){
    fit <- glmFit(normGEdisp, data_list$design, robust = T)
  }
  if(parameters$glm=="qlf"){
    fit <- glmQLFit(normGEdisp, data_list$design, robust = T)
    #plotQLDisp(fit)
  }
  
  sum<-data.frame(row.names = rownames(fit))
  if(length(data_list$contrast)==1){
    contrast<-makeContrasts(contrasts = data_list$contrast, levels = data_list$design)
    colnames(contrast)<-colnames(data_list$contrast)
    if(parameters$glm=="lrt"){ 
      glm_test<-glmLRT(fit, contrast=contrast) 
    }
    if(parameters$glm=="qlf"){ 
      glm_test<-glmQLFTest(fit, contrast=contrast)
    }
    sum[,colnames(contrast)]<-decideTestsDGE(glm_test, adjust.method = parameters$p_adj_method, lfc=parameters$threshold_logFC)
    AskoStats(glm_test, fit, colnames(contrast), asko_list,normGEdisp,parameters)
  }
  else{
      for (contrast in colnames(data_list$contrast)){
        if(parameters$glm=="lrt"){
          glm_test<-glmLRT(fit, contrast=data_list$contrast[,contrast])
        }
        if(parameters$glm=="qlf"){
          glm_test<-glmQLFTest(fit, contrast=data_list$contrast[,contrast])
        }
        sum[,contrast]<-decideTestsDGE(glm_test, adjust.method = parameters$p_adj_method, lfc=parameters$threshold_logFC)
        AskoStats(glm_test, fit, contrast, asko_list,normGEdisp,parameters)
      }
  }

  cat("\n\nCreating HeatMap\n\n")
  #####heatmap cpm value per sample #####
  if(nrow(norm_GE) <= 100000){
    cpm_norm <- cpm(norm_GE, log = F)
    dat.n <- scale(t(cpm_norm))
    dim(dat.n)
    dat.tn <- t(dat.n)
    d1 <- dist(dat.n,method = "euclidean", diag = FALSE, upper = FALSE)
    d2 <- dist(dat.tn,method = "euclidean", diag = FALSE, upper = TRUE)
    c1 <- hclust(d1, method = "complete", members = NULL)
    c2 <- hclust(d2, method = "complete", members = NULL)
    my_palette <- colorRampPalette(c("green","black","red"), interpolate = "linear")
    ctime<-format(Sys.time(), "%d-%m-%Y_%Hh%Mm%Ss")
    png(paste0(image_dir, parameters$analysis_name, "_heatmap_",ctime,".png"))
    heatmap.2(dat.tn,                     # Tidy, normalised data
              Colv = as.dendrogram(c1),     # Experiments clusters in cols
              Rowv = as.dendrogram(c2),     # Protein clusters in rows
              density.info="histogram",   # Plot histogram of data and colour key
              trace = "none",               # Turn of trace lines from heat map
              col = my_palette,           # Use my colour scheme
              cexRow = 0.1,
              cexCol = 1,
              xlab = "samples",
              #lmat = rbind(c(1, 3), c(1,3), c(1,3)),
              #lwid = c(2, 10),
              lhei = c(1, 5, 2),
              sepwidth = c(0, 0),
              ColSideColors = norm_GE$samples$color,
              main = paste0("gene 1 to ",nrow(norm_GE)),
              margins = c(15,0))     # Amend row and column label fonts
    
    dev.off()
  }
  
  #####heatmap per condition#####
  n_count<-.NormCountsMean(fit, ASKOlist = asko_list)
  cpm_condition<-cpm(n_count)
  cat("\n\nCPM conditions\n")
  print(head(cpm_condition))

  cat("\n\nCreate summary file\n\n")
  # create summary file with annotations (if available) and contrast value for each gene
  if(is.null(data_list$annot)==FALSE)
  {
    #---------------------------------
    # Ne gère qu'une colonne ... à modifier pour :
    #          1) un nombre de colonnes de 1 à N 
    #          2) place NA pour les gènes non présents dans le fichier "annotation"
    #*********************************
    # rnames<-row.names(sum)                        # récupère les noms des gènes DE
    # annDE<-as.matrix(data_list$annot[rnames,])    # récupère les annotations correspondante
    # rownames(annDE)<-rnames                       # replace les idgenes en nom de lignes
    # colnames(annDE)<-c("Description")             # annotation défini comme description
    # SumMat<-cbind(annDE,sum)                      # merge les deux matrices
    #---------------------------------

    # print(colnames(sum))
    # print(colnames(data_list$annot)) 
    SumMat<-sum
    # It is ugly but the R merge function uses too much memory 
    for (i in row.names(sum)){
      for (j in colnames(data_list$annot)) {
        SumMat[i,j]<-data_list$annot[i,j]
      }
    }
    ctime<-format(Sys.time(), "%d-%m-%Y_%Hh%Mm%Ss")
    write.table(SumMat, paste0(study_dir,parameters$analysis_name,"_summary_DE_",ctime,".csv"), col.names=T, row.names = T, quote=F, sep='\t')
  }
  else
  {
    write.table(sum, paste0(study_dir,parameters$analysis_name,"_summary_DE_",ctime,".csv"), col.names=T, row.names = T, quote=F, sep='\t')
  }
  
  #return(glm_test) 
  return(sum)
  #VD(sum, parameters, asko_list)
}


#### Enrichment Analysis ####
#############################
# Function : Loop goSTAG 
#---------------------------------------
loopGoStag<-function(gene_list,go_list,lvl,nameGo){
  study_dir = paste0(parameters$dir_path, "/", parameters$analysis_name, "/") 
  image_dir = paste0(study_dir, "images/")
  
  # Generating the Enrichment Matrix
  cat("1st step: create a matrix of GO enrichment scores\n\n")
  enrichment_matrix<-matrix()
  try(enrichment_matrix <- performGOEnrichment(gene_list, go_list,
                                               filter_method = parameters$GO_filt_meth,
                                               significance_threshold = parameters$GO_threshold,
                                               p.adjust_method = parameters$GO_padj_meth) 
  )
  if (nrow(enrichment_matrix)==1) {return(NULL)} 
  # Hierarchical Clustering
  cat("2nd step: cluster the GO terms\n")
  # Il peut prendre plus de parametres comme : distance_method="euclidean",clustering_method="complete" (cf doc goSTAG)
  try(hclust_results <- performHierarchicalClustering(enrichment_matrix)) 
  print(hclust_results)
  
  # Grouping the Clusters
  cat("3rd step: group the GO tems into clusters\n\n")
  try(clusters <- groupClusters(hclust_results))
  print(clusters)
  
  # Annotating the Clusters
  cat("4th step: annotate each of the clusters\n\n")
  try(cluster_labels <- annotateClusters(clusters))
  print(cluster_labels)
  
  # Plotting a Heatmap
  GOtitle=""
  if(nameGo=="MF"){ GOtitle="Molecular Function" }
  if(nameGo=="BP"){ GOtitle="Biological Process" }
  if(nameGo=="CC"){ GOtitle="Cellular Component" }
  png(paste0(image_dir, parameters$analysis_name,"_",nameGo,"_enrich_heatmap_",lvl,".png"), width=1500, height=1200)
  try(plotHeatmap(enrichment_matrix,
                  hclust_results, 
                  clusters,
                  cluster_labels,
                  min_num_terms = parameters$GO_min_num_terms,
                  dendrogram_width=0.5,
                  cluster_label_width=0.8,
                  cluster_label_cex=2,
                  sample_label_cex=3,
                  dendrogram_lwd=0.5,
                  header_lwd=0.5,
                  header_height=0.3, 
                  heatmap_colors = "extra"))
  dev.off()
  
  # save all in matrix
  matrix<-list()
  matrix[[paste0(as.character(nameGo),"_enrich_",as.character(lvl))]] <- enrichment_matrix
  matrix[[paste0(as.character(nameGo),"_clusters_",as.character(lvl))]] <- clusters
  matrix[[paste0(as.character(nameGo),"_label_cluster_",as.character(lvl))]] <- cluster_labels
  matrix[[paste0(as.character(nameGo),"_hclust_results_",as.character(lvl))]] <- hclust_results
  matrix[[paste0(as.character(nameGo),"_gene_list_",as.character(lvl))]] <- gene_list
  matrix[[paste0(as.character(nameGo),"_go_list_",as.character(lvl))]] <- go_list
  return(matrix)
}

# Function : Running goSTAG 
#--------------------------------------- 
runGoStag<-function(summaryDEG, asko_list, data_go, nameGo){
  # classe GO pour les titres
  GOtitle=""
  if(nameGo=="MF"){ GOtitle="Molecular Function" }
  if(nameGo=="BP"){ GOtitle="Biological Process" }
  if(nameGo=="CC"){ GOtitle="Cellular Component" }
  
  # GO terms lists
  #----------------------------------------
  all_genes <- as.character(unique(data_go[[1]]))
  go_all<-as.character(unique(data_go[[2]]))
  go_list<-list()
  for(n in 1:length(go_all)){
    gene <- data_go[1][data_go[2]==go_all[n]]
    if(length(gene)!=0){
      go_list$g <- as.character(gene)
      Go_name <- as.character(go_all[n])
      names(go_list)[names(go_list)=="g"]<-Go_name
    }
  }
  go_list$all <- as.character(all_genes)
  names(go_list)[names(go_list)=="all"]<-"ALL"
  
  # create Gene DE lists for each contrast
  #----------------------------------------
  gene_listUP<-list()
  gene_listDOWN<-list() 
  matrixUP<-list()
  matrixDOWN<-list()
  
  if(parameters$GO=="both" | parameters$GO=="up"){
    # retrieve data and place it in list of genes by contrast
    for(n in 1:ncol(summaryDEG)){
      # contrast
      contrast_name<-colnames(summaryDEG[n])    
      contrastDE_name <- asko_list$contrast$Contrast[rownames(asko_list$contrast)==contrast_name]
      # all genes up for this contrast
      DEGup<-rownames(summaryDEG[n])[summaryDEG[n]==-1] 
      if(length(DEGup)!=0){
        gene_listUP$up<-DEGup
        names(gene_listUP)[names(gene_listUP)=="up"]<-contrastDE_name
      }
    }
    # Make enrichment analysis for DEGenes
    cat("\n")
    cat(GOtitle,": Enrichment Analysis for DEGenes UP\n")
    matrixUP<-loopGoStag(gene_listUP,go_list,"up",nameGo)
  }
  
  if(parameters$GO=="both" | parameters$GO=="down"){
    # retrieve data and place it in list of genes by contrast
    for(n in 1:ncol(summaryDEG)){
      # contrast
      contrast_name<-colnames(summaryDEG[n])    
      contrastDE_name <- asko_list$contrast$Contrast[rownames(asko_list$contrast)==contrast_name]
      # all genes up for this contrast
      DEGdown<-rownames(summaryDEG[n])[summaryDEG[n]==1] 
      if(length(DEGdown)!=0){
        gene_listDOWN$down<-DEGdown
        names(gene_listDOWN)[names(gene_listDOWN)=="down"]<-contrastDE_name
      }
    }
    # Make enrichment analysis for DEGenes
    cat("\n")
    cat(GOtitle,": Enrichment Analysis for DEGenes DOWN\n")
    matrixDOWN<-loopGoStag(gene_listDOWN,go_list,"down",nameGo)
  }
  
  # save and return all results
  matrixAll<-append(matrixUP,matrixDOWN)
  return(matrixAll)
}
