asko3c <- function(data_list){
  asko<-list()
  
  ######### Condition ############ 
  
  condition<-unique(data_list$samples$condition)                                                 # retrieval of different condition's names
  col1<-which(colnames(data_list$samples)=="condition")                                          # determination of number of the column "condition"
  colcol<-which(colnames(data_list$samples)=="color") 
  if(is.null(parameters$fileofcount)){
    col2<-which(colnames(data_list$samples)=="file")                                          # determination of number of the column "replicate"
    column_name<-colnames(data_list$samples[,c(-col1,-col2,-colcol)])    # retrieval of column names needful to create the file condition
  }else{column_name<-colnames(data_list$samples[,c(-col1,-colcol)])}
  condition_asko<-data.frame(row.names=condition)                                           # initialization of the condition's data frame
  #level<-list()                                                                             # initialization of the list will contain the level
                                                                                            # of each experimental factor
  for (name in column_name){                                                                # for each experimental factor :
    # if(str_detect(name, "condition")){                                                      # for the column of conditions, the level is fixed to 0 because
    #   level<-append(level, 0)                                                               # "condition" must be the first column of the data frame
    # }else{                                                                                  #
    #   level<-append(level, length(levels(data_list$samples[,name])))                             # adding to the list the level of other experimental factors
    # }
    # 
    condition_asko$n<-NA                                                                    # initialization of new column in the condition's data frame
    colnames(condition_asko)[colnames(condition_asko)=="n"]<-name                           # to rename the new column with with the name of experimental factor
    for(condition_name in condition){                                                       # for each condition's names
      condition_asko[condition_name,name]<-as.character(unique(data_list$samples[data_list$samples$condition==condition_name, name]))
    }                                                                                       # filling the condition's data frame
  }
  # order_level<-order(unlist(level))                                                         # list to vector
  # condition_asko<-condition_asko[,order_level]                                              # order columns according to their level
  #asko$condition<-condition_asko                                                            # adding data frame of conditions to asko object
  
  #print(condition_asko)
  
  
  #############contrast + context##################  
  i=0
  
  contrast_asko<-data.frame(row.names = colnames(data_list$contrast))           # initialization of the contrast's data frame
  contrast_asko$Contrast<-NA                                                    # all columns are created et initialized with
  contrast_asko$context1<-NA                                                    # NA values
  contrast_asko$context2<-NA                                                    #
  
  list_context<-list()                                                          # initialization of context and condition lists 
  list_condition<-list()                                                        # will be used to create the context data frame
  if(parameters$mk_context==TRUE){
    for (contrast in colnames(data_list$contrast)){                               # for each contrast :
    i=i+1                                                                       # contrast data frame will be filled line by line
    #print(contrast)
    set_cond1<-row.names(data_list$contrast)[data_list$contrast[,contrast]>0]  # retrieval of 1st set of condition's names implicated in a given contrast
    set_cond2<-row.names(data_list$contrast)[data_list$contrast[,contrast]<0] # retrieval of 2nd set of condition's names implicated in a given contrast
    parameters<-colnames(condition_asko)                                        # retrieval of names of experimental factor
    print(paste("set_cond1 : ", set_cond1, sep = ""))
    # print(length(set_cond1))
    print(paste("set_cond2 : ", set_cond2, sep = ""))
    # print(length(set_cond2))
    if(length(set_cond1)==1){complex1=F}else{complex1=T}# to determine if we have complex contrast (multiple conditions
    if(length(set_cond2)==1){complex2=F}else{complex2=T}# compared to multiple conditions) or not
    #print(complex1)
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
      for (param_names in parameters){                                          # for each experimental factor 
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
      for (param_names in parameters){                                          # for each experimental factor:
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
        contx<-str_replace_all(contx, "NULL", "")}else{contx<-common_factor}   # and all common factor are concatenated to become the name of context
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
      for (param_names in parameters){                                          # for each experimental factor:
        print(w)
        w=w+1
        facteur<-unique(c(condition_asko[,param_names]))                        # retrieval of possible values for the experimental factor
        print(paste("facteur : ", facteur, sep=""))
        for(value in facteur){                                                  # for each possible values:
          print(value)
          #print(class(value))
          #print(set_cond1)
          verif1<-unique(str_detect(set_cond1, value))                          # verification of the presence of values in each condition 
                                                                                # contained in the 1st context
          verif2<-unique(str_detect(set_cond2, value))                          # verification of the presence of values in each condition
                                                                                # contained in the 2nd context
          
          if(length(verif1)==1 && verif1==TRUE){m=m+1;common_factor1[m]<-value} # if verif=only TRUE, value of experimental factor is added as common factor
          if(length(verif2)==1 && verif2==TRUE){n=n+1;common_factor2[n]<-value} # if verif=only TRUE, value of experimental factor is added as common factor
        }
      }
      print(paste("common_factor1 : ",common_factor1,sep=""))
      print(paste("common_factor2 : ",common_factor2,sep=""))
      
      if(length(common_factor1)>1){                                             # if there are several common factor for conditions in the 1st context 
        common_factor1<-toString(common_factor1)                                # conversion list to string
        contx1<-str_replace(common_factor1,", ","")}else{contx1<-common_factor1}# all common factor are concatenated to become the name of context
      contx1<-str_replace_all(contx1, "NULL", "")
      print(paste("contx1 : ", contx1, sep=""))
      if(length(common_factor2)>1){                                             # if there are several common factor for conditions in the 2nd context 
        common_factor2<-toString(common_factor2)                                # conversion list to string
        contx2<-str_replace(common_factor2,", ","")}else{contx2<-common_factor2}# all common factor are concatenated to become the name of context
      contx2<-str_replace_all(contx2, "NULL", "")
      print(paste("contx2 : ", contx2, sep=""))
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
#        print(contexts[1])
 #       print(cond1)
        list_context<-append(list_context, contexts[1])
        list_condition<-append(list_condition, cond1)
      }
      for (cond2 in set_cond2){
        list_context<-append(list_context, contexts[2])
        list_condition<-append(list_condition, cond2)
      }
    }
  }
  
  list_context<-unlist(list_context)   # conversion list to vector
  list_condition<-unlist(list_condition)                                                                    # conversion list to vector
#  print(list_condition)
#  print(list_context)
  context_asko<-data.frame(list_context,list_condition)                                                     # creation of the context data frame 
  context_asko<-unique(context_asko)
  colnames(context_asko)[colnames(context_asko)=="list_context"]<-"context"                                 # header formatting for askomics
  colnames(context_asko)[colnames(context_asko)=="list_condition"]<-"condition"                             # header formatting for askomics
  #asko$contrast<-contrast_asko                                                                              # adding context data frame to asko object
  #asko$context<-context_asko                                                                                # adding context data frame to asko object
  asko<-list("condition"=condition_asko, "contrast"=contrast_asko, "context"=context_asko)
  colnames(context_asko)[colnames(context_asko)=="context"]<-"Context"                                      # header formatting for askomics
  colnames(context_asko)[colnames(context_asko)=="condition"]<-"has@Condition"                              # header formatting for askomics
  colnames(contrast_asko)[colnames(contrast_asko)=="context1"]<-paste("context1_of", "Context", sep="@")    # header formatting for askomics
  colnames(contrast_asko)[colnames(contrast_asko)=="context2"]<-paste("context2_of", "Context", sep="@")    # header formatting for askomics
  
  ######## Files creation ########
  
  write.table(condition_asko, paste0(parameters$out_dir,"/condition.asko.txt"), sep = parameters$sep, row.names = F, quote=F)            # creation of condition file for asko 
  write.table(context_asko,  paste0(parameters$out_dir,"/context.asko.txt"), sep=parameters$sep, col.names = T, row.names = F,quote=F)            # creation of context file for asko
  write.table(contrast_asko,  paste0(parameters$out_dir,"/contrast.asko.txt"), sep=parameters$sep, col.names = T, row.names = F, quote=F)          # creation of contrast file for asko
  return(asko)
}

.NormCountsMean <- function(glmfit, ASKOlist, context){
  
  lib_size_norm<-glmfit$samples$lib.size*glmfit$samples$norm.factors                          # normalization computation of all library sizes 
  set_condi<-ASKOlist$context$condition[ASKOlist$context$context==context]                    # retrieval of condition names associated to context
  
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
    ASKOlist$stat.table$mean<-mean_counts                                                     # subset integration in the glm_result table 
    colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="mean"]<-paste(context,condition,sep = "/")             
  }                                                                                           # to rename the column with the context name
  return(ASKOlist$stat.table)                                                                 # return the glm object
}

AskoStats <- function (glm_test, fit, contrast, ASKOlist, dge,parameters){   
  contrasko<-ASKOlist$contrast$Contrast[row.names(ASKOlist$contrast)==contrast]         # to retrieve the name of contrast from Asko object
  contx1<-ASKOlist$contrast$context1[row.names(ASKOlist$contrast)==contrast]            # to retrieve the name of 1st context from Asko object 
  contx2<-ASKOlist$contrast$context2[row.names(ASKOlist$contrast)==contrast]            # to retrieve the name of 2nd context from Asko object
  
  ASKO_stat<-glm_test$table
  ASKO_stat$Test_id<-paste(contrasko, rownames(ASKO_stat), sep = "_")                   # addition of Test_id column = unique ID
  ASKO_stat$contrast<-contrasko                                                         # addition of the contrast of the test
  ASKO_stat$gene <- row.names(ASKO_stat)                                                # addition of gene column = gene ID
  ASKO_stat$FDR<-p.adjust(ASKO_stat$PValue, method=parameters$p_adj_method)                                # computation of False Discovery Rate
  
  ASKO_stat$Significance=0                                                              # Between context1 and context2 :
  ASKO_stat$Significance[ASKO_stat$logFC< 0 & ASKO_stat$FDR<=parameters$threshold_FDR] = -1       # Significance values = -1 for down regulated genes
  ASKO_stat$Significance[ASKO_stat$logFC> 0 & ASKO_stat$FDR<=parameters$threshold_FDR] = 1         # Significance values = 1 for up regulated genes
  
  if(parameters$Expression==TRUE){
    ASKO_stat$Expression=NA                                                             # addition of column "expression" 
    ASKO_stat$Expression[ASKO_stat$Significance==-1]<-paste(contx1, contx2, sep="<")    # the value of attribute "Expression" is a string
    ASKO_stat$Expression[ASKO_stat$Significance==1]<-paste(contx1, contx2, sep=">")     # this attribute is easier to read the Significance
    ASKO_stat$Expression[ASKO_stat$Significance==0]<-paste(contx1, contx2, sep="=")     # of expression between two contexts
  }
  if(parameters$logFC==T){cola="logFC"}else{cola=NULL}                                             #
  if(parameters$FC==T){colb="FC";ASKO_stat$FC <- 2^abs(ASKO_stat$logFC)}else{colb=NULL}            # computation of Fold Change from log2FC
  if(parameters$Sign==T){colc="Significance"}                                                      #
  if(parameters$logCPM==T){cold="logCPM"}else{cold=NULL}                                           #
  if(parameters$LR==T){cole="LR"}else{cole=NULL}                                                   #
  if(parameters$FDR==T){colf="FDR"}else{colf=NULL}
  
  ASKOlist$stat.table<-ASKO_stat[,c("Test_id","contrast","gene",cola,colb,"PValue",     # adding table "stat.table" to the ASKOlist
                                    "Expression",colc,cold,cole,colf)]
  if(parameters$mean_counts==T){                                                                   # computation of the mean of normalized counts for conditions
    ASKOlist$stat.table<-.NormCountsMean(fit, ASKOlist, contx1)                       # in the 1st context
    ASKOlist$stat.table<-.NormCountsMean(fit, ASKOlist, contx2)                       # in the 2nd context
  }
  print(table(ASKO_stat$Expression))
  colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="gene"] <- paste("is", "gene", sep="@")                  # header formatting for askomics
  colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="contrast"] <- paste("measured_in", "Contrast", sep="@") # header formatting for askomics
  o <- order(ASKOlist$stat.table$FDR)                                                                                   # ordering genes by FDR value
  ASKOlist$stat.table<-ASKOlist$stat.table[o,]                                                                          #
  write.table(ASKOlist$stat.table,paste0(parameters$out_dir,"/", parameters$organism, contrasko, ".txt"),                                    #
              sep=parameters$sep, col.names = T, row.names = F, quote=FALSE)
  
  if(parameters$heatmap==TRUE){
    cpm_gstats<-cpm(dge, log=TRUE)[o,][1:parameters$numhigh,]
    heatmap.2(cpm_gstats, cexRow=0.5, cexCol=0.8, scale="row", labCol=dge$samples$Name, xlab=contrast, Rowv = FALSE, dendrogram="col")
  }
  
  return(ASKOlist)
  
}

loadData <- function(parameters){
  
  #####samples#####
  samples<-read.table(parameters$sample_file, header=TRUE, sep="\t", row.names=1, comment.char = "#")       #prise en compte des r?sultats de T2
  if(is.null(parameters$select_sample)==FALSE){
    if(parameters$regex==TRUE){
      selected<-c()
      for(sel in parameters$select_sample){
        select<-grep(sel, rownames(samples))
        if(is.null(selected)){selected=select}else{selected<-append(selected, select)}
      }
      samples<-samples[selected,]
    }else{samples<-samples[parameters$select_sample,]}
  }
  
  if(is.null(parameters$rm_sample)==FALSE){
    if(parameters$regex==TRUE){
      for(rm in parameters$rm_sample){
        removed<-grep(rm, rownames(samples))
        print(removed)
        if(length(removed!=0)){samples<-samples[-removed,]}
      }
    }else{
      for (rm in parameters$rm_sample) {
        rm2<-match(rm, rownames(samples))
        samples<-samples[-rm2,]
      }
    }
  }
  condition<-unique(samples$condition)
  #print(condition)
  color<-brewer.pal(length(condition), parameters$palette)
  #print(color)
  samples$color<-NA
  j=0
  for(name in condition){
    j=j+1
    samples$color[samples$condition==name]<-color[j]
  }
  #print(samples)
  
  
  #####counts#####
  if(is.null(parameters$fileofcount)){
    dge<-readDGE(samples$file, labels=rownames(samples), columns=c(parameters$col_genes,parameters$col_counts), header=TRUE, comment.char="#")
    dge<-DGEList(counts=dge$counts, samples=samples)
    #  dge$samples=samples
    #countT<-dge$counts
    # if(is.null(parameters$select_sample)==FALSE){
    #   slct<-grep(parameters$select_sample, colnames(countT))
    #   dge$counts<-dge$counts[,slct]
    #   dge$samples<-dge$samples[,slct]
    # }
    # if(is.null(parameters$rm_sample)==FALSE){
    #   rmc<-grep(parameters$rm_count, colnames(dge$counts))
    #   dge$counts<-dge$counts[,-rmc]
    #   print(ncol(dge$counts))
    #   rms<-grep(parameters$rm_sample, row.names(dge$samples))
    #   dge$samples<-dge$samples[-rms,]
    # }
  }else {
    if(grepl(".csv", parameters$fileofcount)==TRUE){
      print(parameters$col_genes)
      count<-read.csv(parameters$fileofcount, header=TRUE, sep = "\t", row.names = parameters$col_genes)
      print(count)
      print(row.names(count))
    }
    if(grepl(".txt", parameters$fileofcount)==TRUE){
      count<-read.table(parameters$fileofcount, header=TRUE, sep = "\t", row.names = parameters$col_genes)
    }
    select_counts<-row.names(samples)
    print(select_counts)
    #countT<-count[,c(parameters$col_counts:length(colnames(count)))]
    countT<-count[,select_counts]
    dge<-DGEList(counts=countT, samples=samples) 
    # if(is.null(parameters$select_sample)==FALSE){
    #   slct<-grep(parameters$select_sample, colnames(countT))
    #   countT<-countT[,slct] 
    # }
    # if(is.null(parameters$rm_count)==FALSE){
    #   rms<-grep(parameters$rm_count, colnames(countT))
    #   #print(rms)
    #   countT<-countT[,-rms]
    #   
    # }
    #print(nrow(samples))
    #print(ncol(countT))
  }

  #####design#####
  Group<-factor(samples$condition)
  designExp<-model.matrix(~0+Group)
  rownames(designExp) <- row.names(samples)
  colnames(designExp) <- levels(Group)

  #####contrast#####
  contrastab<-read.table(parameters$contrast_file, sep="\t", header=TRUE, row.names = 1,  comment.char="#", stringsAsFactors = FALSE)

  rmcol<-list()
  for(condition_name in row.names(contrastab)){
    test<-match(condition_name, colnames(designExp),nomatch = 0)
    if(test==0){
      print(condition_name)
      rm<-grep("0", contrastab[condition_name,], invert = T)
      print(rm)
      if(is.null(rmcol)){rmcol=rm}else{rmcol<-append(rmcol, rm)}
    }
  }
  if (length(rmcol)> 0){
    rmcol<-unlist(rmcol)
    rmcol<-unique(rmcol)
    rmcol=-rmcol
    contrastab<-contrastab[,rmcol]
  }
  
  ord<-match(colnames(designExp),row.names(contrastab), nomatch = 0)
  contrast_table<-contrastab[ord,]
  colnum<-c()
 
  for(contrast in colnames(contrast_table)){
    set_cond1<-row.names(contrast_table)[contrast_table[,contrast]=="+"]
    #print(set_cond1)
    set_cond2<-row.names(contrast_table)[contrast_table[,contrast]=="-"]
    #print(set_cond2)
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

    #####annotation#####
  #annotation <- read.csv(parameters$annotation_file, header = T, sep = '\t', quote = "", row.names = 1)
  
  #data<-list("counts"=countT, "samples"=samples, "contrast"=contrast_table, "annot"=annotation, "design"=designExp)
  #print(countT)
  rownames(dge$samples)<-rownames(samples) # replace the renaming by files              
  data<-list("dge"=dge, "samples"=samples, "contrast"=contrast_table, "design"=designExp)
  return(data)
}

GEfilt <- function(dge_list, parameters){
  cpm<-cpm(dge_list)
  logcpm<-cpm(dge_list, log=TRUE)
  colnames(logcpm)<-rownames(dge_list$samples)
  nsamples <- ncol(dge_list)                                                                    # cr?ation nouveau plot 
  plot(density(logcpm[,1]),
       col=as.character(dge_list$samples$color[1]),      # plot exprimant la densit? de chaque g?ne   
       lwd=1,
       ylim=c(0,0.21),
       las=2,
       main="A. Raw data",
       xlab="Log-cpm")        # en fonction de leurs valeurs d'expression
  abline(v=0, lty=3)
  for (i in 2:nsamples){                                                        # on boucle sur chaque condition restante
    den<-density(logcpm[,i])                                                    # et les courbes sont rajout?es dans le plot
    lines(den$x, col=as.character(dge_list$samples$color[i]), den$y, lwd=1)   #
  }
  legend("topright", rownames(dge_list$samples), 
         text.col=as.character(dge_list$samples$color), 
         bty="n",
         text.width=6,
         cex=0.5)
                                                          # rowSums compte le nombre de score (cases) pour chaque colonne Sup ? 0.5
  keep.exprs <- rowSums(cpm>parameters$threshold_cpm)>=parameters$replicate_cpm      # en ajoutant >=3 cela donne un test conditionnel
  filtered_counts <- dge_list[keep.exprs,,keep.lib.sizes=F] # si le comptage respecte la condition alors renvoie TRUE
  filtered_cpm<-cpm(filtered_counts$counts, log=TRUE)

  plot(density(filtered_cpm[,1]),
       col=as.character(dge_list$samples$color[1]),
       lwd=2,
       ylim=c(0,0.21),
       las=2,
       main="B. Filtered data", xlab="Log-cpm")
  abline(v=0, lty=3) 
  for (i in 2:nsamples){
    den <- density(filtered_cpm[,i])
    lines(den$x,col=as.character(dge_list$samples$color[i]), den$y, lwd=1)
  } 
  legend("topright", rownames(dge_list$samples),
         text.col=as.character(dge_list$samples$col),
         bty="n",
         text.width=6,
         cex=0.5)
  return(filtered_counts)
}

GEnorm <- function(filtered_GE, parameters){
  filtered_cpm <- cpm(filtered_GE, log=TRUE)   #nouveau calcul Cpm sur donn?es filtr?es, si log=true alors valeurs cpm en log2 
  colnames(filtered_cpm)<-rownames(filtered_GE$samples)
  boxplot(filtered_cpm,
          col=filtered_GE$samples$color,         #boxplot des scores cpm non normalis?s
          main="A. Before normalization",
          cex.axis=0.5,
          las=2,
          ylab="Log-cpm")
  
  norm_GE<-calcNormFactors(filtered_GE, method = parameters$normal_method)                      # normalisation de nos comptages par le methode TMM, estimation du taux de production d'un ARN                                                                      # en estimant l'?chelle des facteurs entre echantillons -> but : pouvoir comparer nos ech entre eux
  
  logcpm_norm <- cpm(norm_GE, log=TRUE)
  colnames(logcpm_norm)<-rownames(filtered_GE$samples)
  boxplot(logcpm_norm,
          col=filtered_GE$samples$color, 
          main="B. After normalization",
          cex.axis=0.5,
          las=2,
          ylab="Log-cpm")

  return(norm_GE)
}

GEcorr <- function(dge, parameters){
  lcpm<-cpm(dge, log=TRUE)
  colnames(lcpm)<-rownames(dge$samples)
  cormat<-cor(lcpm)
 # color<- colorRampPalette(c("yellow", "white", "green"))(20)
  color<-colorRampPalette(c("black","red","yellow","white"),space="rgb")(28)
  heatmap(cormat, col=color, symm=TRUE,RowSideColors =as.character(dge$samples$color), ColSideColors = as.character(dge$samples$color))
  #MDS
  mds <- cmdscale(dist(t(lcpm)),k=3, eig=TRUE)
  eigs<-round((mds$eig)*100/sum(mds$eig[mds$eig>0]),2)

  mds1<-ggplot(as.data.frame(mds$points), aes(V1, V2, label = rownames(mds$points))) + labs(title="MDS Axes 1 and 2") + geom_point(color =as.character(dge$samples$color) ) + xlab(paste('dim 1 [', eigs[1], '%]')) +ylab(paste('dim 2 [', eigs[2], "%]")) + geom_label_repel(aes(label = rownames(mds$points)), color = 'black',size = 3.5)
  print(mds1)
  #ggsave("mds_corr1-2.tiff")
  #ggtitle("MDS Axes 2 and 3")
  mds2<-ggplot(as.data.frame(mds$points), aes(V2, V3, label = rownames(mds$points))) + labs(title="MDS Axes 2 and 3") + geom_point(color =as.character(dge$samples$color) ) + xlab(paste('dim 2 [', eigs[2], '%]')) +ylab(paste('dim 3 [', eigs[3], "%]")) + geom_label_repel(aes(label = rownames(mds$points)), color = 'black',size = 3.5)
  print(mds2)
  # ggtitle("MDS Axes 1 and 3")
  #ggsave("mds_corr2-3.tiff")
  mds3<-ggplot(as.data.frame(mds$points), aes(V1, V3, label = rownames(mds$points))) + labs(title="MDS Axes 1 and 3") + geom_point(color =as.character(dge$samples$color) ) + xlab(paste('dim 1 [', eigs[1], '%]')) +ylab(paste('dim 3 [', eigs[3], "%]")) + geom_label_repel(aes(label = rownames(mds$points)), color = 'black',size = 3.5)
  print(mds3)
  #ggsave("mds_corr1-3.tiff")
}
  
DEanalysis <- function(norm_GE, data_list, asko_list, parameters){
  
  normGEdisp <- estimateDisp(norm_GE, data_list$design)
  if(parameters$glm=="lrt"){
    fit <- glmFit(normGEdisp, data_list$design, robust = T)
  }
  if(parameters$glm=="qlf"){
    fit <- glmQLFit(normGEdisp, data_list$design, robust = T)
    plotQLDisp(fit)
  }
  
  #plotMD.DGEGLM(fit)     
  #plotBCV(norm_GE)
  
  #sum<-norm_GE$genes
  for (contrast in colnames(data_list$contrast)){
    print(asko_list$contrast$Contrast[contrast])
    if(parameters$glm=="lrt"){
      glm_test<-glmLRT(fit, contrast=data_list$contrast[,contrast])
    }
    if(parameters$glm=="qlf"){
      glm_test<-glmQLFTest(fit, contrast=data_list$contrast[,contrast])
    }
    
    #sum[,contrast]<-decideTestsDGE(glm, adjust.method = parameters$p_adj_method, lfc=1)
    #print(table(sum[,contrast]))
    AskoStats(glm_test, fit, contrast, asko_list,normGEdisp,parameters)
  }
} 

Asko_start <-function(){
  library(limma)
  library(statmod)
  library(edgeR)
  library(ggplot2)
  library(RColorBrewer)
  library(ggrepel)
  library(gplots)
  library(stringr)
  library(optparse)
  option_list = list(
    make_option(c("-o", "--out"), type="character", default="out.pdf",dest="output_pdf",
                help="output file name [default= %default]", metavar="character"),
    make_option(c("-d", "--dir"), type="character", default=".",dest="dir_path",
                help="data directory path [default= %default]", metavar="character"),
    make_option("--outdir", type="character", default=".",dest="out_dir",
                help="outputs directory [default= %default]", metavar="character"),
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
    make_option(c("-a", "--annotation"), type="character", default="annotation.txt", dest="annotation_file",
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
                 help="number of genes in the heatmap [default= %default]", metavar="integer")
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
