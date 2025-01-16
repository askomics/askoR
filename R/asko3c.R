#' @title asko3c
#'
#' @description Create contrast/condition/context file in format readable by Askomics Software.
#'
#' @param data_list list contain all data and metadata (DGEList, samples descriptions, contrast, design and annotations)
#' @param parameters list that contains all arguments charged in Asko_start
#' @return list of data.frame contain condition, contrast and context information
#'
#' @examples
#' \dontrun{
#'    asko<-asko3c(data, parameters)
#' }
#'
#' @note Remember to read the Wiki section in \url{https://github.com/askomics/askoR/wiki}
#' @export
asko3c <- function(data_list, parameters){
  study_dir = paste0(parameters$dir_path,"/", parameters$analysis_name, "/")

  # Askomics directory
  asko_dir = paste0(study_dir, "DEanalysis/AskoTables/")
  asko<-list()

  # Condition
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
  i=0
  contrast_asko<-data.frame(row.names = colnames(data_list$contrast))           # initialization of the contrast's data frame
  contrast_asko$Contrast<-NA                                                    # all columns are created et initialized with
  contrast_asko$context1<-NA                                                    # NA values
  contrast_asko$context2<-NA                                                    #

  list_context<-list()                                                          # initialization of context and condition lists
  list_condition<-list()                                                        # will be used to create the context data frame
  # if(parameters$mk_context==TRUE){
  #   for (contrast in colnames(data_list$contrast)){                             # for each contrast :
  #     i=i+1                                                                       # contrast data frame will be filled line by line
  #     set_cond1<-row.names(data_list$contrast)[data_list$contrast[,contrast]>0]   # retrieval of 1st set of condition's names implicated in a given contrast
  #     set_cond2<-row.names(data_list$contrast)[data_list$contrast[,contrast]<0]   # retrieval of 2nd set of condition's names implicated in a given contrast
  #     set_condition<-colnames(condition_asko)                                     # retrieval of names of experimental factor
  #
  #
  #     if(length(set_cond1)==1){complex1=FALSE}else{complex1=TRUE}                        # to determine if we have complex contrast (multiple conditions
  #     if(length(set_cond2)==1){complex2=FALSE}else{complex2=TRUE}                        # compared to multiple conditions) or not
  #     if(complex1==FALSE && complex2==FALSE){                                             # Case 1: one condition against one condition
  #       contrast_asko[i,"context1"]<-set_cond1                                    # filling contrast data frame with the name of the 1st context
  #       contrast_asko[i,"context2"]<-set_cond2                                    # filling contrast data frame with the name of the 2nd context
  #       contrast_name<-paste(set_cond1,set_cond2, sep = "vs")                     # creation of contrast name by associating the names of contexts
  #       #contrast_name <- paste0(contrast_asko[i,"context1"],"vs",contrast_asko[i,"context2"])
  #       contrast_asko[i,"Contrast"]<-contrast_name                                # filling contrast data frame with contrast name
  #       list_context<-append(list_context, set_cond1)                             #
  #       list_condition<-append(list_condition, set_cond1)                         # adding respectively to the lists "context" and "condition" the context name
  #       list_context<-append(list_context, set_cond2)                             # and the condition name associated
  #       list_condition<-append(list_condition, set_cond2)                         #
  #     }
  #
  #     if(complex1==FALSE && complex2==TRUE){                                             # Case 2: one condition against multiple condition
  #       contrast_asko[i,"context1"]<-set_cond1                                    # filling contrast data frame with the name of the 1st context
  #       list_context<-append(list_context, set_cond1)                             # adding respectively to the lists "context" and "condition" the 1st context
  #       list_condition<-append(list_condition, set_cond1)                         # name and the condition name associated
  #       l=0
  #       # "common_factor" will contain the common experimental factors shared by
  #       common_factor=list()                                                      # conditions belonging to the complex context
  #       for (param_names in set_condition){                                       # for each experimental factor
  #         facteur<-unique(c(condition_asko[,param_names]))                        # retrieval of possible values for the experimental factor
  #         l=l+1                                                                   #
  #         for(value in facteur){                                                  # for each possible values
  #           verif<-unique(stringr::str_detect(set_cond2, value))                           # verification of the presence of values in each condition contained in the set
  #           if(length(verif)==1 && verif==TRUE){common_factor[l]<-value}          # if verif contains only TRUE, value of experimental factor
  #         }                                                                       # is added as common factor
  #       }
  #       if(length(common_factor)>1){                                              # if there are several common factor
  #         common_factor<-toString(common_factor)                                  # the list is converted to string
  #         contx<-stringr::str_replace(common_factor,", ","")
  #         contx<-stringr::str_replace_all(contx, "NULL", "")}else{contx<-common_factor}    # and all common factor are concatenated to become the name of context
  #       contrast_asko[i,"context2"]<-contx                                        # filling contrast data frame with the name of the 2nd context
  #       contrast_name<-paste(set_cond1,contx, sep = "vs")                         # concatenation of context names to make the contrast name
  #       contrast_asko[i,"Contrast"]<-contrast_name                                # filling contrast data frame with the contrast name
  #       for(j in length(set_cond2)){                                            # for each condition contained in the complex context (2nd):
  #         list_context<-append(list_context, contx)                               # adding condition name with the context name associated
  #         list_condition<-append(list_condition, set_cond2[j])                    # to their respective list
  #       }
  #     }
  #     if(complex1==TRUE && complex2==FALSE){                                             # Case 3: multiple conditions against one condition
  #       contrast_asko[i,"context2"]<-set_cond2                                    # filling contrast data frame with the name of the 2nd context
  #       list_context<-append(list_context, set_cond2)                             # adding respectively to the lists "context" and "condition" the 2nd context
  #       list_condition<-append(list_condition, set_cond2)                         # name and the 2nd condition name associated
  #       l=0
  #       # "common_factor" will contain the common experimental factors shared by
  #       common_factor=list()                                                      # conditions belonging to the complex context
  #       for (param_names in set_condition){                                       # for each experimental factor:
  #         facteur<-unique(c(condition_asko[,param_names]))                        # retrieval of possible values for the experimental factor
  #         l=l+1
  #         for(value in facteur){                                                  # for each possible values:
  #           verif<-unique(stringr::str_detect(set_cond1, value))                           # verification of the presence of values in each condition contained in the set
  #           if(length(verif)==1 && verif==TRUE){common_factor[l]<-value}          # if verif contains only TRUE, value of experimental factor
  #         }                                                                       # is added as common factor
  #       }
  #       if(length(common_factor)>1){                                              # if there are several common factor
  #         common_factor<-toString(common_factor)                                  # the list is converted to string
  #         contx<-stringr::str_replace(common_factor,", ","")
  #         contx<-stringr::str_replace_all(contx, "NULL", "")}else{contx<-common_factor}    # and all common factor are concatenated to become the name of context
  #       contrast_asko[i,"context1"]<-contx                                        # filling contrast data frame with the name of the 1st context
  #       contrast_name<-paste(contx,set_cond2, sep = "vs")                         # concatenation of context names to make the contrast name
  #       contrast_asko[i,"Contrast"]<-contrast_name                                # filling contrast data frame with the contrast name
  #       for(j in length(set_cond1)){                                            # for each condition contained in the complex context (1st):
  #         list_context<-append(list_context, contx)                               # adding condition name with the context name associated
  #         list_condition<-append(list_condition, set_cond1[j])                    # to their respective list
  #       }
  #     }
  #     if(complex1==TRUE && complex2==TRUE){                                             # Case 4: multiple conditions against multiple conditions
  #       m=0                                                                       #
  #       n=0                                                                       #
  #       common_factor1=list()                                                     # list of common experimental factors shared by conditions of the 1st context
  #       common_factor2=list()                                                     # list of common experimental factors shared by conditions of the 2nd context
  #       w=1
  #       for (param_names in set_condition){                                       # for each experimental factor:
  #         print(w)
  #         w=w+1
  #         facteur<-unique(c(condition_asko[,param_names]))                        # retrieval of possible values for the experimental factor
  #
  #         for(value in facteur){                                                  # for each possible values:
  #           verif1<-unique(stringr::str_detect(set_cond1, value))                          # verification of the presence of values in each condition contained in the 1st context
  #           verif2<-unique(stringr::str_detect(set_cond2, value))                          # verification of the presence of values in each condition contained in the 2nd context
  #
  #           if(length(verif1)==1 && verif1==TRUE){m=m+1;common_factor1[m]<-value} # if verif=only TRUE, value of experimental factor is added as common factor
  #           if(length(verif2)==1 && verif2==TRUE){n=n+1;common_factor2[n]<-value} # if verif=only TRUE, value of experimental factor is added as common factor
  #         }
  #       }
  #       if(length(common_factor1)>1){                                             # if there are several common factor for conditions in the 1st context
  #         common_factor1<-toString(common_factor1)                                # conversion list to string
  #         contx1<-stringr::str_replace(common_factor1,", ","")}else{contx1<-common_factor1}# all common factor are concatenated to become the name of context
  #       contx1<-stringr::str_replace_all(contx1, "NULL", "")
  #       if(length(common_factor2)>1){                                             # if there are several common factor for conditions in the 2nd context
  #         common_factor2<-toString(common_factor2)                                # conversion list to string
  #         contx2<-stringr::str_replace(common_factor2,", ","")}else{contx2<-common_factor2}# all common factor are concatenated to become the name of context
  #       contx2<-stringr::str_replace_all(contx2, "NULL", "")
  #       contrast_asko[i,"context1"]<-contx1                                       # filling contrast data frame with the name of the 1st context
  #       contrast_asko[i,"context2"]<-contx2                                       # filling contrast data frame with the name of the 2nd context
  #       contrast_asko[i,"Contrast"]<-paste(contx1,contx2, sep = "vs")             # filling contrast data frame with the name of the contrast
  #       for(j in seq_len(set_cond1)){                                            # for each condition contained in the complex context (1st):
  #         list_context<-append(list_context, contx1)                              # verification of the presence of values in each condition
  #         list_condition<-append(list_condition, set_cond1[j])                    # contained in the 1st context
  #       }
  #       for(j in seq_len(set_cond2)){                                            # for each condition contained in the complex context (2nd):
  #         list_context<-append(list_context, contx2)                              # verification of the presence of values in each condition
  #         list_condition<-append(list_condition, set_cond2[j])                    # contained in the 1st context
  #       }
  #     }
  #   }
  # }
  #else{
  for (contrast in colnames(data_list$contrast)){
    i=i+1                                                                       # contrast data frame will be filled line by line
    set_cond1<-row.names(data_list$contrast)[data_list$contrast[,contrast]>0]   # retrieval of 1st set of condition's names implicated in a given contrast
    set_cond2<-row.names(data_list$contrast)[data_list$contrast[,contrast]<0]   # retrieval of 2nd set of condition's names implicated in a given contrast
    set_condition<-colnames(condition_asko)                                     # retrieval of names of experimental factor

    if(length(set_cond1)==1){complex1=FALSE}else{complex1=TRUE}                        # to determine if we have complex contrast (multiple conditions
    if(length(set_cond2)==1){complex2=FALSE}else{complex2=TRUE}                        # compared to multiple conditions) or not


    if(complex1==FALSE && complex2==FALSE){
      contexts = c(set_cond1,set_cond2)
    }

    if(complex1==FALSE && complex2==TRUE){
      pre_contexts=limma::strsplit2(contrast,"vs")
      contexts = c(set_cond1,pre_contexts[2])
    }

    if(complex1==TRUE && complex2==FALSE){
      pre_contexts=limma::strsplit2(contrast,"vs")
      contexts = c(pre_contexts[1],set_cond2)
    }

    if(complex1==TRUE && complex2==TRUE){
      contexts=limma::strsplit2(contrast,"vs")
    }


    if(parameters$projectName!="DEprj" && stringr::str_replace_all(parameters$projectName, " ", "")!=""){
      contrast_asko[i,"Contrast"]<-paste0(parameters$projectName,"_",contrast)
      contrast_asko[i,"Project"]<-parameters$projectName
      contrast_asko[i,"Name"]<-contrast
    }else{
      contrast_asko[i,"Contrast"]<-contrast
    }
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
  # }

  list_context<-unlist(list_context)
  list_condition<-unlist(list_condition)                                                                    # conversion list to vector
  context_asko<-data.frame(list_context,list_condition)                                                     # creation of the context data frame
  context_asko<-unique(context_asko)
  colnames(context_asko)[colnames(context_asko)=="list_context"]<-"context"                                 # header formatting for askomics
  colnames(context_asko)[colnames(context_asko)=="list_condition"]<-"condition"                             # header formatting for askomics
  asko<-list("condition"=condition_asko, "contrast"=contrast_asko, "context"=context_asko)                  # adding context data frame to asko object
  colnames(context_asko)[colnames(context_asko)=="context"]<-"Context"                                      # header formatting for askomics
  colnames(context_asko)[colnames(context_asko)=="condition"]<-"has@Condition"                              # header formatting for askomics
  if(parameters$projectName!="DEprj" && stringr::str_replace_all(parameters$projectName, " ", "")!=""){
    colnames(contrast_asko)[colnames(contrast_asko)=="Project"]<-paste("from", "Project", sep="@")          # header formatting for askomics
  }
  colnames(contrast_asko)[colnames(contrast_asko)=="context1"]<-paste("context1_of", "Context", sep="@")    # header formatting for askomics
  colnames(contrast_asko)[colnames(contrast_asko)=="context2"]<-paste("context2_of", "Context", sep="@")    # header formatting for askomics

  # Files creation
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
