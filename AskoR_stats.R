AskoStats<-function (glm_result, glmfit, constrast, ASKOlist, organism, logFC=F, FC=T, logCPM=F, FDR=T, LR=F, Sign=T, Expression=T, mean_counts=T, csv=F) {   
  
  contrasko<-ASKOlist$contrast$Contrast[row.names(ASKOlist$contrast)==contrast]         # to retrieve the name of contrast from Asko object
  contx1<-ASKOlist$contrast$context1[row.names(ASKOlist$contrast)==contrast]            # to retrieve the name of 1st context from Asko object 
  contx2<-ASKOlist$contrast$context2[row.names(ASKOlist$contrast)==contrast]            # to retrieve the name of 2nd context from Asko object
  
  ASKO_stat<-glm_result$table
  ASKO_stat$Test_id<-paste(contrasko, rownames(ASKO_stat), sep = "_")                   # addition of Test_id column = unique ID
  ASKO_stat$contrast<-contrasko                                                         # addition of the contrast of the test
  ASKO_stat$gene <- row.names(ASKO_stat)                                                # addition of gene column = gene ID
                                                                           
  ASKO_stat$Significance=0                                                              # Between context1 and context2 :
  ASKO_stat$Significance[ASKO_stat$logFC< 0 & ASKO_stat$FDR<=threshold_FDR] = -1        # Significance values = -1 for down regulated genes
  ASKO_stat$Significance[ASKO_stat$logFC> 0 & ASKO_stat$FDR<=threshold_FDR] = 1         # Significance values = 1 for up regulated genes

  if(Expression==TRUE){
    ASKO_stat$Expression=NA                                                             # addition of column "expression" 
    ASKO_stat$Expression[ASKO_stat$Significance==-1]<-paste(contx1, contx2, sep="<")    # the value of attribute "Expression" is a string
    ASKO_stat$Expression[ASKO_stat$Significance==1]<-paste(contx1, contx2, sep=">")     # this attribute is easier to read the Significance
    ASKO_stat$Expression[ASKO_stat$Significance==0]<-paste(contx1, contx2, sep="=")     # of expression between two contexts
  }
  if(logFC==T){cola="logFC"}else{cola=NULL}                                             #
  if(FC==T){colb="FC";ASKO_stat$FC <- 2^abs(ASKO_stat$logFC)}else{colb=NULL}            # computation of Fold Change from log2FC
  if(Sign==T){colc="Significance"}                                                      #
  if(logCPM==T){cold="logCPM"}else{cold=NULL}                                           #
  if(LR==T){cole="LR"}else{cole=NULL}                                                   #
  if(FDR==T){
    colf="FDR";ASKO_stat$FDR<-p.adjust(ASKO_stat$PValue, method="BH")}else{colf=NULL}   # computation of False Discovery Rate

  ASKOlist$stat.table<-ASKO_stat[,c("Test_id","contrast","gene",cola,colb,"PValue",     # adding table "stat.table" to the ASKOlist
                                        "Expression",colc,cold,cole,colf)]
  if(mean_counts==T){                                                                   # computation of the mean of normalized counts for conditions
    ASKOlist$stat.table<-NormCountsMean(glmfit, ASKOlist, contx1)                       # in the 1st context
    ASKOlist$stat.table<-NormCountsMean(glmfit, ASKOlist, contx2)                       # in the 2nd context
    }

  colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="gene"] <- paste("is", "gene", sep="@")                  # header formatting for askomics
  colnames(ASKOlist$stat.table)[colnames(ASKOlist$stat.table)=="contrast"] <- paste("measured_in", "Contrast", sep="@") # header formatting for askomics
  o <- order(ASKOlist$stat.table$FDR)                                                                                   # ordering genes by FDR value
  ASKOlist$stat.table<-ASKOlist$stat.table[o,]                                                                          #
  return(ASKOlist)
  write.table(ASKOlist$stat.table,paste(organism, contrasko, "_test.txt", sep = ""),                                    #
              sep="\t", col.names = T, row.names = F)                                                                   #
  if(csv==T){
    write.csv(ASKOlist$stat.table,paste(organism, contrasko, "_test.txt", sep = ""),                                    #
              sep="\t", col.names = T, row.names = F)
  }
  #write.table(asko, paste("result_for_asko_",organism, contrasko, "_test.txt", sep = ""),                              #
  #            sep="\t", col.names = T, row.names = F)                                                                  #

}
