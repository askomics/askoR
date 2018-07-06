#askor_path<-"/home/flegeai/workspace/askoR"
askor_path<-"/home/flegeai/local/askoR"
source(paste0(askor_path,"/AskoR.R"))
#source("~/")
##############################################
##                Parameters                ##  
##############################################
setwd("/home/flegeai/workspace/askoR/Test")
parameters<-Asko_start()

# parameters$analysis_name = "" 
# parameters$dir_path = "/"
# parameters$organism = ""
parameters$fileofcount = "mrna.counts.txt"
#parameters$annotation_file = "annotation.txt"
parameters$sample_file = "Samples.txt"
parameters$contrast_file = "Contrasts.txt"
# parameters$sep = "\t"
# parameters$col_genes = 1
# parameters$col_counts = 5
# parameters$select_sample = c("")                                        #
# parameters$rm_sample = c("")                                       #
# parameters$regex = T
# parameters$mk_context = T
# 
# parameters$threshold_cpm = 0.5                                                  #
# parameters$replicate_cpm= 3                                                     #
# parameters$threshold_FDR = 0.05                                                 #
# parameters$normal_method = "TMM"                                                #
# parameters$p_adj_method = "fdr"                                                 #
# parameters$glm = "lrt"                                                          #
# parameters$logFC = T
# parameters$FC = T
# parameters$logCPM = F
# parameters$FDR = T
# parameters$LR = F
# parameters$Sign = T
# parameters$Expression = T
# parameters$mean_counts = T
# 
#parameters$palette = "Set3"
# parameters$heatmap = T
# parameters$numhigh = 50
# 
# parameters$VD = "both"
# parameters$compaVD=c("")
# 
# parameters$GSEA = "both"                                                        #
# parameters$GSEA_filt_meth = "p.adjust"                                          #
# parameters$GSEA_padj_meth = "BH"                                                #
# parameters$GSEA_threshold = 0.05                                                 #
# parameters$GSEA_min_num_terms = 5                                               #

########################################
##  Loading the data from the samples ##
########################################
#####load data#####
data<-loadData(parameters)

getdata$samples
data$contrast
data$design
data$dge$counts

cat("Total number of genes : ", dim(data$dge$counts)[1], "\n")
cat("Total number of samples : ", dim(data$dge$counts)[2], "\n")
cat("summary of CPM by samples\n")
summary(cpm(data$dge))
#####asko files#####
asko_data<-asko3c(data)
asko_data$condition
asko_data$contrast
asko_data$context
#####filtering#####
cat("Filtering genes with more than ", parameters$threshold_cpm, " CPM in ",parameters$replicate_cpm,"samples\n")
asko_filt<-GEfilt(data, parameters)
cat("Total number of filtered genes : ", dim(asko_filt$counts)[1], "\n")
#####normalization#####
asko_norm<-GEnorm(asko_filt,parameters)
#####correlation#####
GEcorr(asko_norm,parameters)
#####DEG analysis#####
cat("Statistical analysis\n")
sum_table<-DEanalysis(asko_norm, data, asko_data, parameters)
#####Venn diagram#####
VD(sum_table, parameters, asko_data)
#####Enrichment Analysis#####
mat<-GSEA(summaryDEG = sum_table, asko_list = asko_data, data_list = data)

