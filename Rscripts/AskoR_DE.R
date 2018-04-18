askor_path<-Sys.getenv("ASKOR_PATH")
source(paste0(askor_path,"/AskoR.R"))

##############################################
##                Parameters                ##  
##############################################

parameters<-Asko_start()
setwd(parameters$dir_path)
# source("/home/flegeai/local/askoR/askoStart.R")
#parameters$col_genes=1
#parameters$col_counts=7
# parameters$regex=FALSE
#parameters$rm_sample=list("T0_4", "T1K_4", "T1A_4", "T2A_4", "T2K_4", "T3K_4", "T3A_4")
# #parameters$select_sample=c("T0_4", "T1K_4", "T1A_4", "T2A_4", "T2K_4", "T3K_4", "T3A_4")
# #parameters$rm_sample=list("_4")
# parameters$organism = "Ap"
# parameters$fileofcount = NULL
# parameters$annotation_file = "annotation.txt"
# parameters$sample_file = "Samples.txt"
# parameters$contrast_file = "Contrasts.txt"
# parameters$mk_context="manual"
# parameters$glm="qlf"


########################################
##  Loading the data from the samples ##+
########################################

data<-loadData(parameters)
cat("Total number of genes : ", dim(data$dge$counts)[1], "\n")
cat("Total number of samples : ", dim(data$dge$counts)[2], "\n")
cat("summary of CPM by samples\n")
summary(cpm(data$dge, normalized.lib.sizes=FALSE))
pdf(parameters$output_pdf)
asko_data<-asko3c(data)
cat("Filtering genes with more than ", parameters$threshold_cpm, " CPM in ",parameters$replicate_cpm,"samples\n")
asko_filt<-GEfilt(data$dge, parameters)
cat("Total number of filtered genes : ", dim(asko_filt$counts)[1], "\n")
asko_norm<-GEnorm(asko_filt,parameters)
GEcorr(asko_norm,parameters)
cat("Statistical analysis\n")
DEanalysis(asko_norm,data, asko_data,parameters)
dev.off()


