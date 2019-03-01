# Removes all objects from the current workspace (R memory) 
rm(list=ls())
askor_path<-Sys.getenv("ASKOR_PATH")
source(paste0(askor_path,"/AskoR.R"))

##############################################
##                Parameters                ##  
##############################################
setwd("~/PATH_TO_YOUR_WORKSPACE/")
parameters<-Asko_start()

# Data and input files descriptions
#--------------------------------------------------------------------------
# WARNING: All the input files must be in the same folder 
#          called "input" (case sensitive)!
#--------------------------------------------------------------------------
# parameters$analysis_name = "DE_analysis "           # output directory name (default DE_analysis, do not put space!)
# parameters$dir_path = "/"                           # workspace directory (default ".")
# parameters$organism = "Asko"                        # output files prefix (default Asko, do not put space!)
# parameters$fileofcount = ""                         # matrix of count for all samples/conditions
# parameters$GO_MF = ""                               # GO Molecular Function : 2 columns idgene et GO (one GO per line, idgene can be repeated)
# parameters$GO_BP = ""                               # GO Biological Process : 2 columns idgene et GO (one GO per line, idgene can be repeated)
# parameters$GO_CC = ""                               # GO Cellular Component : 2 columns idgene et GO (one GO per line, idgene can be repeated)
# parameters$annotation = ""                          # file containing the functional annotations of each gene 
parameters$sample_file = "Samples.csv"                # file describing the samples
parameters$contrast_file = "Contrasts.csv"            # matrix of different contrasts desired
parameters$sep = "\t"                                 # field separator for files counts or matrix counts
parameters$col_genes = 1                              # column with the gene names (default 1) for count files
parameters$col_counts = 7                             # ol of counts in count files (default 7)
# parameters$select_sample = c("")                    # selected samples
# parameters$rm_sample = c("")                        # removed samples
# parameters$regex = T                                # use regex when selecting/removing samples (default FALSE)
# parameters$mk_context = T                           # generate automatically the context names  (default FALSE)
# parameters$norm_mean = T                            # generate file with mormalized mean for each condition/sample, in Askomics format (default FALSE)


# Options for data processing and their analyzes
#--------------------------------------------------------------------------
parameters$threshold_cpm = 0.5                        # CPM's threshold (default 0.5)
parameters$replicate_cpm = 3                          # Minimum number of replicates (default 3) 
parameters$threshold_FDR = 0.05                       # FDR threshold (default 0.05)
parameters$threshold_logFC = 1                        # logFC threshold (default 1)
parameters$normal_method = "TMM"                      # normalization method (TMM/RLE/upperquartile/none) (default TMN)
parameters$p_adj_method = "BH"                        # p-value adjust method (holm/hochberg/hommel/bonferroni/BH/BY/fdr/none) (default fdr)
parameters$glm = "lrt"                                # GLM method (lrt/qlf) (default qlf)
parameters$logFC = T                                  # logFC in the summary table (default TRUE)
parameters$FC = T                                     # FC in the summary table (default TRUE)
parameters$logCPM = F                                 # logCPm in the summary table (default FALSE)
parameters$FDR = T                                    # FDR in the summary table (default TRUE)
parameters$LR = F                                     # LR in the summary table (default FALSE)
parameters$Sign = T                                   # Significance (1/0/-1) in the summary table (default TRUE)
parameters$Expression = T                             # Significance expression in the summary table (default TRUE)
parameters$mean_counts = T                            # Mean counts in the summary table (default TRUE)

# for hierarchical clustering
#-----------------------------------
# parameters$distcluts="euclidiean"                   # The distance measure to be used : "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski" (default euclidean)
# parameters$hclust="complete"                        # The agglomeration method to be used : "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid" (default complete)

# for heatmap plot 
#-----------------------------------
# parameters$palette = "Set3"                         # Color palette (ggplot) (default Set2)
# parameters$heatmap = T                              # generation of the expression heatmap (default TRUE)
# parameters$numhigh = 50                             # number of genes in the heatmap (default 50)

# Visualization of results from differential expression analyses
#-----------------------------------
# parameters$plotMD = T                               # Mean-Difference Plot of Expression Data (aka MA plot) (default FALSE)
# parameters$plotVO = T                               # Volcano plot for a specified coefficient/contrast of a linear model (default FALSE)
# parameters$glimMD = T                               # Glimma - Interactif Mean-Difference Plot of Expression Data (aka MA plot) (default FALSE)
# parameters$glimVO = T                               # Glimma - Interactif Volcano plot for a specified coefficient/contrast of a linear model (default FALSE)

# Parameters for Venn diagrams
#----------------------------------------------------------------------
# Plot VennDiagram, precise type of comparison: all, down, up, both (default NULL)
# "all"  : Create VennDiagrams for all differentially expressed genes
# "up"   : Create VennDiagrams for gene expressed UP
# "down" : Create VennDiagrams for gene expressed DOWN
# "both" : Create VennDiagrams for gene expressed UP and DOWN (in the same graph) 
#  NULL  : Not display Venn diagrams
# parameters$VD = NULL 

# Contrast comparison list to display in VennDiagram
# WARNING: if VD is not null, compaVD must NOT BE EMPTY!

# For VD = "all", "up" or "down" -- accepts up to 5 items separated by a dash and accepts lists too
# parameters$compaVD=c("Cond1vsCond2-Cond1vsCond3-Cond6vsCond3")                  # this create 1 venn diagram
# parameters$compaVD=c("E1Cond1vsE1Cond2-E1Cond1vsE1Cond3-E1Cond6vsE1Cond3",      # this create 3 venn diagrams
#                      "E2Cond1vsE2Cond2-E2Cond1vsE2Cond3-E2Cond6vsE2Cond3",
#                      "E3Cond1vsE3Cond2-E3Cond1vsE3Cond3-E3Cond6vsE3Cond3")

# For VD = "both" -- accept only 2 items and accepts lists too
# parameters$compaVD=c("Cond1vsCond2-Cond1vsCond3")           # this create one venn diagram
# parameters$compaVD=c("Cond1vsCond2-Cond1vsCond3",           # this create 3 venn diagrams
#                      "Cond1vsCond2-Cond1vsCond6",
#                      "Cond1vsCond6-Cond1vsCond3") 

# Parameters for GO enrichment
#----------------------------------------------------------------------
parameters$GO = "both"                    # gene set chosen for analysis 'up', 'down', 'both', or NULL (default NULL)  
parameters$GO_filt_meth = "p.adjust"      # Use 'pval' to filter on nominal p-value or 'p.adjust' to filter on adjusted p-value (default p.adjust)
parameters$GO_padj_meth = "BH"            # correction method used to adjust p-values; available option : 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none' (default BH)
parameters$GO_threshold = 0.05            # the significant threshold used to filter p-values (default 0.05)
parameters$GO_min_num_terms = 1           # the minimum number of GO terms required to plot a cluster label (default 10)

########################################
##  Loading the data from the samples ##
########################################
##### load data #####
data<-loadData(parameters)

cat("\n\nChecking the data content:\n")
data$samples
data$contrast
data$design
head(data$dge$counts)

cat("Total number of genes : ", dim(data$dge$counts)[1], "\n")
cat("Total number of samples : ", dim(data$dge$counts)[2], "\n\n")
cat("Summary of CPM by samples\n")
summary(cpm(data$dge))
cat("\n")

##### asko files #####
asko_data<-asko3c(data)
cat("\nChecking Asko Data : condition, contrast, context.\n")
asko_data$condition ; cat("\n")
asko_data$contrast  ; cat("\n")
asko_data$context   ; cat("\n")

##### filtering #####
cat("\nFiltering genes with more than ", parameters$threshold_cpm, " CPM in ",parameters$replicate_cpm,"samples\n")
asko_filt<-GEfilt(data, parameters)
cat("Total number of filtered genes : ", dim(asko_filt$counts)[1], "\n\n")

##### normalization #####
asko_norm<-GEnorm(asko_filt, asko_data, parameters)

##### correlation #####
GEcorr(asko_norm,parameters)

##### DGE analysis #####
cat("\n\nDifferential expressions analysis\n")
sum_table<-DEanalysis(asko_norm, data, asko_data, parameters)

##### Venn diagram #####
VD(sum_table, parameters, asko_data)

##### Enrichment Analysis #####
if(is.null(data$GO_MF)==FALSE){ matMF<-runGoStag(sum_table, asko_data, data$GO_MF, "MF") }
if(is.null(data$GO_BP)==FALSE){ matBP<-runGoStag(sum_table, asko_data, data$GO_BP, "BP") }
if(is.null(data$GO_CC)==FALSE){ matCC<-runGoStag(sum_table, asko_data, data$GO_CC, "CC") }





