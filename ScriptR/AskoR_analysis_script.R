# Removes all objects from the current workspace (R memory)
rm(list=ls())
askor_path<-"<path AskoR folder" # !!!!!!!!!!!!!!!!!!!! change path !!!!!!!!!!!!!!!!!!!!
source(paste0(askor_path,"/AskoR.R"))

##############################################
##                Parameters                ##  
##############################################
setwd("path to workspace") # !!!!!!!!!!!!!!!!!!!! change path !!!!!!!!!!!!!!!!!!!!
parameters<-Asko_start()

# Data and input files descriptions
#--------------------------------------------------------------------------
# WARNING: All the input files must be in the same folder 
#          called "input" (case sensitive)!
#--------------------------------------------------------------------------
parameters$analysis_name = "TEST"            # output directory name (default DE_analysis, do not put space!)
# parameters$dir_path = "/"                  # workspace directory (default ".")
# parameters$organism = ""                   # output files prefix (default Asko, do not put space!)
parameters$fileofcount = "Counts.txt"        # matrix of count for all samples/conditions
parameters$annotation = "annot_genes.csv"    # file containing the functional annotations of each gene 
parameters$sample_file = "Samples.tsv"       # file describing the samples
parameters$contrast_file = "Contrasts.tsv"   # matrix of different contrasts desired
parameters$sep = "\t"                        # field separator for files counts or matrix counts
# parameters$col_genes = 1                   # column with the gene names (default 1) in the counting matrix or count files
# parameters$col_counts = 7                  # ol of counts in count files (default 7)
# parameters$select_sample = c("")           # selected samples
# parameters$rm_sample = c("")               # removed samples
# parameters$regex = T                       # use regex when selecting/removing samples (default FALSE)
# parameters$mk_context = T                  # generate automatically the context names  (default FALSE)
# parameters$norm_mean = T                   # generate file with mormalized mean for each condition/sample, in Askomics format (default FALSE)

# Options for data processing and their analyzes 
#--------------------------------------------------------------------------
parameters$threshold_cpm = 0.5               # CPM's threshold (default 0.5)
parameters$replicate_cpm = 3                 # Minimum number of replicates (default 3) 
parameters$threshold_FDR = 0.05              # FDR threshold (default 0.05)
parameters$threshold_logFC = 0               # logFC threshold (default 1)
parameters$normal_method = "TMM"             # normalization method (TMM/RLE/upperquartile/none) (default TMN)
parameters$p_adj_method = "BH"               # p-value adjust method (holm/hochberg/hommel/bonferroni/BH/BY/fdr/none) (default fdr)
parameters$glm = "lrt"                       # GLM method (lrt/qlf) (default qlf)
parameters$glm_disp = F                      # Estimate Common, Trended and Tagwise Dispersion for Negative Binomial GLMs (default FALSE)
parameters$logFC = T                         # logFC in the summary table (default TRUE)
parameters$FC = T                            # FC in the summary table (default TRUE)
parameters$logCPM = F                        # logCPm in the summary table (default FALSE)
parameters$FDR = T                           # FDR in the summary table (default TRUE)
parameters$LR = F                            # LR in the summary table (default FALSE)
parameters$Sign = T                          # Significance (1/0/-1) in the summary table (default TRUE)
parameters$Expression = T                    # Significance expression in the summary table (default TRUE)
parameters$mean_counts = T                   # Mean counts in the summary table (default TRUE)

# for legend of density plot 
#-----------------------------------
# parameters$densbotmar=20                   # Set bottom margin of density plot to help position the legend (default 20)
# parameters$densinset=0.45                  # Set position the legend in bottom density graphe (default 0.45)

# for hierarchical clustering
#-----------------------------------
# parameters$distcluts="euclidiean"          # The distance measure to be used : "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski" (default euclidean)
# parameters$hclust="complete"               # The agglomeration method to be used : "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid" (default complete)

# for heatmap plot 
#-----------------------------------
# parameters$palette = "Set3"                # Color palette (ggplot) (default Set2)
# parameters$heatmap = T                     # generation of the expression heatmap (default TRUE)
# parameters$numhigh = 50                    # number of genes in the heatmap (default 50)

# Visualization of results from differential expression analyses
#-----------------------------------
# parameters$plotMD = T                      # Mean-Difference Plot of Expression Data (aka MA plot) (default FALSE)
# parameters$plotVO = T                      # Volcano plot for a specified coefficient/contrast of a linear model (default FALSE)
# parameters$glimMD = T                      # Glimma - Interactif Mean-Difference Plot of Expression Data (aka MA plot) (default FALSE)
# parameters$glimVO = T                      # Glimma - Interactif Volcano plot for a specified coefficient/contrast of a linear model (default FALSE)

# Parameters for GO enrichment
#----------------------------------------------------------------------
# parameters$GO = "both"                     # gene set chosen for analysis 'up', 'down', 'both', or NULL (default NULL)  
# parameters$GO_filt_meth = "p.adjust"       # Use 'pval' to filter on nominal p-value or 'p.adjust' to filter on adjusted p-value (default p.adjust)
# parameters$GO_padj_meth = "BH"             # correction method used to adjust p-values; available option : 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none' (default BH)
# parameters$GO_threshold = 0.05             # the significant threshold used to filter p-values (default 0.05)
# parameters$GO_min_num_terms = 1            # the minimum number of GO terms required to plot a cluster label (default 10)


########################################
##  Loading the data from the samples ##
########################################
##### load data #####
data<-loadData(parameters)

cat("\n\nChecking Data content:\n")
data$samples
data$contrast
data$design
head(data$dge$counts)

cat("Total number of genes : ", dim(data$dge$counts)[1], "\n")
cat("Total number of samples : ", dim(data$dge$counts)[2], "\n\n")
cat("summary of CPM by samples\n")
summary(cpm(data$dge))
cat("\n")

##### asko files #####
asko_data<-asko3c(data, parameters)
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
resDEG<-DEanalysis(asko_norm, data, asko_data, parameters)

###### UpsetR Graphs #####
#---------------------------------------------------------------------------------
# UpsetR graphs parameters 
#---------------------------------------------------------------------------------
# Display upset charts for all contrasts
#---------------------------------------------------------------------------------
# Need precise type of comparison: all, down, up, mixed.
# "all"   : Create chart for all differentially expressed genes 
# "up"    : Create chart for gene expressed UP
# "down"  : Create chart for gene expressed DOWN
# "mixed" : Create chart for gene expressed UP and DOWN (in the same graph) 
# parameters$upset_basic = "all"
#---------------------------------------------------------------------------------
# Contrast comparison list to display
#---------------------------------------------------------------------------------
# Precise type of comparison: all, down, up, mixed.
# parameters$upset_type = "all" 
# Give a list of contrast, for example:
# parameters$upset_list = c("Ctrast1-Ctrast2-Ctrast3")   # this create 1 graphs
# parameters$upset_list = c("Ctrast1-Ctrast2-Ctrast3",   # this create 3 graphs
#                           "Ctrast4-Ctrast5-Ctrast6",
#                           "Ctrast1-Ctrast2-Ctrast3-Ctrast4-Ctrast5-Ctrast6")
#---------------------------------------------------------------------------------

# My list
parameters$upset_list = c("Ctrast1-Ctrast2-Ctrast3",
                          "Ctrast4-Ctrast5-Ctrast6",
                          "Ctrast1-Ctrast2-Ctrast3-Ctrast4-Ctrast5-Ctrast6")

# graphs type "all"
parameters$upset_basic = "all"
parameters$upset_type = "all"
UpSetGraph(resDEG, data, parameters)
 
# graphs type "mixed"
parameters$upset_basic = "mixed"
parameters$upset_type = "mixed"
UpSetGraph(resDEG, data, parameters)
 
# graphs type "up"
parameters$upset_basic = "up"
parameters$upset_type = "up"
UpSetGraph(resDEG, data, parameters)
 
# graphs type "down"
parameters$upset_basic = "down"
parameters$upset_type = "down"
UpSetGraph(resDEG, data, parameters)

##### Venn diagram #####
#---------------------------------------------------------------------------------
# Parameters for Venn diagrams
#---------------------------------------------------------------------------------
# Plot VennDiagram, precise type of comparison: all, down, up, both (default NULL)
# "all"  : Create VennDiagrams for all differentially expressed genes
# "up"   : Create VennDiagrams for gene expressed UP
# "down" : Create VennDiagrams for gene expressed DOWN
# "both" : Create VennDiagrams for gene expressed UP and DOWN (in the same graph) 
#  NULL  : Not display Venn diagrams
# parameters$VD = "all" 
#---------------------------------------------------------------------------------
# Contrast comparison list to display in VennDiagram
# WARNING: if VD is not null, compaVD must NOT BE EMPTY!
# For VD = "all", "up" or "down" -- accepts up to 5 items separated by a dash and accepts lists too
# parameters$compaVD=c("Ctrast1-Ctrast2-Ctrast3")   # this create 1 venn diagram
# parameters$compaVD=c("Ctrast1-Ctrast2-Ctrast3",   # this create 3 venn diagrams
#                      "Ctrast4-Ctrast5-Ctrast6",
#                      "Ctrast7-Ctrast8-Ctrast9")
#
# For VD = "both" -- accept only 2 items and accepts lists too
# parameters$compaVD=c("Ctrast1-Ctrast2")      # this create 1 venn diagram
# parameters$compaVD=c("Ctrast1-Ctrast2",      # this create 3 venn diagrams
#                      "Ctrast1-Ctrast3",
#                      "Ctrast2-Ctrast3") 
#---------------------------------------------------------------------------------

# My list
parameters$compaVD=c("Ctrast1-Ctrast2-Ctrast3",  
                     "Ctrast4-Ctrast5-Ctrast6",
                     "Ctrast7-Ctrast8-Ctrast9")

# graph type "all"
parameters$VD = "all"
VD(resDEG, parameters, asko_data)

# graph type "up"
parameters$VD = "up"
VD(resDEG, parameters, asko_data)

# graph type "down"
parameters$VD = "down"
VD(resDEG, parameters, asko_data)

# graph type "both"
parameters$compaVD=c("Ctrast1-Ctrast2", 
                     "Ctrast1-Ctrast3",
                     "Ctrast2-Ctrast3") 
parameters$VD = "both"
VD(resDEG, parameters, asko_data)

##### Enrichment Analysis #####
# Parameters
#----------------------------------------------------------------------
parameters$geneID2GO_file = "Plasmo_GO.tsv"  # GO annotation files
parameters$GO_threshold = 0.05               # the significant threshold used to filter p-values (default 0.05)
parameters$GO_max_top_terms = 5              # the maximum number of GO terms plot (default 10)
parameters$GO_min_num_genes = 5              # the minimum number of genes for each GO terms (default 10)
parameters$GO = "both"                       # gene set chosen for analysis 'up', 'down', 'both', or NULL (default NULL)
parameters$GO_algo = "classic"               # algorithms which are accessible via the runTest function
parameters$GO_stats = "fisher"               # statistical tests which are accessible via the runTest function

GOenrichment(resDEG, parameters)

# Parameters for Co-expression analysis
#----------------------------------------------------------------------
# parameters$coseq_model = "kmeans"        # (default kmeans)
# parameters$coseq_transformation = "clr"  # (default clr)
# parameters$coseq_ClustersNb = 2:12       # (default : automatic selection between 2 to 12, you can fix the number of clusters to be build only with the "Normal" model)
# parameters$coseq_ContrastsThreshold = 1  # (default 1)
ClustAndGO(asko_norm,resDEG,parameters)

