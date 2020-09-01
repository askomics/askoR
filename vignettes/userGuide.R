## ---- echo=FALSE--------------------------------------------------------------
knitr::opts_chunk$set(echo=TRUE, warning=FALSE, message=FALSE, fig.align="left", fig.show="hold", fig.keep='all')

## ----counts_files, eval=FALSE-------------------------------------------------
#  parameters$col_genes=1
#  parameters$col_counts=7
#  parameters$sep="\t"

## ----run_askor1, eval=FALSE---------------------------------------------------
#  # Path to askoR file
#  library(askoR)
#  
#  # Sets defaults parameters
#  parameters<-Asko_start()

## ----run_askor2, eval=FALSE---------------------------------------------------
#  # output directory name (default DE_analysis)
#  parameters$analysis_name="DEG_test"
#  
#  # input files:
#  # matrix of different contrasts desired
#  parameters$contrast_file = "Contrasts.txt"
#  # file containing the functional annotations for each gene
#  parameters$annotation = "Genes_annotations.txt"
#  # GO annotation files
#  parameters$geneID2GO_file = "GO_annotations.txt"

## ----run_askor3a, eval=FALSE--------------------------------------------------
#  # matrix of count for all samples/conditions
#  parameters$fileofcount = "CountsMatrix.txt"
#  # file describing all samples
#  parameters$sample_file = "Samples_CountsMatrix.txt"

## ----run_askor3b, eval=FALSE--------------------------------------------------
#  # file describing all samples
#  parameters$sample_file = "Samples_CountsFiles.txt"
#  # column with the gene names (default 1)
#  parameters$col_genes = 1
#  # column with the counts values (default 7)
#  parameters$col_counts = 7
#  # field separator (default "\t")
#  parameters$sep = "\t"

## ----run_askor3c, eval=FALSE--------------------------------------------------
#  # delete sample AC3R2
#  parameters$rm_sample = c("AC3R2","BC3R3")

## ----run_askor4, eval=FALSE---------------------------------------------------
#  data<-loadData(parameters)

## ----exec1, echo=FALSE, message=FALSE, warning=FALSE--------------------------
rm(list=ls())
library(askoR)
parameters<-Asko_start()
parameters$dir_path="../inst/extdata/"
parameters$analysis_name="DEG_test"
parameters$fileofcount = "CountsMatrix.txt"  
parameters$sample_file = "Samples_CountsMatrix.txt"  
parameters$contrast_file = "Contrasts.txt"    
parameters$annotation = "Genes_annotations.txt"  
parameters$geneID2GO_file = "GO_annotations.txt"
parameters$rm_sample = c("AC3R2","BC3R3")
data<-loadData(parameters)

## ----run_askor4b--------------------------------------------------------------
# Displays all samples recorded
data$samples  

# Displays all contrast recorded
data$contrast 

# Displays design experiment
data$design   

# Displays the first 5 lines and 8 columns of counts table.
data$dge$counts[1:5,1:8] 

# Total number of genes:
dim(data$dge$counts)[1]

# Total number of samples:
dim(data$dge$counts)[2]

## ----run_askor5---------------------------------------------------------------
asko_data<-asko3c(data, parameters)

## ----run_askor6a--------------------------------------------------------------
# CPM's threshold 
parameters$threshold_cpm = 0.5  
# minimum of sample which are upper to cpm threshold 
parameters$replicate_cpm = 3 # we have 3 replicates

## ----run_askor6b--------------------------------------------------------------
# run filtering 
asko_filt<-GEfilt(data, parameters)

# Total number of filtered genes: 
dim(asko_filt$counts)[1]

## ----filters1a, echo=FALSE, out.width='43%', fig.cap="Filtering Data", out.extra='style="display:inline-block;margin:0;padding:0;"'----
knitr::include_graphics("../inst/extdata/DEG_test/images/DEG_test_boxplot_logcpm_before_filtering.png")
knitr::include_graphics("../inst/extdata/DEG_test/images/DEG_test_boxplot_logcpm_after_filtering.png")
knitr::include_graphics("../inst/extdata/DEG_test/images/DEG_test_barplot_logcpm_before_filtering.png")
knitr::include_graphics("../inst/extdata/DEG_test/images/DEG_test_barplot_logcpm_after_filtering.png")

## ----filters1b, echo=FALSE, out.width='43%', fig.cap="Density Graphes", out.extra='style="display:inline-block;margin:0;padding:0;"'----
knitr::include_graphics("../inst/extdata/DEG_test/images/DEG_test_raw_data_1.png")
knitr::include_graphics("../inst/extdata/DEG_test/images/DEG_test_filtered_data_1.png")

## ----param1-------------------------------------------------------------------
# Set position the legend in bottom density graphe
parameters$densinset = 0.20
# Set numbers of column for legends
parameters$legendcol = 8
# run filtering
asko_filt<-GEfilt(data, parameters)

## ----filters2, echo=FALSE, out.width='43%', fig.cap="Density Graphes Corrected", out.extra='style="display:inline-block;margin:0 0 -50px 0;padding:0;"'----
knitr::include_graphics("../inst/extdata/DEG_test/images/DEG_test_raw_data.png")
knitr::include_graphics("../inst/extdata/DEG_test/images/DEG_test_filtered_data.png")

## ----run_askor7---------------------------------------------------------------
# run normalization
asko_norm<-GEnorm(asko_filt, asko_data, data, parameters)

## ----norm, fig.cap="Normalization graphs (barplot)", echo=FALSE, out.width='43%', fig.align='center', out.extra='style="display:inline-block;margin:0;padding:0;"'----
knitr::include_graphics("../inst/extdata/DEG_test/images/DEG_test_boxplot_logcpm_after_norm.png")

## ----norm2, fig.cap="Normalization graphs (heatmap)", echo=FALSE, out.width='43%', out.extra='style="display:inline-block;margin:0;padding:0;"'----
knitr::include_graphics("../inst/extdata/DEG_test/images/DEG_test_heatmap_CPMcounts_per_sample.png")
knitr::include_graphics("../inst/extdata/DEG_test/images/DEG_test_heatmap_meanCounts_per_condi.png")

## ----run_askor8---------------------------------------------------------------
GEcorr(asko_norm,parameters)

## ----corr, fig.cap="MDS and PCA plots - axis 1 and 2", echo=FALSE, out.width='43%', out.extra='style="display:inline-block;margin:0;padding:0;"'----
knitr::include_graphics("../inst/extdata/DEG_test/images/DEG_test_MDS_corr_axe1_2.png")
knitr::include_graphics("../inst/extdata/DEG_test/images/DEG_test_PCA_axe1_2.png")

## ----run_askor9a--------------------------------------------------------------
# FDR threshold
parameters$threshold_FDR = 0.05
# logFC threshold
parameters$threshold_logFC = 0
# normalization method
parameters$normal_method = "TMM"
# p-value adjust method
parameters$p_adj_method = "BH"
# GLM method
parameters$glm = "lrt"

## ----run_askor9b--------------------------------------------------------------
# Mean-Difference Plot of Expression Data (aka MA plot)
parameters$plotMD = TRUE
# Volcano plot for a specified coefficient/contrast of a linear model
parameters$plotVO = TRUE

## ----run_askor9c--------------------------------------------------------------
# run differential expression analysis
resDEG<-DEanalysis(asko_norm, data, asko_data, parameters)

## ----dge, fig.cap="DEanalysis plots", echo=FALSE, out.width='43%', out.extra='style="display:inline-block;margin:0;padding:0;"'----
knitr::include_graphics("../inst/extdata/DEG_test/images/AC3vsBC3_MeanDifference_of_ExpressionData.png")
knitr::include_graphics("../inst/extdata/DEG_test/images/AC3vsBC3_VolcanoPlot.png")

## ----dge2, fig.cap="DEanalysis plots", echo=FALSE, out.width='50%', fig.align='center', out.extra='style="display:inline-block;margin:0;padding:0;"'----
knitr::include_graphics("../inst/extdata/DEG_test/images/AC3vsBC3_topDGE_heatmap.png")

## ----exemple1-----------------------------------------------------------------
# this create 1 venn diagram
parameters$compaVD=c("Ctrast1-Ctrast2-Ctrast3") 
# this create 3 venn diagrams
parameters$compaVD=c("Ctrast1-Ctrast2-Ctrast3", 
                     "Ctrast4-Ctrast5-Ctrast6",
                     "Ctrast7-Ctrast8-Ctrast9")

## ----exemple2-----------------------------------------------------------------
# this create 1 venn diagram
parameters$compaVD=c("Ctrast1-Ctrast2") 
# this create 3 venn diagrams
parameters$compaVD=c("Ctrast1-Ctrast2", 
                     "Ctrast1-Ctrast3",
                     "Ctrast2-Ctrast3")

## ----venn1, eval=FALSE--------------------------------------------------------
#  parameters$compaVD = c("AC1vsAC2-AC1vsAC3-AC2vsAC3",
#                         "BC1vsBC2-BC1vsBC3-BC2vsBC3",
#                         "AC1vsBC1-AC2vsBC2-AC3vsBC3")
#  
#  # graph type "all"
#  parameters$VD = "all"
#  VD(resDEG, parameters, asko_data)
#  
#  # graph type "up"
#  parameters$VD = "up"
#  VD(resDEG, parameters, asko_data)
#  
#  # graph type "down"
#  parameters$VD = "down"
#  VD(resDEG, parameters, asko_data)

## ----venn2, eval=FALSE--------------------------------------------------------
#  # graph type "both"
#  parameters$compaVD = c("AC1vsBC1-AC2vsBC2",
#                         "AC1vsBC1-AC3vsBC3",
#                         "AC2vsBC2-AC3vsBC3")
#  parameters$VD = "both"
#  VD(resDEG, parameters, asko_data)

## ----venndiagram, fig.cap="Venn Diagrams", echo=FALSE, out.width='43%', out.extra='style="display:inline-block;margin:0;padding:0;"'----
knitr::include_graphics("../inst/extdata/DEG_test/vennDiagram/AC1vsBC1-AC2vsBC2-AC3vsBC3_all.png")
knitr::include_graphics("../inst/extdata/DEG_test/vennDiagram/AC1vsBC1-AC2vsBC2-AC3vsBC3_up.png")
knitr::include_graphics("../inst/extdata/DEG_test/vennDiagram/AC1vsBC1-AC2vsBC2-AC3vsBC3_down.png")
knitr::include_graphics("../inst/extdata/DEG_test/vennDiagram/AC1vsBC1-AC2vsBC2_mixed.png")

## ----exemple3-----------------------------------------------------------------
# Precise type of comparison: all, down, up, mixed.
parameters$upset_type = "all"

# Give a list of contrast, for example:
# this create 1 graphs
parameters$upset_list = c("Ctrast1-Ctrast2-Ctrast3")   
# this create 3 graphs
parameters$upset_list = c("Ctrast1-Ctrast2-Ctrast3",  
                          "Ctrast4-Ctrast5-Ctrast6",
                          "Ctrast1-Ctrast2-Ctrast3-Ctrast4-Ctrast5")

## ----upset, eval=FALSE--------------------------------------------------------
#  parameters$upset_list = c("AC1vsAC2-AC1vsAC3-AC2vsAC3",
#                            "BC1vsBC2-BC1vsBC3-BC2vsBC3",
#                            "AC1vsBC1-AC2vsBC2-AC3vsBC3")
#  
#  # graphs type "all"
#  parameters$upset_basic = "all" # all contrast
#  parameters$upset_type = "all"  # list of contrast
#  UpSetGraph(resDEG, data, parameters)
#  
#  # graphs type "mixed"
#  parameters$upset_basic = "mixed" # all contrast
#  parameters$upset_type = "mixed"  # list of contrast
#  UpSetGraph(resDEG, data, parameters)
#  
#  # graphs type "up"
#  parameters$upset_basic = "up" # all contrast
#  parameters$upset_type = "up"  # list of contrast
#  UpSetGraph(resDEG, data, parameters)
#  
#  # graphs type "down"
#  parameters$upset_basic = "down" # all contrast
#  parameters$upset_type = "down"  # list of contrast
#  UpSetGraph(resDEG, data, parameters)

## ----upsetR, fig.cap="UpsetR Graphs", echo=FALSE, out.width='43%',  out.extra='style="display:inline-block;margin:0;padding:0;"'----
knitr::include_graphics("../inst/extdata/DEG_test/UpSetR_graphs/Subset_upset/DEG_test_UpSetR_AC1vsBC1-AC2vsBC2-AC3vsBC3_allDEG.png")
knitr::include_graphics("../inst/extdata/DEG_test/UpSetR_graphs/Subset_upset/DEG_test_UpSetR_AC1vsBC1-AC2vsBC2-AC3vsBC3_upDEG.png")
knitr::include_graphics("../inst/extdata/DEG_test/UpSetR_graphs/Subset_upset/DEG_test_UpSetR_AC1vsBC1-AC2vsBC2-AC3vsBC3_downDEG.png")
knitr::include_graphics("../inst/extdata/DEG_test/UpSetR_graphs/Subset_upset/DEG_test_UpSetR_AC1vsBC1-AC2vsBC2-AC3vsBC3_mixedDEG.png")

## ----GO, eval=FALSE-----------------------------------------------------------
#  # Parameters
#  parameters$GO_threshold = 0.05
#  parameters$GO_max_top_terms = 10
#  parameters$GO_min_num_genes = 10
#  parameters$GO = "both"
#  parameters$GO_algo = "weight01"
#  parameters$GO_stats = "fisher"
#  parameters$Ratio_threshold = 1
#  
#  # run analysis
#  GOenrichment(resDEG, data, parameters)

## ----GOgraph, fig.cap="GO Enrichment graphs", echo=FALSE, out.width='43%',  out.extra='style="display:inline-block;margin:0;padding:0;"'----
knitr::include_graphics("../inst/extdata/DEG_test/GO_images/AC1vsBC1_Pvalue_BUBBLESgraph.png")
knitr::include_graphics("../inst/extdata/DEG_test/GO_images/AC1vsBC1_Ratio_BUBBLESgraph.png")

