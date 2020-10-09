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
  # Specify desired options in a list
  option_list = list(
    optparse::make_option(c("-o", "--out"), type="character", default="DE_analysis",dest="analysis_name",
                help="output directory name [default= %default]", metavar="character"),
    optparse::make_option(c("-d", "--dir"), type="character", default=".",dest="dir_path",
                help="data directory path [default= %default]", metavar="character"),
    optparse::make_option(c("-O", "--org"), type="character", default="Asko", dest="organism",
                help="Organism name [default= %default]", metavar="character"),
    optparse::make_option(c("-f", "--fileofcount"), type="character", default=NULL, dest="fileofcount",
                help="file of counts [default= %default]", metavar="character"),
    optparse::make_option(c("-G", "--col_genes"), type="integer", default=1, dest="col_genes",
                help="col of genes ids in count files [default= %default]", metavar="integer"),
    optparse::make_option(c("-C", "--col_counts"), type="integer", default=7,dest="col_counts",
                help="col of counts in count files [default= %default (featureCounts) ]", metavar="integer"),
    optparse::make_option(c("-t", "--sep"), type="character", default="\t", dest="sep",
                help="field separator for count files or count matrix [default= %default]", metavar="character"),
    optparse::make_option(c("-a", "--annotation"), type="character", default=NULL, dest="annotation",
                help="annotation file [default= %default]", metavar="character"),
    optparse::make_option(c("-s", "--sample"), type="character", default="Samples.txt", dest="sample_file",
                help="Samples file [default= %default]", metavar="character"),
    optparse::make_option(c("-c", "--contrasts"), type="character", default="Contrasts.txt",dest="contrast_file",
                help="Contrasts file [default= %default]", metavar="character"),
    optparse::make_option(c("-k", "--mk_context"), type="logical", default=FALSE,dest="mk_context",
                help="generate automatically the context names [default= %default]", metavar="logical"),
    optparse::make_option(c("--palette"), type="character", default="Set2", dest="palette",
                help="Color palette (ggplot)[default= %default]", metavar="character"),
    optparse::make_option(c("-R", "--regex"), type="logical", default=FALSE, dest="regex",
                help="use regex when selecting/removing samples [default= %default]", metavar="logical"),
    optparse::make_option(c("-S", "--select"), type="character", default=NULL, dest="select_sample",
                help="selected samples [default= %default]", metavar="character"),
    optparse::make_option(c("-r", "--remove"), type="character", default=NULL, dest="rm_sample",
                help="removed samples [default= %default]", metavar="character"),
    optparse::make_option(c("--th_cpm"), type="double", default=0.5, dest="threshold_cpm",
                help="CPM's threshold [default= %default]", metavar="double"),
    optparse::make_option(c("--rep"), type="integer", default=3, dest="replicate_cpm",
                help="Minimum number of replicates [default= %default]", metavar="integer"),
    optparse::make_option(c("--th_FDR"), type="double", default=0.05, dest="threshold_FDR",
                help="FDR threshold [default= %default]", metavar="double"),
    optparse::make_option(c("-n", "--normalization"), type="character", default="TMM", dest="normal_method",
                help="normalization method (TMM/RLE/upperquartile/none) [default= %default]", metavar="character"),
    optparse::make_option(c("--adj"), type="character", default="fdr", dest="p_adj_method",
                help="p-value adjust method (holm/hochberg/hommel/bonferroni/BH/BY/fdr/none) [default= %default]", metavar="character"),
    optparse::make_option("--glm", type="character", default="qlf", dest="glm",
                help="GLM method (lrt/qlf) [default= %default]", metavar="character"),
    optparse::make_option("--glmDisp", type="logical", default=FALSE, dest="glm_disp",
                help="Estimate Common, Trended and Tagwise Negative Binomial dispersions GLMs (TRUE/FALSE) [default= %default]", metavar="logical"),
    optparse::make_option(c("--lfc"), type="logical", default=TRUE, dest="logFC",
                help="logFC in the summary table [default= %default]", metavar="logical"),
    optparse::make_option(c("--th_lfc"), type="double", default=1, dest="threshold_logFC",
                help="logFC threshold [default= %default]", metavar="double"),
    optparse::make_option("--fc", type="logical", default=TRUE, dest="FC",
                help="FC in the summary table [default= %default]", metavar="logical"),
    optparse::make_option(c("--lcpm"), type="logical", default=FALSE, dest="logCPM",
                help="logCPm in the summary table [default= %default]", metavar="logical"),
    optparse::make_option("--fdr", type="logical", default=TRUE, dest="FDR",
                help="FDR in the summary table [default= %default]", metavar="logical"),
    optparse::make_option("--lr", type="logical", default=FALSE, dest="LR",
                help="LR in the summary table [default= %default]", metavar="logical"),
    optparse::make_option(c("--sign"), type="logical", default=TRUE, dest="Sign",
                help="Significance (1/0/-1) in the summary table [default= %default]", metavar="logical"),
    optparse::make_option(c("--expr"), type="logical", default=TRUE, dest="Expression",
                help="Significance expression in the summary table [default= %default]", metavar="logical"),
    optparse::make_option(c("--mc"), type="logical", default=TRUE, dest="mean_counts",
                help="Mean counts in the summary table [default= %default]", metavar="logical"),
    optparse::make_option(c("--dclust"), type="character", default="euclidean", dest="distcluts",
                help="The distance measure to be used : euclidean, maximum, manhattan, canberra, binary or minkowski [default= %default]", metavar="character"),
    optparse::make_option(c("--hclust"), type="character", default="complete", dest="hclust",
                help="The agglomeration method to be used : ward.D, ward.D2, single, complete, average, mcquitty, median or centroid [default= %default]", metavar="character"),
    optparse::make_option(c("--hm"), type="logical", default=TRUE, dest="heatmap",
                help="generation of the expression heatmap [default= %default]", metavar="logical"),
    optparse::make_option(c("--nh"), type="integer", default="50", dest="numhigh",
                help="number of genes in the heatmap [default= %default]", metavar="integer"),
    optparse::make_option(c("--norm_factor"), type="logical", default=FALSE, dest="norm_factor",
                help="generate file with normalize factor for each condition/sample [default= %default]", metavar="logical"),
    optparse::make_option(c("--norm_counts"), type="logical", default=FALSE, dest="norm_counts",
                help="Generate files with mormalized counts [default= %default]", metavar="logical"),
    optparse::make_option(c("--VD"), type="character", default=NULL, dest="VD",
                help="Plot VennDiagram, precise type of comparison: all, down, up or both [default=%default]", metavar = "character"),
    optparse::make_option(c("--compaVD"), type="character", default=NULL, dest="compaVD",
                help="Contrast comparison list to display in VennDiagram", metavar="character"),
    optparse::make_option(c("--GO"), type="character", default=NULL, dest="GO",
                help="GO enrichment analysis for gene expressed 'up', 'down' or 'both', NULL for no GO enrichment.", metavar="character"),
    optparse::make_option(c("--ID2GO"), type="character", default=NULL, dest="geneID2GO_file",
                help="GO annotation file [default= %default]", metavar="character"),
    optparse::make_option(c("--GO_threshold"), type="numeric", default="0.05", dest="GO_threshold",
                help="the significant threshold used to filter p-values", metavar="double"),
    optparse::make_option(c("--GO_min_num_genes"), type="integer", default="10", dest="GO_min_num_genes",
                help="the minimum number of genes for each GO terms", metavar="integer"),
    optparse::make_option(c("--GO_max_top_terms"), type="integer", default="10", dest="GO_max_top_terms",
                help="the maximum number of GO terms plot", metavar="integer"),
    optparse::make_option(c("--GO_algo"), type="character", default="weight01", dest="GO_algo",
                help="algorithms which are accessible via the runTest function: shown by the whichAlgorithms() function, [default=%default]", metavar="character"),
    optparse::make_option(c("--GO_stats"), type="character", default="fisher", dest="GO_stats",
                help = "statistical tests which are accessible via the runTest function: shown by the whichTests() function, [default=%default]", metavar = "character"),
    optparse::make_option(c("--Ratio_threshold"), type="numeric", default="0", dest="Ratio_threshold",
                help="the minimum ratio value to display GO in graph", metavar="double"),
    optparse::make_option(c("--plotMD"),type="logical", default=FALSE, dest="plotMD", metavar="logical",
                help="Mean-Difference Plot of Expression Data (aka MA plot) [default= %default]"),
    optparse::make_option(c("--plotVO"),type="logical", default=FALSE, dest="plotVO", metavar="logical",
                help="Volcano plot for a specified coefficient/contrast of a linear model [default= %default]"),
    optparse::make_option(c("--glimMD"),type="logical", default=FALSE, dest="glimMD", metavar="logical",
                help="Glimma - Interactif Mean-Difference Plot of Expression Data (aka MA plot) [default= %default]"),
    optparse::make_option(c("--glimVO"),type="logical", default=FALSE, dest="glimVO", metavar="logical",
                help="Glimma - Interactif Volcano plot for a specified coefficient/contrast of a linear model [default= %default]"),
    optparse::make_option(c("--dens_bottom_mar"), type="integer", default="20", dest="densbotmar", metavar="integer",
                help="Set bottom margin of density plot to help position the legend [default= %default]"),
    optparse::make_option(c("--dens_inset"), type="double", default="0.45", dest="densinset", metavar="double",
                help="Set position the legend in bottom density graphe [default= %default]"),
    optparse::make_option(c("--legend_col"), type="integer", default="6", dest="legendcol", metavar="integer",
                help="Set numbers of column for density plot legends [default= %default]"),
    optparse::make_option(c("--upset_basic"),type="character", default=NULL, dest="upset_basic",
                help="Display UpSetR charts for all contrasts, precise type of comparison: all, down, up, mixed [default=%default].", metavar = "character"),
    optparse::make_option(c("--upset_type"),type="character", default=NULL, dest="upset_type",
                help="Display UpSetR charts for list of contrasts, precise type of comparison: all, down, up, mixed [default=%default].", metavar = "character"),
    optparse::make_option(c("--upset_list"), type="character", default=NULL, dest="upset_list",
                help="Contrast comparison list to display in UpSetR chart. See documentation. [default=%default]", metavar="character"),
    optparse::make_option("--coseq_model", type="character", default="kmeans", dest="coseq_model",
                help="Coseq model (kmeans, Normal) [default= %default]", metavar="character"),
    optparse::make_option("--coseq_transformation", type="character", default="clr", dest="coseq_transformation",
                help="Coseq tranformation (arcsin, logit, logclr, clr, alr, ilr, none) [default= %default]", metavar="character"),
    optparse::make_option("--coseq_ClustersNb", type="double", default=2:12, dest="coseq_ClustersNb",
                help="Coseq : number of clusters desired (2:12 (auto), number from 2 to 12) [default= %default]", metavar="double"),
    optparse::make_option("--coseq_ContrastsThreshold", type="integer", default="1", dest="coseq_ContrastsThreshold",
                help="Coseq : number of contrasts in which DE genes are found for clustering [default= %default]", metavar="integer"),
    optparse::make_option("--coseq_normFactors", type="character", default="none", dest="coseq_normFactors",
                help="Coseq normalization factor (TC , UQ, Med , DESeq, TMM , none) [default= %default]", metavar="character")
  )
  # Get command line options
  opt_parser = optparse::OptionParser(option_list=option_list)
  parameters = optparse::parse_args(opt_parser)

  if(is.null(parameters$rm_sample) == FALSE ) {
    stringr::str_replace_all(parameters$rm_sample, " ", "")
    parameters$rm_sample<-limma::strsplit2(parameters$rm_sample, ",")
  }

  if(is.null(parameters$select_sample) == FALSE ) {
    stringr::str_replace_all(parameters$select_sample, " ", "")
    parameters$select_sample<-limma::strsplit2(parameters$select_sample, ",")
  }

  return(parameters)
}
