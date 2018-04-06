##############################################
##                Parameters                ##  
##############################################
option_list = list(
  make_option(c("-o", "--out"), type="character", default="out.pdf",dest="output_pdf",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default=".",dest="dir_path",
              help="data directory path [default= %default]", metavar="character"),
  make_option(c("-o", "--org"), type="character", default="Asko", dest="organism",
              help="Organism name [default= %default]", metavar="character"),
  make_option(c("-f", "--fileofcount"), type="character", default=NULL, dest="fileofcount",
              help="file of counts [default= %default]", metavar="character"),
  make_option(c("-cg", "--col_genes"), type="integer", default=NULL, dest="col_genes",
              help="col of ids in count files [default= %default]", metavar="integer"),
  make_option(c("-cc", "--col_counts"), type="integer", default=NULL,dest="col_counts",
              help="col of counts in count files [default= %default]", metavar="integer"),
  make_option(c("-sep", "--sep"), type="character", default="\t", dest="sep",
              help="col separator [default= %default]", metavar="character"),
  make_option(c("-ann", "--annotation"), type="character", default="annotation.txt", dest="annotation_file",
              help="annotation file [default= %default]", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default="Samples.txt", dest="sample_file",
              help="Samples file [default= %default]", metavar="character"),
  make_option(c("-c", "--contrasts"), type="character", default="Contrasts.txt",dest="contrast_file",
              help="Contrasts file [default= %default]", metavar="character"),
  make_option(c("-pal", "--palette"), type="character", default="Set2", dest="palette",
              help="Contrasts file [default= %default]", metavar="character"),
  make_option(c("-sel", "--select"), type="character", default=NULL, dest="select_sample",
              help="selected samples [default= %default]", metavar="character"),
  make_option(c("-rm", "--remove"), type="character", default=NULL, dest="remove_sample",
              help="removed samples [default= %default]", metavar="character"),
  make_option(c("-th_cpm", "--threshold_cpm"), type="double", default=0.5, dest="threshold_cpm",
              help="CPM's threshold [default= %default]", metavar="double"),
  make_option(c("-rep", "--replicate_count"), type="integer", default=3, dest="replicate_cpm",
              help="Minimum number of replicates [default= %default]", metavar="integer"),
  make_option(c("-th_FDR", "--threshold_FDR"), type="double", default=0.5, dest="threshold_FDR",
              help="FDR threshold [default= %default]", metavar="double"),
  make_option(c("-n", "--normalization"), type="character", default="TMM", dest="normal_method",
              help="normalization method (TMM/RLE/upperquartile/none) [default= %default]", metavar="character"),
  make_option(c("-adj", "--adjust"), type="character", default="fdr", dest="p_adj_method",
              help="p-value adjust method (holm/hochberg/hommel/bonferroni/BH/BY/fdr/none) [default= %default]", metavar="character"),
  make_option(c("-glm", "--glm"), type="character", default="qlf", dest="glm",
              help=" GLM method (lrt/qlf) [default= %default]", metavar="character"),
  make_option(c("-lfc", "--logFC"), type="logical", default="TRUE", dest="logFC",
              help="logFC in the summary table [default= %default]", metavar="logical"),
  make_option(c("-fc", "--fc"), type="logical", default="TRUE", dest="FC",
              help="FC in the summary table [default= %default]", metavar="logical"),
  make_option(c("-lcpm", "--log_cpm"), type="logical", default="FALSE", dest="logCPM",
              help="logCPm in the summary table [default= %default]", metavar="logical"),
  make_option(c("-fdr", "--fdr"), type="logical", default="TRUE", dest="FDR",
              help="FDR in the summary table [default= %default]", metavar="logical"),
  make_option(c("-lr", "--lr"), type="logical", default="FALSE", dest="LR",
              help="LR in the summary table [default= %default]", metavar="logical"),
  make_option(c("-sign", "--significance"), type="logical", default="TRUE", dest="Sign",
              help="Significance (1/0/-1) in the summary table [default= %default]", metavar="logical"),
  make_option(c("-expr", "--expression"), type="logical", default="TRUE", dest="Expression",
              help="Significance expression in the summary table [default= %default]", metavar="logical"),
  make_option(c("-mc", "--mean_counts"), type="logical", default="TRUE", dest="mean_counts",
              help="Mean counts in the summary table [default= %default]", metavar="logical"),
  make_option(c("-hm", "--heatmap"), type="logical", default="TRUE", dest="heatmap",
              help="generation of the expression heatmap [default= %default]", metavar="logical"),
  make_option(c("-hm", "--heatmap"), type="logical", default="TRUE", dest="heatmap",
              help="generation of the expression heatmap [default= %default]", metavar="logical"),
  make_option(c("-nh", "--numhigh"), type="integer", default="50", dest="numhigh",
              help="number of genes in the heatmap [default= %default]", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
parameters = parse_args(opt_parser);


parameters<-list(
  ####Inputs Outputs####
  output_pdf = "out.pdf",
  dir_path = "",
  organism = "Asko",
  fileofcount = NULL,
  col_genes = NULL,
  col_counts = NULL,
  sep="\t",
  annotation_file = "annotation.txt",
  sample_file = "Samples.txt",
  contrast_file = "Contrasts.txt",
  #### Color palette ####
  palette ="Set2",
  ####Selection####
  select_sample=NULL,
  rm_sample=NULL,
  #### Filtration #####
  threshold_cpm = 0.5,
  replicate_cpm = 3,
  threshold_FDR = 0.05,
  #### Normalisation ####
  normal_method = "TMM",
  #### Test ####
  p_adj_method = "BH",
  mk_context ="auto",
  glm = "lrt",
  ####result_table_option####
  logFC = TRUE,
  FC = TRUE, 
  logCPM = FALSE, 
  FDR = TRUE, 
  LR = FALSE, 
  Sign = TRUE,
  Expression = TRUE, 
  mean_counts = TRUE, 
  heatmap = TRUE,
  numhigh = 50
)
