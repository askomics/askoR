##############################################
##                Parameters                ##  
##############################################
# option_list = list(
#   make_option(c("-o", "--out"), type="character", default="out.pdf", 
#               help="output file name [default= %default]", metavar="character"),
#   make_option(c("-d", "--dir"), type="character", default=".", 
#               help="data directory path [default= %default]", metavar="character"),
#   make_option(c("-o", "--org"), type="character", default="Asko", 
#               help="Organism name [default= %default]", metavar="character"),
#   make_option(c("-f", "--fileofcount"), type="character", default=NULL, 
#               help="file of counts [default= %default]", metavar="character"),
#   make_option(c("-cg", "--col_genes"), type="character", default=NULL, 
#               help="col of ids in count files [default= %default]", metavar="integer"),
#   make_option(c("-cc", "--col_counts"), type="character", default=NULL, 
#               help="col of counts in count files [default= %default]", metavar="integer"),
#   make_option(c("-sep", "--sep"), type="character", default="\t", 
#               help="col separator [default= %default]", metavar="character"),
#   
# ); 
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);


parameters<-list(
  ####Inputs Outputs####
  title_output_pdf = "",
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
  context ="auto",
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
  csv = FALSE,
  heatmap = TRUE,
  numhigh = 50
)