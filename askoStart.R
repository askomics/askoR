##############################################
##                Parameters                ##  
##############################################
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