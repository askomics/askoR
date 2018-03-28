##############################################
##                Parameters                ##  
##############################################
parameters<-list(
  title_output_pdf = "DGE analysis on T2 Pb samples without T2MD6TB.pdf",
  dir_path = "C:\\Users\\Sylvin\\Desktop\\Stage_INRA\\DEG_analysis",
  organism = "Pb",
  fileofcount = "SousSetComptagesPlasmo.txt",
  col_genes = NULL,
  col_counts = NULL,
  select_sample="T2",
  rm_sample=list("T2MD6TBR"),
  annotation_file = "annotation_Pb.txt",
  sample_file = "Sample_Herniome_T2.txt",
  contrast_file = "real_contrast_T2.txt",
  threshold_cpm = 0.5,
  replicate_cpm = 3,
  threshold_FDR = 0.05,
  normal_method = "TMM",
  p_adj_method = "BH",
  context ="auto",
  palette ="Set1",
  glm = "lrt",
  ########result_table_option######
  logFC = TRUE,
  FC = TRUE, 
  logCPM = FALSE, 
  FDR = TRUE, 
  LR = FALSE, 
  Sign = TRUE,
  Expression = TRUE, 
  mean_counts = TRUE, 
  csv = FALSE,
  heatmap = TRUE
)