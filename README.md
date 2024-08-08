## AskoR - RNAseq data analysis tool

Authors: Susete Alves Carvalho, Kévin Gazengel (co-first), Stéphanie Robin, Sylvin Masanelli, Stéphanie Daval, Fabrice Legeai <br/>

**AskoR** is an R pipeline for gene expression analysis from RNAseq data. It allows, through the articulation of well-known R packages (EdgeR, UpSetR, topGO, Coseq), to perform a complete and automatic transcriptomic analysis: data filtering, normalisation of filtered data, descriptive analysis of data quality, differential gene expression analysis, enrichment in terms of "Gene Ontology", clustering genes by co-expression profile, and enhanced with multiple graphic representations.

### Sample data
You'll find a test set in the **inst/extdata/input** folder. It'll be used for the vignettes documentation.<br/>

  - Count matrix file: CountsMatrix.txt **OR** Counts files per samples: files in "counts" directory _(counts.tgz)_
  - Samples file: Samples_CountsMatrix.txt **OR** Samples_CountsFiles.txt 
  - Contrasts file: Contrasts.txt
  - Genes annotations file: Genes_annotations.txt (optional)
  - GO annotations file: GO_annotations.txt (optional)

**IMPORTANT :** All input files must be in a folder named **input** _(case sensitive)_.

### Vignettes

* [User Guide](https://github.com/askomics/askoR/wiki/Pipeline-askoR:-User-Guide)
* [Parameters Table](https://github.com/askomics/askoR/wiki/Pipeline-askoR:-Parameters-Table)

### Download
**NOTE: Currently it is preferable to use this method which is the most up to date.**<br/>
If you don't want to install it: You can use the AskoR.R file in [ScriptR](https://github.com/asusete/askoR/tree/master/ScriptR) instead and "run file" in the same directory.  
Just source it in your R script file (see _AskoR_analysis_script.R_):  
```
source("/directory/where/you/downloaded/the/file/AskoR.R")
```
OR

**NOTE: This method is not up to date, we will update it asap !**<br/>
Download the latest development code of AskoR from GitHub using [devtools](https://cran.r-project.org/package=devtools) with
```
install.packages("devtools")
library(devtools)
devtools::install_github("askomics/askoR")
```
For Windows users only: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) or check that it is already installed (needed to build the package).
 
### License

The coseq package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at http://www.r-project.org/Licenses/GPL-3.
