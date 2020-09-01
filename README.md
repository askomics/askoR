## AskoR 

Authors: Fabrice Legeai, Kevin Gazengel, Susete Alves Carvalho


**AskoR** is a pipeline for the analysis of gene expression data, using edgeR.
Several steps are performed: data filters (cpm method), normalize this filtered data, look at the correlation of our data, run differential expression analysis, compare contrast, GO enrichment and co-expression.

### Sample data
You'll find a test set in the **inst/extdata/input** folder. It'll be used for the vignettes documentation.<br/>

  - Count matrix file: CountsMatrix.txt **OR** Counts files per samples: files in "counts" directory _(counts.tgz)_
  - Samples file: Samples_CountsMatrix.txt **OR** Samples_CountsFiles.txt 
  - Contrasts file: Contrasts.txt
  - Genes annotations file: Genes_annotations.txt (optional)
  - GO annotations file: GO_annotations.txt (optional)

**IMPORTANT :** All input files must be in a folder named **input** _(case sensitive)_.

### Vignettes

* [User Guide]()
* [Parameters Table]()

### Download
Download the latest development code of AskoR from GitHub using [devtools](https://cran.r-project.org/package=devtools) with
```
install.packages("devtools")
library(devtools)
devtools::install_github("asusete/askoR")
```
### License

The coseq package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at http://www.r-project.org/Licenses/GPL-3.
