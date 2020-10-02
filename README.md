
<!-- README.md is generated from README.Rmd. Please edit that file -->
SCDC: Bulk Gene Expression Deconvolution by Multiple Single-Cell RNA Sequencing References
==========================================================================================

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/meichendong/SCDC.svg?branch=master)](https://travis-ci.org/meichendong/SCDC) [![CRAN status](https://www.r-pkg.org/badges/version/SCDC)](https://CRAN.R-project.org/package=SCDC) <!-- badges: end -->

SCDC is a [deconvolution](https://en.wikipedia.org/wiki/Deconvolution) method for bulk [RNA-seq](https://en.wikipedia.org/wiki/RNA-Seq) that leverages cell-type specific gene expressions from **multiple [scRNA-seq](https://en.wikipedia.org/wiki/Single_cell_sequencing) reference datasets**. SCDC adopts an ENSEMBLE method to integrate deconvolution results from different scRNA-seq datasets that are produced in different laboratories and at different times, implicitly addressing the [batch-effect](http://www.molmine.com/magma/global_analysis/batch_effect.html) confounding.

![SCDC framework](framework.PNG)

Installation
------------

You can install the released version of SCDC from [GitHub](https://github.com/) with:

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("meichendong/SCDC")
```

Vignettes
---------

Please see the [vignettes page](https://meichendong.github.io/SCDC/articles/SCDC.html).

The SCDC paper is published at [Briefings In Bioinformatics](https://doi.org/10.1093/bib/bbz166).


Note of the current repository
---------

This repository is a clone of https://github.com/crhisto/SCDC. It contains the following modifications: 
   - Compatibility with sparse matrices using: `dgCMatrix` objects in R.
   - Routines with parallelization
   - Dynamic threshold for markers selection
   - Improvements in logs and so on.

This has been done as part of the project: https://github.com/crhisto/thymus_NPM-ALK_notebook. 

If you want to install the SCDC library with these modifications you can use: 

``` r
if("SCDC" %in% rownames(installed.packages())){
  library(xbioc)
}else{
  devtools::install_github( repo = "crhisto/SCDC")
  library(SCDC)
}
```

Also, you must use the following libraries with support of sparse matrix (`dgCMatrix`) for large scRNA-seq datasets: 
- Biobase: [https://github.com/crhisto/Biobase](https://github.com/crhisto/Biobase)
- xbioc: [https://github.com/crhisto/xbioc](https://github.com/crhisto/xbioc)
