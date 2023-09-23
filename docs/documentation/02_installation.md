---
layout: page
title: Installation
description: ~
---

`CARD` is implemented as an R package, which can be installed from GitHub by:

### Dependencies 
* R version >= 4.0.0.
* R packages: SingleCellExperiment, SummarizedExperiment, concaveman, sp, Matrix, methods, ggplot2, ggcorrplot, MuSiC, fields, MCMCpack, dplyr, sf, RANN, stats, reshape2, RColorBrewe, scatterpie, grDevices, stats, nnls, pbmcapply, spatstat, gtools, RcppML, NMF


#### 1. Install `devtools` if necessary
```r
install.packages('devtools')
```

#### 2. Install `CARD`
```r
devtools::install_github('YMa-lab/CARD')
```
#### 3. Load package
```r
library(CARD)
```

This package is supported for Windows 10, MAC and Linux. The package has been tested on the following systems:
- Windows 10: Home (1903)
- MAC: OSX (10.14.1)
- Linux: Ubuntu (16.04.6)

#### 4. Some possible issues when installing the package, especially on the MacOS system
(1) Cannot find tools necessary when using R in the MacOS system.
```r
Error: Failed to install 'CARD' from GitHub:
  Could not find tools necessary to compile a package
Call `pkgbuild::check_build_tools(debug = TRUE)` to diagnose the problem.
``` 
possible solution: in R, type the code ``` options(buildtools.check = function(action) TRUE )```, see the discussion about this error, [link](https://stackoverflow.com/questions/37776377/error-when-installing-an-r-package-from-github-could-not-find-build-tools-neces)

(2) library not found for -lgfortran when using R in the MacOS system.
```r
ld: library not found for -lgfortran
```
It seems the gfortran is not well installed on the MacOS system. Please check this [link](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/) for the gfortran installation to see if it helps. 


