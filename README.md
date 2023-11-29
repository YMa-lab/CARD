# CARD

## Spatially Informed Cell Type Deconvolution for Spatial Transcriptomics 

![CARD\_pipeline](Overview.jpg)
We developed a statistical method for spatially informed cell type deconvolution for spatial transcriptomics. Briefly,CARD is a reference-based deconvolution method that estimates cell type composition in spatial transcriptomics based on cell type specific expression information obtained from a reference scRNA-seq data. A key feature of CARD is its ability to accommodate spatial correlation in the cell type composition across tissue locations, enabling accurate and spatially informed cell type deconvolution as well as refined spatial map construction. CARD relies on an efficient optimization algorithm for constrained maximum likelihood estimation and is scalable to spatial transcriptomics with tens of thousands of spatial locations and tens of thousands of genes. CARD is implemented as an open-source R package, freely available at www.xzlab.org/software.html. 

Installation
------------
You can install the released version of CARD from Github with the following code, for more installation details or solutions that might solve related issues (specifically MacOS system) see the [link](https://yma-lab.github.io/CARD/documentation/02_installation.html).

## Dependencies 
* R version >= 4.0.0.
* R packages: SingleCellExperiment, SummarizedExperiment, concaveman, sp, Matrix, methods, ggplot2, ggcorrplot, MuSiC, fields, MCMCpack, dplyr, sf, RANN, stats, reshape2, RColorBrewe, scatterpie, grDevices, stats, nnls, pbmcapply, spatstat, gtools, RcppML, NMF

``` r
# install devtools if necessary
install.packages('devtools')

# install the CARD package
devtools::install_github('YingMa0107/CARD')

# load package
library(CARD)

```
The R package has been installed successfully on Operating systems: 
* macOS Catalina 10.15, macOS Monterey 12.3.1
* Ubuntu 18.04.5 LTS (Bionic Beaver) 
* Windows 10

# Issues
All feedback, bug reports and suggestions are warmly welcomed! Please make sure to raise issues with a detailed and reproducible exmple and also please provide the output of your sessionInfo() in R! 

How to cite `CARD`
-------------------
Ma, Y., Zhou, X. Spatially informed cell-type deconvolution for spatial transcriptomics. Nat Biotechnol 40, 1349–1359 (2022). https://doi.org/10.1038/s41587-022-01273-7

How to use `CARD`
-------------------
Details in [Tutorial](https://yma-lab.github.io/CARD/)
