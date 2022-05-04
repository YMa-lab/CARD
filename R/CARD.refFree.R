####################################################################################################
## Package : CARD
## Version : 1.0.1
## Date    : 2021-1-7 09:10:08
## Modified: 2021-12-13 16:18:07
## Title   : Spatially Informed Cell Type Deconvolution for Spatial Transcriptomics by CARD.
## Authors : Ying Ma
## Contacts: yingma@umich.edu
##           University of Michigan, Department of Biostatistics
####################################################################################################

#' Extension of CARD into a reference-free version of deconvolution: CARDfree. 
#'
#' @param CARDfree_object CARDfree object create by the createCARDfreeObject function.
#' 
#' @importFrom Rcpp sourceCpp
#' @importFrom fields rdist

#' @return Returns a CARD object with estimated cell type proportion stored in object@Proportion_CARD. Because this is a reference-free version, the columns of estimated proportion is not cell type but cell type cluster
#'
#' @export
#'
CARD_refFree = function(CARDfree_object){
### load in spatial transcriptomics data stored in CARDfree_object
spatial_countMat = CARDfree_object@spatial_countMat
spatial_location = CARDfree_object@spatial_location
### load in markerList
## number of cell type clusters
markerList = CARDfree_object@markerList
numK = length(markerList)
marker = unique(unlist(markerList))
#marker = toupper(marker)

cat(paste0("## Number of unique marker genes: ",length(marker)," for ",numK," cell types ...\n"))
commonGene = intersect(toupper(rownames(spatial_countMat)),toupper(marker))
#### remove mitochondrial and ribosomal genes
commonGene  = commonGene[!(commonGene %in% commonGene[grep("mt-",commonGene)])]
if(length(commonGene) < numK * 10){
  cat(paste0("## STOP! The average number of unique marker genes for each cell type is less than 20 ...\n"))
}
Xinput = spatial_countMat[order(rownames(spatial_countMat)),]
Xinput = Xinput[order(rownames(Xinput)),]
Xinput = Xinput[toupper(rownames(Xinput)) %in% commonGene,]
Xinput = Xinput[rowSums(Xinput) > 0,]
Xinput = Xinput[,colSums(Xinput) > 0]
Xinput_norm = sweep(Xinput,2,colSums(Xinput),"/")

#### initialization
if(ncol(Xinput_norm) < 5000){
suppressMessages(require(NMF))
sink("NUL")
set.seed(20200107)
NMFout <- invisible(nmf(as.matrix(Xinput_norm),numK))
sink()
B = NMFout@fit@W
Vint1 = as.matrix(t(NMFout@fit@H))
rownames(Vint1) = colnames(Xinput_norm)
}else{
suppressMessages(require(RcppML))
sink("NUL")
set.seed(20200107)
NMFout <- invisible(RcppML::nmf(as.matrix(Xinput_norm),numK))
sink()
B = NMFout$w
rownames(B) = rownames(Xinput_norm)
Vint1 = as.matrix(t(NMFout$h))
rownames(Vint1) = colnames(Xinput_norm)
}
spatial_location = spatial_location[rownames(spatial_location) %in% colnames(Xinput_norm),]
spatial_location = spatial_location[match(colnames(Xinput_norm),rownames(spatial_location)),]
norm_cords = spatial_location[ ,c("x","y")]
norm_cords$x = norm_cords$x - min(norm_cords$x)
norm_cords$y = norm_cords$y - min(norm_cords$y)
scaleFactor = max(norm_cords$x,norm_cords$y)
norm_cords$x = norm_cords$x / scaleFactor
norm_cords$y = norm_cords$y / scaleFactor
### Euclidiean distance
ED <- rdist(as.matrix(norm_cords))##Euclidean distance matrix
b = rep(0,ncol(B))
###### parameters that need to be set
isigma = 0.1 ####construct Gaussian kernel with the default scale /length parameter to be 0.1
epsilon = 1e-04  #### convergence epsion 
phi = c(0.01,0.1,0.3,0.5,0.7,0.9,0.99) #### grided values for phi
kernel_mat <- exp(-ED^2 / (2 * isigma^2))
diag(kernel_mat) <- 0

rm(ED)
rm(Xinput)
rm(norm_cords)
gc()
###### scale the Xinput_norm and B to speed up the convergence. 
mean_X = mean(Xinput_norm)
mean_B = mean(B)
Xinput_norm = Xinput_norm * 1e-01 / mean_X
B = B * 1e-01 / mean_B
gc()
ResList = list()
Obj = c()
for(iphi in 1:length(phi)){
res = CARDfree(
  XinputIn = as.matrix(Xinput_norm),
  UIn = as.matrix(B),
  WIn = kernel_mat, 
  phiIn = phi[iphi],
  max_iterIn =1000,
  epsilonIn = epsilon,
  initV = Vint1,
  initb = rep(0,ncol(B)),
  initSigma_e2 = 0.1, 
  initLambda = rep(10,ncol(B)))
rownames(res$V) = colnames(Xinput_norm)
colnames(res$V) = paste0("CT",1:ncol(B))
ResList[[iphi]] = res
Obj = c(Obj,res$Obj)
}
Optimal = which(Obj == max(Obj))
Optimal = Optimal[length(Optimal)] #### just in case if there are two equal objective function values
OptimalPhi = phi[Optimal]
OptimalRes = ResList[[Optimal]]
cat(paste0("## Deconvolution Finish! ...\n"))
CARDfree_object@info_parameters$phi = OptimalPhi
CARDfree_object@Proportion_CARD = sweep(OptimalRes$V,1,rowSums(OptimalRes$V),"/")
CARDfree_object@algorithm_matrix = list(B = OptimalRes$B, Xinput_norm = Xinput_norm * mean_X / 1e-01, Res = OptimalRes)
CARDfree_object@spatial_location = spatial_location
return(CARDfree_object)
}

