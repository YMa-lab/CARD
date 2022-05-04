####################################################################################################
## Package : CARD
## Version : 1.0.1
## Date    : 2021-1-7 09:10:08
## Modified: 2021-5-20 15:25:07
## Title   : Spatially Informed Cell Type Deconvolution for Spatial Transcriptomics by CARD.
## Authors : Ying Ma
## Contacts: yingma@umich.edu
##           University of Michigan, Department of Biostatistics
####################################################################################################

#' Construct the mean gene expression basis matrix (B), this is the faster version 
#'
#' @param x S4 class for storing data from single-cell experiments. This format is usually created by the package SingleCellExperiment with stored counts, along with the usual metadata for genes and cells.
#' @param ct.select vector of cell type names that you are interested in to deconvolute, default as NULL. If NULL, then use all cell types provided by single cell dataset;
#' @param ct.varname character, the name of the column in metaData that specifies the cell type annotation information
#' @param sample.varname character,the name of the column in metaData that specifies the sample information. If NULL, we just use the whole as one sample.
#'
#' @importFrom SummarizedExperiment assays colData
#' @importFrom wrMisc rowGrpMeans
#' @return Return a list of basis (B) matrix
#'
#' @export
#'
createscRef <- function(x, ct.select = NULL, ct.varname, sample.varname = NULL){
if (is.null(ct.select)) {
        ct.select <- unique(colData(x)[, ct.varname])
}
ct.select <- ct.select[!is.na(ct.select)]
countMat <- as(assays(x)$counts,"sparseMatrix")
ct.id <- droplevels(as.factor(colData(x)[, ct.varname]))
if(is.null(sample.varname)){
  colData(x)$sampleID = "Sample"
  sample.varname = "sampleID"
}
sample.id <- as.character(colData(x)[, sample.varname])
ct_sample.id <- paste(ct.id, sample.id, sep = "$*$")
colSums_countMat <- colSums(countMat)
colSums_countMat_Ct = aggregate(colSums_countMat ~ ct.id + sample.id, FUN = 'sum')
colSums_countMat_Ct_wide = reshape(colSums_countMat_Ct, idvar = "sample.id", timevar = "ct.id", direction = "wide")
colnames(colSums_countMat_Ct_wide) = gsub("colSums_countMat.","",colnames(colSums_countMat_Ct_wide))
rownames(colSums_countMat_Ct_wide) = colSums_countMat_Ct_wide$sample.id
colSums_countMat_Ct_wide$sample.id <- NULL
tbl <- table(sample.id,ct.id)
colSums_countMat_Ct_wide = colSums_countMat_Ct_wide[,match(colnames(tbl),colnames(colSums_countMat_Ct_wide))]
colSums_countMat_Ct_wide = colSums_countMat_Ct_wide[match(rownames(tbl),rownames(colSums_countMat_Ct_wide)),]

S_JK <- colSums_countMat_Ct_wide / tbl
S_JK <- as.matrix(S_JK)
S_JK[S_JK == 0] = NA
S_JK[!is.finite(S_JK)] = NA
S = colMeans(S_JK, na.rm = TRUE)
S = S[match(unique(ct.id),names(S))]
#library("wrMisc")
# Theta_S <- sapply(unique(ct_sample.id),ćfunction(ict.sample) {
#             y = countMat[, ct_sample.id %in% ict.sample,drop = F]
#             rsm = rowSums(y)
#             rsm / sum(rsm)
#             })
Theta_S_rowMean <- rowGrpMeans(as.matrix(countMat), grp = ct_sample.id, na.rm = TRUE)
tbl_sample = table(ct_sample.id)
tbl_sample = tbl_sample[match(colnames(Theta_S_rowMean),names(tbl_sample))]
Theta_S_rowSums <- sweep(Theta_S_rowMean,2,tbl_sample,"*")
Theta_S <- sweep(Theta_S_rowSums,2,colSums(Theta_S_rowSums),"/")
grp <- sapply(strsplit(colnames(Theta_S),split="$*$",fixed = TRUE),"[",1)
Theta = rowGrpMeans(Theta_S, grp = grp, na.rm = TRUE)
Theta = Theta[,match(unique(ct.id),colnames(Theta))]
S = S[match(colnames(Theta),names(S))]
basis = sweep(Theta,2,S,"*")
colnames(basis) = colnames(Theta)
rownames(basis) = rownames(Theta)
return(list(basis = basis))
}

#' Select Informative Genes used in the deconvolution
#'
#' @param Basis Reference basis matrix.
#' @param sc_eset scRNAseq data along with meta data stored in the S4 class format (SingleCellExperiment).
#' @param commonGene common genes between scRNAseq count data and spatial resolved transcriptomics data.
#' @param ct.select vector of cell type names that you are interested in to deconvolute, default as NULL. If NULL, then use all cell types provided by single cell dataset;
#' @param ct.varname character, the name of the column in metaData that specifies the cell type annotation information
#' 
#' @importFrom SummarizedExperiment assays colData
#' @importFrom stats quantile
#' @return a vector of informative genes selected
#'
#' @export
#' 
selectInfo <- function(Basis,sc_eset,commonGene,ct.select,ct.varname){
#### log2 mean fold change >0.5
gene1 = lapply(ct.select,function(ict){
rest = rowMeans(Basis[,colnames(Basis) != ict])
FC = log((Basis[,ict] + 1e-06)) - log((rest + 1e-06))
rownames(Basis)[FC > 1.25 & Basis[,ict] > 0]
})
gene1 = unique(unlist(gene1))
gene1 = intersect(gene1,commonGene)
counts = assays(sc_eset)$counts
counts = counts[rownames(counts) %in% gene1,]
sd_within = sapply(ct.select,function(ict){
  temp = counts[,colData(sc_eset)[,ct.varname] == ict]
  apply(temp,1,var) / apply(temp,1,mean)
  })
##### remove the outliers that have high dispersion across cell types
gene2 = rownames(sd_within)[apply(sd_within,1,mean,na.rm = T) < quantile(apply(sd_within,1,mean,na.rm = T),prob = 0.99)]
return(gene2)
}

#' Spatially Informed Cell Type Deconvolution for Spatial Transcriptomics by CARD
#'
#' @param CARD_object CARD object create by the createCARDObject function
#' 
#' @importFrom Rcpp sourceCpp
#' @importFrom MCMCpack rdirichlet
#' @importFrom fields rdist
#' @return Returns a CARD object with estimated cell type proportion stored in object@Proportion_CARD.
#'
#' @export
#' 
CARD_deconvolution <- function(CARD_object){
ct.select = CARD_object@info_parameters$ct.select
ct.varname = CARD_object@info_parameters$ct.varname
sample.varname = CARD_object@info_parameters$sample.varname
cat(paste0("## create reference matrix from scRNASeq...\n"))
sc_eset = CARD_object@sc_eset
Basis_ref = createscRef(sc_eset, ct.select, ct.varname, sample.varname)
Basis = Basis_ref$basis
Basis = Basis[,colnames(Basis) %in% ct.select]
Basis = Basis[,match(ct.select,colnames(Basis))]
spatial_count = CARD_object@spatial_countMat
commonGene = intersect(rownames(spatial_count),rownames(Basis))
#### remove mitochondrial and ribosomal genes
commonGene  = commonGene[!(commonGene %in% commonGene[grep("mt-",commonGene)])]
cat(paste0("## Select Informative Genes! ...\n"))
common = selectInfo(Basis,sc_eset,commonGene,ct.select,ct.varname)
Xinput = spatial_count
rm(spatial_count)
B = Basis
rm(Basis)
##### match the common gene names
Xinput = Xinput[order(rownames(Xinput)),]
B = B[order(rownames(B)),]
B = B[rownames(B) %in% common,]
Xinput = Xinput[rownames(Xinput) %in% common,]
##### filter out non expressed genes or cells again
Xinput = Xinput[rowSums(Xinput) > 0,]
Xinput = Xinput[,colSums(Xinput) > 0]
##### normalize count data
colsumvec = colSums(Xinput)
Xinput_norm = sweep(Xinput,2,colsumvec,"/")
B = B[rownames(B) %in% rownames(Xinput_norm),]    
B = B[match(rownames(Xinput_norm),rownames(B)),]
#### spatial location
spatial_location = CARD_object@spatial_location
spatial_location = spatial_location[rownames(spatial_location) %in% colnames(Xinput_norm),]
spatial_location = spatial_location[match(colnames(Xinput_norm),rownames(spatial_location)),]

##### normalize the coordinates without changing the shape and relative position
norm_cords = spatial_location[ ,c("x","y")]
norm_cords$x = norm_cords$x - min(norm_cords$x)
norm_cords$y = norm_cords$y - min(norm_cords$y)
scaleFactor = max(norm_cords$x,norm_cords$y)
norm_cords$x = norm_cords$x / scaleFactor
norm_cords$y = norm_cords$y / scaleFactor
##### initialize the proportion matrix
ED <- rdist(as.matrix(norm_cords))##Euclidean distance matrix
cat(paste0("## Deconvolution Starts! ...\n"))
set.seed(20200107)
Vint1 = as.matrix(rdirichlet(ncol(Xinput_norm), rep(10,ncol(B))))
colnames(Vint1) = colnames(B)
rownames(Vint1) = colnames(Xinput_norm)
b = rep(0,length(ct.select))
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
res = CARDref(
  XinputIn = as.matrix(Xinput_norm),
  UIn = as.matrix(B),
  WIn = kernel_mat, 
  phiIn = phi[iphi],
  max_iterIn =1000,
  epsilonIn = epsilon,
  initV = Vint1,
  initb = rep(0,ncol(B)),
  initSigma_e2 = 0.1, 
  initLambda = rep(10,length(ct.select)))
rownames(res$V) = colnames(Xinput_norm)
colnames(res$V) = colnames(B)
ResList[[iphi]] = res
Obj = c(Obj,res$Obj)
}
Optimal = which(Obj == max(Obj))
Optimal = Optimal[length(Optimal)] #### just in case if there are two equal objective function values
OptimalPhi = phi[Optimal]
OptimalRes = ResList[[Optimal]]
cat(paste0("## Deconvolution Finish! ...\n"))
CARD_object@info_parameters$phi = OptimalPhi
CARD_object@Proportion_CARD = sweep(OptimalRes$V,1,rowSums(OptimalRes$V),"/")
CARD_object@algorithm_matrix = list(B = B * mean_B / 1e-01, Xinput_norm = Xinput_norm * mean_X / 1e-01, Res = OptimalRes)
CARD_object@spatial_location = spatial_location
return(CARD_object)
}



