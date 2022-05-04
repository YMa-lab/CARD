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

#' Make new spatial locations on unmeasured tissue through grids. 
#'
#' @param location Data frame, spatial location data frame of the original spatial resolved transcriptomics dataset, stored in the CARD_object@spatial_location
#' @param numSample Numeric, approximate number of cells in grid within the shape of the spatial location data frame
#' @param concavity Numeric, a relative measure of concavity. The default is 2.0, which can prodecure detailed enough shapes. Infinity results in a convex hull while 1 results in a more detailed shape.
#' 
#' @import concaveman
#' @import sp 
#' @importFrom dplyr '%>%'
#' @importFrom sf st_as_sf
#' @return Return a list of data frame with newly grided points
#'
#' @export
#'

sampleGridWithin <- function(location,numSample,concavity = 2){
###### recognise the outline of the shape and extract the coordinates of the polygons
df = location[,1:2]
pnts <- df %>% st_as_sf(coords = c("x", "y"))
polygon <- concaveman(pnts,concavity) #### tune the number, default is 2.0, is enough. 
poly_coords = as.data.frame(as.matrix(polygon$polygons[[1]]))
colnames(poly_coords) <- c("x","y")
coordinates(poly_coords) <- ~x + y
###### create the spatial polygon data frame
crds <- poly_coords
# str(crds)
Pl <- Polygon(crds)
# str(Pl)
ID <- "21*15"
Pls <- Polygons(list(Pl), ID=ID)
# str(Pls)
SPls <- SpatialPolygons(list(Pls))
# str(SPls)
df <- data.frame(value=1, row.names=ID)
# str(df)
SPDF <- SpatialPolygonsDataFrame(SPls, df)
###### random make the grids within the polygon
pts <- makegrid(SPDF, numSample)
pts1 <- SpatialPoints(pts)
###### extract the points within the locations
spgrdWithin <- SpatialPixels(pts1[SPDF,])
data = as.data.frame(spgrdWithin@coords)
data$id = 1:nrow(data)
colnames(data) <- c("x","y")
return(data)
}

#' Normalize the new spatial locations without changing the shape and relative positions
#'
#' @param locationOrig Data frame, spatial location data frame of the original spatial resolved transcriptomics dataset, stored in the CARD_object@spatial_location
#' @param trainInd Vector, Index of the original spatial locations 
#' @param testInd Vector, Index of the newly grided spatial locations
#' 
#' @return Return the normalized spatial location data frame
#'
#' @export
#'
normCoordsTrainTest <- function(locationOrig,trainInd,testInd){ ### normalize to 0-1 scale
   norm_coords_train = as.data.frame(locationOrig[trainInd,1:2])
   locatiooFactor_x = min(norm_coords_train$x)
   locatiooFactor_y = min(norm_coords_train$y)
   norm_coords_train$x = norm_coords_train$x - locatiooFactor_x
   norm_coords_train$y = norm_coords_train$y - locatiooFactor_y
   scaleFactor = max(norm_coords_train$x,norm_coords_train$y)
   norm_coords_train$x = norm_coords_train$x / scaleFactor
   norm_coords_train$y = norm_coords_train$y / scaleFactor
   norm_coords_test = as.data.frame(locationOrig[testInd,1:2])
   norm_coords_test[,1] = (norm_coords_test[,1] - locatiooFactor_x) / scaleFactor 
   norm_coords_test[,2] = (norm_coords_test[,2] - locatiooFactor_y) / scaleFactor
   locationTemp = rbind(norm_coords_test,norm_coords_train)
   return(locationTemp)
}

#' Calculate the variance covariance matrix used in the imputation of the new grided locations
#'
#' @param locationOrig Data frame, spatial location data frame of the original spatial resolved transcriptomics dataset, stored in the CARD_object@spatial_location
#' @param trainInd Vector, index of the original spatial locations 
#' @param testInd Vector, index of the newly grided spatial locations
#' @param OptimalPhi Numeric, the optimal phi value stored in CARD_object
#' @param ineibor Numeric, number of neighbors used in the imputation on newly grided spatial locations, default is 10.
#' 
#' @import Matrix
#' @importFrom RANN nn2 
#' @return Return a list with the imputed Cell type composition Vtest matrix on the newly grided spatial locations and predicted normalized gene expression
#'
#' @export
#'

Sigma <- function(locationOrig,trainInd,testInd,OptimalPhi,ineibor){
norm_cordsTemp = normCoordsTrainTest(locationOrig,trainInd,testInd) 
##### find neiibors
near_data = nn2(norm_cordsTemp[ ,1:2],k = ineibor + 1) ###### nn2 will consider the location itself as a neighbor
neibors = near_data$nn.idx
neibors = neibors[,-1] ##### delete the location itself as the neighbor
Nmat = Matrix(0,nrow = nrow(neibors),ncol = nrow(neibors),sparse = TRUE)
##### speed up the matrices when it is large
for(icol in 1:ncol(neibors)){
edges = data.frame(i = 1:nrow(neibors), j = neibors[,icol])
adjacency = sparseMatrix(i = as.integer(edges$i),
                          j = as.integer(edges$j),
                          x = 1,
                          dims = rep(nrow(neibors), 2),
                          use.last.ij = TRUE)
Nmat = Nmat + adjacency
}
##### find mutual neibors, for non-mutual neighbors, the value will be zero
Nmat = Nmat * t(Nmat)
##### create euclidiean distance
dist_neibors = near_data$nn.dists
isigma = 0.1
kernel_neibors = exp(-dist_neibors^2/(2*isigma^2))
WTemp <- Matrix(0,nrow = nrow(neibors),ncol = nrow(neibors),sparse = TRUE)
##### speed up the matrices when it is large
for(icol in 1:ncol(kernel_neibors)){
temp = near_data$nn.idx
edges = data.frame(i = 1:nrow(temp), j = temp[,icol])

adjacency = sparseMatrix(i = as.integer(edges$i),
                          j = as.integer(edges$j),
                          x = kernel_neibors[,icol],
                          dims = rep(nrow(neibors), 2),
                          use.last.ij = TRUE)
WTemp = WTemp + adjacency
}
##### create euclidiean distance for mutual neibors, for non-mutual neighbors, the value will be zero
WTemp = WTemp * Nmat
rownames(WTemp) = colnames(WTemp) = rownames(norm_cordsTemp)
##### calculate variance matrix that I need
DTemp = Diagonal(x = colSums(WTemp))
In = Diagonal(nrow(WTemp ))
W22 = WTemp[(length(testInd)+1):nrow(WTemp),(length(testInd)+1):ncol(WTemp)]
SigmaTemp = DTemp - OptimalPhi * WTemp
Sigma11 =  SigmaTemp[1:length(testInd),1:length(testInd)]
Sigma12 =  SigmaTemp[1:length(testInd),(length(testInd)+1):ncol(SigmaTemp)]
Sigma21 =  SigmaTemp[(length(testInd)+1):nrow(SigmaTemp),1:length(testInd)]
Sigma22 =  SigmaTemp[(length(testInd)+1):nrow(SigmaTemp),(length(testInd)+1):ncol(SigmaTemp)]
return(list(SigmaTemp = SigmaTemp,Sigma11 = Sigma11,Sigma12 = Sigma12,Sigma21 = Sigma21,Sigma22 = Sigma22,W22 = W22))
}

#' Imputation and Construction of High-Resolution Spatial Maps for Cell Type Composition and Gene Expression by the spatial correlation structure between original spatial locations and new grided spatial locations
#'
#' @param Vtrain Matrix, estimated V matrix from CARD 
#' @param locationOrig Data frame, spatial location data frame of the original spatial resolved transcriptomics dataset, stored in the CARD_object@spatial_location
#' @param trainInd Vector, index of the original spatial locations 
#' @param testInd Vector, index of the newly grided spatial locations
#' @param B Matrix, used in the deconvolution as the reference basis matrix
#' @param Xinput_norm Matrix, used in the deconvolution as the normalized spatial count data
#' @param Optimalb Vector, vector of the intercept for each cel type estimated based on the original spatial resolution
#' @param OptimalPhi Numeric, the optimal phi value stored in CARD_object
#' @param lambda Vector, vector of cell type specific scalar in the CAR model 
#' @param ineibor Numeric, number of neighbors used in the imputation on newly grided spatial locations, default is 10.
#' 
#' @return Return a list with the imputed Cell type composition Vtest matrix on the newly grided spatial locations and predicted normalized gene expression
#'
#' @export
#'
MVN_CV= function(Vtrain,locationOrig,trainInd,testInd,B,Xinput_norm,Optimalb,OptimalPhi,lambda,ineibor){
SigmaList = Sigma(locationOrig,trainInd,testInd,OptimalPhi,ineibor)
Sigma11 = SigmaList$Sigma11
Sigma21 = SigmaList$Sigma21
Sigma12 = SigmaList$Sigma12
Sigma22 = SigmaList$Sigma22
W22 = SigmaList$W22
Mu_cond = NULL
for(ict in 1:length(lambda)){
    temp = Sigma12 %*% (Vtrain[,ict] - Optimalb[1:nrow(Vtrain),ict])
    solver = solve(Sigma11,temp)
    Mu = Optimalb[1:length(testInd),ict] - solver  ### q * k
    Mu_cond = cbind(Mu_cond,Mu)
}
Vtest = Mu_cond
colnames(Vtest) = colnames(Vtrain)
XtestHat = B %*% t(Vtest)
rownames(Vtest) = colnames(XtestHat)
Vtest = sweep(Vtest,1,rowSums(Vtest),"/")
return(list(Vtest = Vtest, XtestHat = XtestHat))
#mean(mse(Xtest,XtestHat))
}

#' Construct an enhanced spatial expression map on the unmeasured tissue locations
#'
#' @param CARD_object CARD Object with estimated cell type compositions on the original spatial resolved transcriptomics data.
#' @param NumGrids Initial number of newly grided spatial locations. The final number of newly grided spatial locations will be lower than this value since the newly grided locations outside the shape of the tissue will be filtered
#' @param ineibor Numeric, number of neighbors used in the imputation on newly grided spatial locations, default is 10.
#' @param exclude Vector, the rownames of spatial location data on the original resolution that you want to exclude. This is to avoid the weird detection of the shape. 
#' 
#' @return Return CARD object with the refined cell type compositions estimated for newly grided spots and the refined predicted gene expression (normalized). 
#'
#' @export
#'
CARD.imputation <- function(CARD_object,NumGrids,ineibor = 10,exclude = NULL){
  ##### extract results from CARD object 
  B = CARD_object@algorithm_matrix$B
  Xinput_norm = CARD_object@algorithm_matrix$Xinput_norm
  Vtrain = CARD_object@algorithm_matrix$Res$V
  location <- CARD_object@spatial_location
  ##### check
  if(sum(rownames(location) == rownames(Vtrain)) == nrow(Vtrain)){
    cat(paste0("## The rownames of locations are matched ...\n"))
  }
  ##### Make the new grids
  cat(paste0("## Make grids on new spatial locations ...\n"))
  if(!is.null(exclude)){
    location_usetoSample <- location[!(rownames(location) %in% exclude),]
  }else{
    location_usetoSample <- location
  }
  data = sampleGridWithin(location_usetoSample,NumGrids,concavity = 2.0)
  rownames(data) = paste0(data$x,"x",data$y)
  ##### delete the newly grided locations that are the same as the original ones
  data = data[!(rownames(data) %in% rownames(location)),]
  location.train = location[,1:2]
  location.test = data[,1:2]
  locationCombine = rbind(data[,1:2],location[,1:2])
  testInd = 1:nrow(data)
  trainInd = (nrow(data)+1):nrow(locationCombine)
  locationOrig = locationCombine
  ##### results from the CARD
  Optimalb = CARD_object@algorithm_matrix$Res$b
  vecOne = rep(1,nrow(locationCombine))
  Optimalb = vecOne %*% t(Optimalb)
  OptimalPhi = CARD_object@info_parameters$phi
  lambda = CARD_object@algorithm_matrix$Res$lambda 
  #rm(res.spatialDeconv)
  imputation = MVN_CV(Vtrain,locationOrig,trainInd,testInd,B,Xinput_norm,Optimalb,OptimalPhi,lambda,ineibor)
  CARD_object@refined_prop = as.matrix(imputation$Vtest)
  CARD_object@refined_expression = as.matrix(imputation$XtestHat)
  return(CARD_object) 
}
