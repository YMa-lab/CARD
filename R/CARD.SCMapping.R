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
####################################################################################################
#' The function to estimate the cell type composition signature for each single cell in the scRNaseq reference data
#'
#' @param sc_eset the sc_eset stored in the CARD object
#' @param ct.varname character, the name of the column in metaData that specifies the cell type annotation information, stored in the CARD object
#' @param ct.select vector of cell type names that you are interested in to deconvolute, default as NULL. stored in the CARD object
#' @param sample.varname character,the name of the column in metaData that specifies the sample information. stored in the CARD object
#' @param B reference basis matrix stored in the CARD object.

#' @importFrom SummarizedExperiment assays
#' @importFrom nnls nnls

#' @return Returns a matrix of the cell type composition signature for each single cell in the scRNaseq reference
#'
getWeightForCell <- function(sc_eset,ct.varname,ct.select,sample.varname,B){
   count = assays(sc_eset)$counts
   count = count[,colSums(count) > 0]
   count = sweep(count,2,colSums(count),"/")
   count = count[rownames(count) %in% rownames(B),]
   count = count[match(rownames(B),rownames(count)),]
   Mean_Cell <- sapply(1:ncol(count),function(icell){
        mod1 = nnls(as.matrix(B),as.matrix(count[,icell]))
        return(mod1$x)}
        )
   Mean_Cell = t(Mean_Cell)
   rownames(Mean_Cell) = colnames(count)
   colnames(Mean_Cell) = colnames(B) 
   return(Mean_Cell)
}
#' The function to sample the spatial location information for each single cell
#'
#' @param Cords The spatial location information in the measure spatial locations, with the first and second columns represent the 2-D x-y coordinate system
#' @param numCell a numeric value indicating the number of single cells in each measured location, we suggest 20 for ST technology, 7 for 10x Viisum and 2 for Slide-seq
#' @param shape a character indicating whether the sampled spatial coordinates for single cells locating in a Square-like region or a Circle-like region. The center of this region is the measured spatial location in the non-single cell resolution spatial transcriptomics data. The default is "Square", the other shape is "Circle"

#' @import spatstat
#' @importFrom fields rdist
#' @importFrom stats runif
#' @return Returns a dataframe with the sampled spatial location information for each single cell 
#'
getHighresCords = function(Cords,numCell,shape = "Square"){ 
    ED <- rdist(as.matrix(Cords))
    n=dim(Cords)[1]    
    dis = c()
    for(i in 1:dim(ED)[1]){
        dis[i] = min(ED[i,-i]) # not count the cell it self, only use distance to its nearest cell
    }
    min_distance = median(dis)/2
    set.seed(20210107) 
    Cords_new = NULL
    for(i in 1:nrow(Cords)){   
        getPointsWithinCircle <- function(ED,i,numCell,min_distance){
            #dis = min(ED[i,-i])
            #min_distance = dis / 2.0
            circle = runifdisc(numCell, radius=min_distance, centre=c(Cords[i,1],Cords[i,2]),nsim=1, drop=TRUE)
            df = data.frame(x = circle$x, y = circle$y)
            return(df)
        }
        getPointsWithinSquare <- function(Cords,i,numCell,min_distance){
        minX = min_distance
        minY = min_distance
        rectangular_x = runif(numCell, min=Cords[i,1] - minX, max=Cords[i,1] + minX)
        rectangular_y = runif(numCell, min=Cords[i,2] - minY, max=Cords[i,2] + minY)
        df = data.frame(x = rectangular_x,y=rectangular_y)
        return(df)
        }
        if(shape == "Square"){
    df = getPointsWithinSquare(Cords,i,numCell,min_distance)
    }else if(shape == "Circle"){
    df = getPointsWithinCircle(ED,i,numCell,min_distance)
    }
    colnames(df) <- c("x","y")   
    df$centerSPOT = paste0(Cords[i,1],"x",Cords[i,2])
    df$centerx = Cords[i,1]
    df$centery = Cords[i,2]
    Cords_new = rbind(Cords_new,df)
}
    Cords_new = Cords_new[!duplicated(paste0(Cords_new$x,"x",Cords_new$y)),]
    rownames(Cords_new) = paste0(Cords_new$x,"x",Cords_new$y)
    return(Cords_new)
}

#' The function to assign the spatial location information for each single cell
#'
#' @param MappintSpotCellCor a mapped correlation matrix indicating the relashionship between each measured spatial location and the single cell in the scRNAseq reference 
#' @param Cords_new output from the function getHighresCords
#' @param numCell a numeric value indicating the number of single cells in each measured location, we suggest 20 for ST technology, 7 for 10x Viisum and 2 for Slide-seq
#' @param sc_eset a single cell experiment object stored in CARD object
#' @param ct.varname character, the name of the column in metaData that specifies the cell type annotation information, stroed in CARD object

#' @importFrom SingleCellExperiment colData
#' @return Return the assigned spatial location information for the mapped single cell 
#'
AssignSCcords <- function(MappintSpotCellCor,Cords_new,numCell,sc_eset,ct.varname){
MapCellCords = NULL
for(ispot in 1:nrow(MappintSpotCellCor)){
    mapCell = MappintSpotCellCor[ispot,]
    mapCell = mapCell[order(mapCell,decreasing = T)]
    centerspot = rownames(MappintSpotCellCor)[ispot]
    ispot_cords = data.frame(x=as.numeric(sapply(strsplit(centerspot,split="x"),"[",1)),y=as.numeric(sapply(strsplit(centerspot,split="x"),"[",2)))
    subCordsCell = Cords_new[Cords_new$centerSPOT == centerspot,]
    ##### calculate Euclidiean distance with this
    EDwithCenter = rdist(as.matrix(subCordsCell[,1:2]),as.matrix(ispot_cords))
    subCordsCell$EDwithCenter = EDwithCenter
    subCordsCell = subCordsCell[order(subCordsCell$EDwithCenter,decreasing = F),]
    mapCellCordsTemp = data.frame(CorwithSpot = mapCell[1:numCell])
    mapCellCordsTemp = cbind(mapCellCordsTemp,subCordsCell)
    mapCellCordsTemp$CT = colData(sc_eset)[rownames(mapCellCordsTemp),ct.varname]
    mapCellCordsTemp$Cell = rownames(mapCellCordsTemp)
    MapCellCords = rbind(MapCellCords,mapCellCordsTemp)
}
return(MapCellCords)
}

#' Extension of CARD into performing single cell Mapping from non-single cell spatial transcriptomics dataset. 
#'
#' @param CARD_object CARD object create by the createCARDObject function.
#' @param shapeSpot a character indicating whether the sampled spatial coordinates for single cells locating in a Square-like region or a Circle-like region. The center of this region is the measured spatial location in the non-single cell resolution spatial transcriptomics data. The default is "Square", the other shape is "Circle"
#' @param numCell a numeric value indicating the number of single cells in each measured location, we suggest 20 for ST technology, 7 for 10x Viisum and 2 for Slide-seq
#' @param ncore a numeric value indicating the number of cores used to accelerating the procedure

#' @importFrom SummarizedExperiment assays colData
#' @importFrom pbmcapply pbmclapply
#' @return Returns a SingleCellExperiment SCE object with the mapped expression at single cell resolution and the spatial location information of each single cell
#'
#' @export
#'
CARD_SCMapping = function(CARD_object,shapeSpot="Square",numCell,ncore = 10){
### load in spatial transcriptomics data stored in CARD_object
sc_eset = CARD_object@sc_eset
B = CARD_object@algorithm_matrix$B
ct.select = CARD_object@info_parameters$ct.select
ct.varname = CARD_object@info_parameters$ct.varname
sample.varname = CARD_object@info_parameters$sample.varname
res_CARD = CARD_object@Proportion_CARD
res_CARD = res_CARD[,order(colnames(res_CARD))]
Mean_Cell = getWeightForCell(sc_eset,ct.varname,ct.select,sample.varname,B)
Mean_Cell = Mean_Cell[,order(colnames(Mean_Cell))]
Cords = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(res_CARD),split="x"),"[",1)),y=as.numeric(sapply(strsplit(rownames(res_CARD),split="x"),"[",2)))
rownames(res_CARD) = paste0(Cords$x,"x",Cords$y)
Cords_new = getHighresCords(Cords,numCell = numCell,shape = shapeSpot)
MappintSpotCellCor = cor(t(res_CARD),t(Mean_Cell))
rownames(MappintSpotCellCor) = paste0(Cords$x,"x",Cords$y)
############ high resolution coordinates
MapCellCords = AssignSCcords(MappintSpotCellCor,Cords_new,numCell,sc_eset,ct.varname)
########### map scRNAseq count data 
count_sc = assays(sc_eset)$counts
count_sc = count_sc[,colSums(count_sc) > 0]
count_CT = NULL
count_CT = pbmclapply(1:nrow(res_CARD),function(ispot){
spot = rownames(res_CARD)[ispot]
    MapCellCords_spot = MapCellCords[MapCellCords$centerSPOT == spot,]
    df = as(count_sc[,MapCellCords_spot$Cell],"sparseMatrix")
    colnames(df) = paste0(MapCellCords_spot$Cell,":",spot)
    colnames(df) = paste0(colnames(df),":",MapCellCords_spot$x,"x",MapCellCords_spot$y)
    df
},mc.cores = ncore,mc.set.seed = F)
count_CT = do.call("cbind",count_CT)
sce <- SingleCellExperiment(list(counts=count_CT),
    colData=as.data.frame(MapCellCords[,c("x","y","centerSPOT","centerx","centery","CT","Cell")]),
    rowData=as.data.frame(rownames(count_CT)))
return(sce)
}

