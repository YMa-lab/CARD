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
#' Each CARD object has a number of slots which store information. Key slots to access
#' are listed below.
#'
#' @slot sc_eset The filtered scRNA-seq data along with meta data stored in the format of SingleCellExperiment.
#' @slot spatial_countMat The filtered spatial count data.
#' @slot spatial_location The weights for combining p-values from multiple kernels.
#' @slot Proportion_CARD The estimated cell type proportion by CARD with each row is a spatial location and each column is a cell type.
#' @slot project The name of the project, default is deconvolution.
#' @slot info_parameters The paramters that are used in model fitting.
#' @slot algorithm_matrix The intermediate matrices that are used in the model fitting step.
#' @slot refined_prop The refined cell type proportion matrix estimated by CARD for the newly grided spatial locations. The number of initial grids are defined by the user.
#' @slot refined_expression The refined predicted expression matrix (normalized) estimated by CARD for the newly grided spatial locations. The number of initial grids are defined by the user.
#'
setClass("CARD", 
	slots = list(
	sc_eset = "ANY",
	spatial_countMat = "ANY",
	spatial_location = "data.frame",
	Proportion_CARD = "matrix",
	project = "character",
	info_parameters = "list",
	algorithm_matrix = "list",
	refined_prop = "matrix",
	refined_expression = "matrix")
	)

#' Quality control of scRNA-seq count data
#'
#' @param counts_in Raw scRNAseq count data, each column is a cell and each row is a gene.
#' @param metaData data frame, metaData with "ct.varname" specify the cell type annotation information and "sample.varname" specify the sample information
#' @param ct.varname character, the name of the column in metaData that specifies the cell type annotation information
#' @param ct.select vector of cell type names that you are interested in to deconvolute, default as NULL. If NULL, then use all cell types provided by single cell dataset;
#' @param sample.varname character,the name of the column in metaData that specifies the sample information. If NULL, we just use the whole as one sample.
#' @param min.cells numeric, we filtered out the non-expressed cells.
#' @param min.genes numeric we filtered out the non-expressed genes
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @return Return the filtered scRNA-seq data and meta data stored in a S4 class (SingleCellExperiment)
#'
#' @export
#'
sc_QC <- function(counts_in,metaData,ct.varname,ct.select,sample.varname = NULL, min.cells = 0,min.genes = 0){
# Filter based on min.features
    coldf = metaData
    counts = counts_in
    if (min.genes >= 0) {
        nfeatures <- colSums(x = counts )
        counts <- counts[, which(x = nfeatures > min.genes)]
        coldf <- coldf[which(x = nfeatures > min.genes),]
    }
    # filter genes on the number of cells expressing
    if (min.cells >= 0) {
        num.cells <- rowSums(x = counts > 0)
        counts <- counts[which(x = num.cells > min.cells), ]
    }
    fdata = as.data.frame(rownames(counts))
    rownames(fdata) = rownames(counts)
    keepCell = as.character(coldf[,ct.varname]) %in% ct.select
    counts = counts[,keepCell]
    coldf = coldf[keepCell,]
    keepGene = rowSums(counts) > 0
    fdata = as.data.frame(fdata[keepGene,])
    counts = counts[keepGene,]
    sce <- SingleCellExperiment(list(counts=counts),
    colData=as.data.frame(coldf),
    rowData=as.data.frame(fdata))
    return(sce)
}

#' Create the CARD object
#'
#' @param sc_count Raw scRNA-seq count data, each column is a cell and each row is a gene.
#' @param sc_meta data frame, with each row representing the cell type and/or sample information of a specific cell. The row names of this data frame should match exactly with the column names of the sc_count data
#' @param spatial_count Raw spatial resolved transcriptomics data, each column is a spatial location, and each row is a gene. 
#' @param spatial_location data frame, with two columns representing the x and y coordinates of the spatial location. The rownames of this data frame should match eaxctly with the columns of the spatial_count.
#' @param ct.varname character, the name of the column in metaData that specifies the cell type annotation information
#' @param ct.select vector of cell type names that you are interested in to deconvolute, default as NULL. If NULL, then use all cell types provided by single cell dataset;
#' @param sample.varname character,the name of the column in metaData that specifies the sample information. If NULL, we just use the whole as one sample.
#' @param minCountGene Minimum counts for each gene 
#' @param minCountSpot Minimum counts for each spatial location
#'
#' @importFrom SummarizedExperiment assays
#' @import methods
#' @return Returns CARD object with filtered spatial count and single cell RNA-seq dataset.
#'
#' @export
#'
createCARDObject <- function(sc_count,sc_meta,spatial_count,spatial_location,ct.varname,ct.select,sample.varname,minCountGene = 100,minCountSpot =5){  

#### QC on scRNASeq dataset
cat(paste0("## QC on scRNASeq dataset! ...\n"))
if(is(sc_count,"matrix")){
	sc_countMat  <- as(as.matrix(sc_count), "sparseMatrix")
	}else if(is(sc_count,"vector")){
		sc_countMat  <- as(t(as.matrix(sc_count)), "sparseMatrix")
		}else if(is(sc_count,"sparseMatrix")){
			sc_countMat  <- sc_count
			}else{
				stop("scRNASeq counts has to be of following forms: vector,matrix or sparseMatrix")
}
if (missing(x = sc_countMat)) {
	stop("Please provide scRNASeq count data")
	} else if (is.null(sample.varname) || missing(sample.varname)) {
		sample.varname = "Sample"
		sc_meta = as.data.frame(sc_meta)
		sc_meta$sampleID = "Sample"
		} else if (any(rownames(x = sc_countMat) == '')) {
			stop("Feature names of sc_count matrix cannot be empty", call. = FALSE)
			} else if(sum(rownames(sc_meta) == colnames(sc_countMat)) != ncol(sc_countMat)){
				stop("Cell name in scRNAseq count data does not match with the rownames of metaData")
				} else if(ncol(sc_countMat)!=nrow(sc_meta)){
					stop("The number of cells in scRNA-seq counts and sc_meta should be consistent! (sc_count -- p x c; sc_meta -- c x 2)")
}
if (is.null(ct.varname)){
	stop("Please provide the column name indicating the cell type information in the meta data of scRNA-seq")
	}else if (is.null(ct.select)){
		cat(paste0("No cell types selected, we will use all the cell types in the scRNA-seq data\n"))
		ct.select <- unique(sc_meta[, ct.varname])
	}
ct.select <- as.character(ct.select[!is.na(ct.select)])
sc_eset = sc_QC(sc_countMat,sc_meta,ct.varname,ct.select,sample.varname)
#### Check the spatial count dataset
#### QC on spatial dataset
cat(paste0("## QC on spatially-resolved dataset! ...\n"))

if(is(spatial_count,"matrix")){
	spatial_countMat  <- as(as.matrix(spatial_count), "sparseMatrix")
	}else if(is(spatial_count,"vector")){
		spatial_countMat  <- as(t(as.matrix(spatial_count)), "sparseMatrix")
		}else if(is(spatial_count,"sparseMatrix")){
			spatial_countMat <- spatial_count
			}else{
				stop("spatial resolved transcriptomic counts has to be of following forms: vector,matrix or sparseMatrix")
}
if (any(rownames(x = spatial_countMat) == '')) {
			stop("Gene names of spatial count matrix cannot be empty", call. = FALSE)
}
commonGene = intersect(rownames(spatial_countMat),rownames(assays(sc_eset)$counts))
if (length(commonGene) == 0) {
			stop("There are no common gene names in spatial count data and single cell RNAseq count data", call. = FALSE)
}
if(is.null(spatial_location)){
	stop("Please provide the matched spatial location data frame")
}
if(ncol(spatial_countMat)!=nrow(spatial_location)){
	stop("The number of spatial locations in spatial_count and spatial_location should be consistent! (spatial_count -- p x n; spatial_location -- n x 2)")
	}# end fi
## check data order should consistent
if(!identical(colnames(spatial_countMat), rownames(spatial_location))){
	stop("The column names of spatial_count and row names of spatial_location should be should be matched each other! (spatial_count -- p x n; spatial_location -- n x 2)")
	}# end fi
#### QC on spatial dataset
spatial_countMat = spatial_countMat[rowSums(spatial_countMat > 0) > minCountSpot,]
spatial_countMat = spatial_countMat[,(colSums(spatial_countMat) >= minCountGene & colSums(spatial_countMat) <= 1e6)]
spatial_location = spatial_location[rownames(spatial_location) %in% colnames(spatial_countMat),]
spatial_location = spatial_location[match(colnames(spatial_countMat),rownames(spatial_location)),]

object <- new(
		Class = "CARD",
		sc_eset = sc_eset,
		spatial_countMat = spatial_countMat,
		spatial_location = spatial_location,
		project = "Deconvolution",
		info_parameters = list(ct.varname = ct.varname,ct.select = ct.select,sample.varname = sample.varname)
		)
return(object)
}

#' Each CARDfree object has a number of slots which store information. Key slots to access
#' are listed below.
#'
#' @slot spatial_countMat The filtered spatial count data.
#' @slot spatial_location The weights for combining p-values from multiple kernels.
#' @slot Proportion_CARD The estimated cell type proportion by CARD with each row is a spatial location and each column is a cell type.
#' @slot estimated_refMatrix The estimated reference matrix by CARDfree with each row represents a gene and each column represents a cell type cluster. 
#' @slot project The name of the project, default is deconvolution.
#' @slot markerList The nlist of cell type specific markers, with each element represents the vector of cell type specific markers
#' @slot info_parameters The paramters that are used in model fitting.
#' @slot algorithm_matrix The intermediate matrices that are used in the model fitting step.
#' @slot refined_prop The refined cell type proportion matrix estimated by CARD for the newly grided spatial locations. The number of initial grids are defined by the user.
#' @slot refined_expression The refined predicted expression matrix (normalized) estimated by CARD for the newly grided spatial locations. The number of initial grids are defined by the user.
#'
setClass("CARDfree", 
	slots = list(
	spatial_countMat = "ANY",
	spatial_location = "data.frame",
	Proportion_CARD = "matrix",
	estimated_refMatrix = "matrix",
	project = "character",
	markerList = "list",
	info_parameters = "list",
	algorithm_matrix = "list",
	refined_prop = "matrix",
	refined_expression = "matrix")
	)
#' Create the CARD object
#'
#' @param markerList a list of marker genes, with each element of the list being the vector of cell type specific marker genes
#' @param spatial_count Raw spatial resolved transcriptomics data, each column is a spatial location, and each row is a gene. 
#' @param spatial_location data frame, with two columns representing the x and y coordinates of the spatial location. The rownames of this data frame should match eaxctly with the columns of the spatial_count.
#' @param minCountGene Minimum counts for each gene 
#' @param minCountSpot Minimum counts for each spatial location
#'
#' @import methods
#' @return Returns CARDfree object with filtered spatial count and marker gene list.
#'
#' @export
#'
createCARDfreeObject <- function(markerList,spatial_count,spatial_location,minCountGene = 100,minCountSpot =5){  

if(is(spatial_count,"matrix")){
	spatial_countMat  <- as(as.matrix(spatial_count), "sparseMatrix")
	}else if(is(spatial_count,"vector")){
		spatial_countMat  <- as(t(as.matrix(spatial_count)), "sparseMatrix")
		}else if(is(spatial_count,"sparseMatrix")){
			spatial_countMat <- spatial_count
			}else{
				stop("spatial resolved transcriptomic counts has to be of following forms: vector,matrix or sparseMatrix")
}# end fi
if (any(rownames(x = spatial_countMat) == '')) {
			stop("Gene names of spatial count matrix cannot be empty", call. = FALSE)
}# end fi
if(is.null(spatial_location)){
	stop("Please provide the matched spatial location data frame")
}# end fi
if(ncol(spatial_countMat)!=nrow(spatial_location)){
	stop("The number of spatial locations in spatial_count and spatial_location should be consistent! (spatial_count -- p x n; spatial_location -- n x 2)")
}# end fi
## check data order should consistent
if(!identical(colnames(spatial_countMat), rownames(spatial_location))){
	stop("The column names of spatial_count and row names of spatial_location should be should be matched each other! (spatial_count -- p x n; spatial_location -- n x 2)")
}# end fi
#### QC on spatial dataset
spatial_countMat = spatial_countMat[rowSums(spatial_countMat > 0) > minCountSpot,]
spatial_countMat = spatial_countMat[,(colSums(spatial_countMat) >= minCountGene & colSums(spatial_countMat) <= 1e6)]
spatial_location = spatial_location[rownames(spatial_location) %in% colnames(spatial_countMat),]
spatial_location = spatial_location[match(colnames(spatial_countMat),rownames(spatial_location)),]
#### check marker gene list
if (missing(x = markerList)) {
	stop("Please provide the marker list for CARDfree")
} 		

object <- new(
		Class = "CARDfree",
		spatial_countMat = spatial_countMat,
		spatial_location = spatial_location,
		project = "Deconvolution (reference-free)",
		markerList = markerList,
		info_parameters = list()
		)
return(object)
}


