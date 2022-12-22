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
#' Visualize the spatial distribution of cell type proportion
#'
#' @param proportion Data frame, cell type proportion estimated by CARD in either original resolution or enhanced resolution.
#' @param spatial_location Data frame, spatial location information.
#' @param ct.visualize Vector of selected cell type names that are interested to visualize
#' @param colors Vector of color names that you want to use, if NULL, we will use the default color scale c("lightblue","lightyellow","red")
#' @param NumCols Numeric, number of columns in the figure panel, it depends on the number of cell types you want to visualize.
#'
#' @import ggplot2 
#' @importFrom reshape2 melt
#' @return Returns a ggplot2 figure. 
#'
#' @export
#'

CARD.visualize.prop <- function(proportion,spatial_location,ct.visualize = ct.visualize,colors = c("lightblue","lightyellow","red"),NumCols){
if(is.null(colors)){
	colors = c("lightblue","lightyellow","red")
}else{
	colors = colors
}
res_CARD = as.data.frame(proportion)
res_CARD = res_CARD[,order(colnames(res_CARD))]
location = as.data.frame(spatial_location)
if(sum(rownames(res_CARD)==rownames(location))!= nrow(res_CARD)){
   stop("The rownames of proportion data does not match with the rownames of spatial location data")
}
ct.select = ct.visualize
res_CARD = res_CARD[,ct.select]
res_CARD_scale = as.data.frame(apply(res_CARD,2,function(x){
    (x - min(x)) / (max(x) - min(x))
} ))
res_CARD_scale$x = as.numeric(location$x)
res_CARD_scale$y = as.numeric(location$y)
mData = melt(res_CARD_scale,id.vars = c("x","y"))
colnames(mData)[3] <- "Cell_Type"
b = c(0,1)
p = suppressMessages(ggplot(mData, aes(x, y)) + 
geom_point(aes(colour = value),size = 3.0) +
scale_color_gradientn(colours = colors) + 
#scale_color_viridis_c(option = 2)+
scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+ 
facet_wrap(~Cell_Type,ncol = NumCols)+ 
coord_fixed()+
theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    #legend.position=c(0.14,0.76),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
    axis.text =element_blank(),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 14,face="bold"),
    legend.text=element_text(size = 11),
    strip.text = element_text(size = 12,face="bold"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')))
return(p)
}

#' Visualize the spatial distribution of two cell type proportions on the same plot 
#'
#' @param proportion Data frame, cell type proportion estimated by CARD in either original resolution or enhanced resolution.
#' @param spatial_location Data frame, spatial location information.
#' @param ct2.visualize Vector of selected two cell type names that are interested to visualize, here we only focus on two cell types
#' @param colors list of color names that you want to use for each cell type, if NULL, we will use the default color scale list
#' list(c("lightblue","lightyellow","red"),c("lightblue","lightyellow","black")
#'
#' @import ggplot2 
#' @importFrom reshape2 melt
#' @return Returns a ggplot2 figure. 
#'
#' @export
#'
CARD.visualize.prop.2CT <- function(proportion,spatial_location,ct2.visualize = ct2.visualize,colors = NULL){
if(is.null(colors)){
	colors = list(c("lightblue","lightyellow","red"),
c("lightblue","lightyellow","black"))
}else{
	colors = colors
}
res_CARD = as.data.frame(proportion)
res_CARD = res_CARD[,order(colnames(res_CARD))]
location = as.data.frame(spatial_location)
if(sum(rownames(res_CARD)==rownames(location))!= nrow(res_CARD)){
   stop("The rownames of proportion data does not match with the rownames of spatial location data")
}
ct.select = ct2.visualize
res_CARD = res_CARD[,ct.select]
res_CARD_scale = as.data.frame(apply(res_CARD,2,function(x){
    (x - min(x)) / (max(x) - min(x))
} ))
res_CARD_scale$x = as.numeric(location$x)
res_CARD_scale$y = as.numeric(location$y)
mData = melt(res_CARD_scale,id.vars = c("x","y"))
colnames(mData)[3] <- "Cell_Type"
b = c(0,1)
p = suppressMessages(ggplot() + 
geom_point(data=mData[mData$Cell_Type == ct2.visualize[1],], aes(x=x, y=y, color=value), shape=21, size=5) +
scale_color_gradientn(colours = colors[[1]]) +
geom_point(data=mData[mData$Cell_Type == ct2.visualize[2],], aes(x=x, y=y, fill = value), color = "white",shape=22, size=2) +
scale_fill_gradientn(colours = colors[[2]]) +
#scale_color_viridis_c(option = 2)+
scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+ 
coord_fixed()+
theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    #legend.position=c(0.14,0.76),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
    axis.text =element_blank(),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 14,face="bold"),
    legend.text=element_text(size = 11),
    strip.text = element_text(size = 12,face="bold"),
    #legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')))
return(p)
}

#' Visualize the spatial distribution of cell type proportion in a geom scatterpie plot
#'
#' @param proportion Data frame, cell type proportion estimated by CARD in either original resolution or enhanced resolution.
#' @param spatial_location Data frame, spatial location information.
#' @param colors Vector of color names that you want to use, if NULL, we will use the color palette "Spectral" from RColorBrewer package. 
#'
#' @import ggplot2 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scatterpie geom_scatterpie
#' @importFrom grDevices colorRampPalette
#' @importFrom gtools mixedsort
#' @return Returns a ggplot2 figure. 
#'
#' @export
#'

CARD.visualize.pie <- function(proportion,spatial_location,colors = NULL){
res_CARD = as.data.frame(proportion)
res_CARD = res_CARD[,mixedsort(colnames(res_CARD))]
location = as.data.frame(spatial_location)
if(sum(rownames(res_CARD)==rownames(location))!= nrow(res_CARD)){
   stop("The rownames of proportion data does not match with the rownames of spatial location data")
}
colorCandidate = c("#1e77b4","#ff7d0b","#ceaaa3","#2c9f2c","#babc22","#d52828","#9267bc",
  "#8b544c","#e277c1","#d42728","#adc6e8","#97df89","#fe9795","#4381bd","#f2941f","#5aa43a","#cc4d2e","#9f83c8","#91675a",
  "#da8ec8","#929292","#c3c237","#b4e0ea","#bacceb","#f7c685",
  "#dcf0d0","#f4a99f","#c8bad8",
  "#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3",
 "#D1D1D1", "#6D1A9C", "#15821E", "#3A84E6", "#997273",
 "#787878", "#DB4C6C", "#9E7A7A", "#554236", "#AF5F3C",
 "#93796C", "#F9BD3F", "#DAB370", "#877F6C", "#268785")
if(is.null(colors)){
	#colors = brewer.pal(11, "Spectral")
	if(ncol(res_CARD) > length(colorCandidate)){
	colors = colorRampPalette(colorCandidate)(ncol(res_CARD))
	}else{
		colors = colorCandidate[sample(1:length(colorCandidate),ncol(res_CARD))]
	}
}else{
	colors = colors
}
data = cbind(res_CARD,location)
ct.select = colnames(res_CARD)
p = suppressMessages(ggplot() + geom_scatterpie(aes(x=x, y=y,r = 0.52),data=data,
                                cols=ct.select,color=NA) + coord_fixed(ratio = 1) + 
                                scale_fill_manual(values =  colors)+
                                theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
    axis.text =element_blank(),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 15),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm'),
    strip.text = element_text(size = 16,face="bold"),
    legend.position="bottom")+
                                guides(fill=guide_legend(title="Cell Type")))
return(p)
}

#' Visualize the cell type proportion correlation
#'
#' @param proportion Data frame, cell type proportion estimated by CARD in either original resolution or enhanced resolution.
#' @param colors Vector of color names that you want to use, if NULL, we will use the default color scale c("#91a28c","white","#8f2c37")
#'
#' @import ggcorrplot 
#' @importFrom stats cor
#' @return Returns a ggcorrplot figure. 
#'
#' @export
#'
CARD.visualize.Cor <- function(proportion,colors = colors){
proportion = proportion[,order(colnames(proportion))]
cor_CARD = cor(as.matrix(proportion))
if(is.null(colors)){
	colors = c("#91a28c","white","#8f2c37")
}else{
	colors = colors
}

p = suppressMessages(ggcorrplot(cor_CARD, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = colors)+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) + 
 coord_fixed()+
 ggtitle("Correlation")+theme(plot.title = element_text(size=22,face="bold")))
return(p)
}

#' Visualize the spatial distribution of cell type proportion
#'
#' @param spatial_expression Data frame, spatial gene expression in either original resolution or enhanced resolution.
#' @param spatial_location Data frame, spatial location information.
#' @param gene.visualize Vector of selected gene names that are interested to visualize
#' @param colors Vector of color names that you want to use, if NULL, we will use the default color scale in virdis palette
#' @param NumCols Numeric, number of columns in the figure panel, it depends on the number of cell types you want to visualize.
#'
#' @import ggplot2 
#' @return Returns a ggplot2 figure. 
#'
#' @export
#'
CARD.visualize.gene <- function(spatial_expression,spatial_location,gene.visualize,colors = colors,NumCols){
expression = as.data.frame(as.matrix(spatial_expression))
expression = sweep(expression,2,colSums(expression),"/")
location = as.data.frame(spatial_location)
if(sum(colnames(expression)==rownames(location))!= nrow(location)){
   stop("The colnames of expression data does not match with the rownames of spatial location data")
}
gene.select = gene.visualize
if(sum(toupper(gene.select) %in% toupper(rownames(spatial_expression))) != length(gene.select)){
stop("There existing selected genes that are not in the expression data!")
}
Data = NULL
for(i in 1:length(gene.select)){
#### load spatial dataset
gene = gene.select[i]
ind = which(toupper(rownames(expression)) == toupper(gene))
df = as.numeric(expression[ind,])
names(df) = colnames(expression)
df = (df - min(df)) / (max(df) - min(df))
d = data.frame(value = df,
x=as.numeric(location$x),
y = as.numeric(location$y))
d$gene = gene
Data = rbind(Data,d)
}
Data$gene = factor(Data$gene,levels = gene.select)
p = suppressMessages(ggplot(Data, aes(x, y)) + 
geom_point(aes(color = value),size = 1.5,shape = 15,position = position_dodge2(width = 0.05))+
scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+
coord_equal()+
facet_wrap(~gene,ncol = NumCols)+
theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    legend.position="bottom",
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
    axis.text =element_blank(),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 15,face="bold"),
    legend.text=element_text(size = 14),
    strip.text = element_text(size = 18,face="bold"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(1.0, 'cm'))+
guides(color = guide_colourbar(title = "Expression")))
if(is.null(colors)){
p <- p + scale_color_viridis_c(labels = c("0","0.25","0.5","0.75","1.0"))
}else{
p <- p + scale_color_gradientn(colours = colors)
}
return(p)
}


