########## Heatmap-plot ##########
# Author: Huitian (Yolanda) Diao
# Dec 3rd, 2018
# For plottinng Heatmaps

########## Libraries ##########
library(pheatmap)
library(RColorBrewer)

########## Self-defined functions ##########
fnameNoPathPdf <- function(namex){
  namex.vec <- unlist(strsplit(namex, "/"))
  namex <- tail(namex.vec, n=1)
  namex <- gsub(".csv", ".pdf", namex)
  return(namex)
}

fnameNoPath <- function(namex){
  namex.vec <- unlist(strsplit(namex, "/"))
  namex <- tail(namex.vec, n=1)
  namex <- gsub(".csv", "", namex)
  return(namex)
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

setRnames <- function(tabx){
  row.names(tabx) <- unlist(tabx[,1])
  tabx[,1] <- NULL
  return(tabx)
}

plot_heatmap_Zscore <- function(infile, col.pal, pdf.l, pdf.h, max){
  outname <- fnameNoPathPdf(infile)
  intab <- read.csv(infile)
  intab <- setRnames(intab)
  myBreaks <- seq(-max, max, by=.01)
  myColor <- colorRampPalette(col.pal)(length(myBreaks))
  
  heat.map <- pheatmap(mat = intab, fontsize = 6, 
                       cluster_cols = TRUE, cluster_rows = TRUE,
                       color=myColor, breaks=myBreaks)
  # clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
  # ,color=myColor, breaks=myBreaks
  #                       cellwidth = 9, cellheight = 9, main = infile,
  #cluster_cols = FALSE, cluster_rows = FALSE
  # show_rownames = TRUE, show_colnames = TRUE, fontsize = 6
  save_pheatmap_pdf(heat.map, outname, pdf.l, pdf.h)
}

########## Main ##########

###----- Plot Z scores
if (FALSE){
  wk.dir <- "/Users/yolandatiao/Desktop/CRF_Screen/2_0_t-test_by_gene/compiled_Amt"
  setwd(wk.dir)
  
  amt.z <- "Amt_normbycontrolZP_t-test.by.geneavg_z-score.csv"
  z.cols <- c("midnightblue", "royalblue1", "white", "orangered", "darkred")
  plot_heatmap_Zscore(amt.z, z.cols, 4, 25, 6)
  
  amt.z <- "Amt_normbycontrolZP_t-test.by.geneavg_z-score_10U.csv"
  z.cols <- c("midnightblue", "royalblue1", "white", "orangered", "darkred")
  plot_heatmap_Zscore(amt.z, z.cols, 3, 25, 6)

  amt.z <- "Amt_normbycontrolZP_t-test.by.geneavg_z-score_100U.csv"
  z.cols <- c("midnightblue", "royalblue1", "white", "orangered", "darkred")
  plot_heatmap_Zscore(amt.z, z.cols, 3, 25, 6)
  
  amt.z <- "Amt_normbycontrolZP_t-test.by.geneavg_z-score_GeoMean.csv"
  z.cols <- c("midnightblue", "royalblue1", "white", "orangered", "darkred")
  plot_heatmap_Zscore(amt.z, z.cols, 3, 25, 6)
  
  amt.z <- "Amt_normbycontrolZP_t-test.by.geneavg_z-score_Percentage.csv"
  z.cols <- c("midnightblue", "royalblue1", "white", "orangered", "darkred")
  plot_heatmap_Zscore(amt.z, z.cols, 3, 25, 6)
  
  amt.z <- "Amt_normbycontrolZP_t-test.by.geneavg_z-score_CD127-CD25_10U.csv"
  z.cols <- c("midnightblue", "royalblue1", "white", "orangered", "darkred")
  plot_heatmap_Zscore(amt.z, z.cols, 3, 25, 6)

  amt.z <- "Amt_normbycontrolZP_t-test.by.geneavg_z-score_CD127-CD25_100U.csv"
  z.cols <- c("midnightblue", "royalblue1", "white", "orangered", "darkred")
  plot_heatmap_Zscore(amt.z, z.cols, 3, 25, 6)  
}






