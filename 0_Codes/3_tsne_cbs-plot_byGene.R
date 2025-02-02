########## CBS-plot ##########
# Author: Huitian (Yolanda) Diao
# Dec 3rd, 2018
# For plottinng color-bar scatter plot

########## Libraries ##########
library(ggplot2)
library(ggrepel)

########## Self-defined functions ##########
build_df <- function(tsne.file, z.score.file, p.value.file, key.col){
  #key.col <- "CD127GeoMean_10U"
  tsnf.df <- read.csv(tsne.file)
  z.score.df <- read.csv(z.score.file)
  p.value.df <- read.csv(p.value.file)
  #colnames(z.score.df) <- c("gene_name", "CD103Percentage_100U", "CD103Percentage_10U", "CD127GeoMean_100U", "CD127GeoMean_10U", "CD25GeoMean_100U", "CD25GeoMean_10U", "CD44GeoMean_100U", "CD44GeoMean_10U", "CD62LPercentage_100U", "CD62LPercentage_10U", "CXCR3GeoMean_100U", "CXCR3GeoMean_10U", "CXCR3Percentage_100U", "CXCR3Percentage_10U", "Lag3GeoMean_100U", "Lag3GeoMean_10U", "PD1GeoMean_100U", "PD1GeoMean_10U", "Tim3GeoMean_100U", "Tim3GeoMean_10U", "Tim3Percentage_100U", "Tim3Percentage_10U")
  #colnames(p.value.df) <- c("gene_name", "CD103Percentage_100U", "CD103Percentage_10U", "CD127GeoMean_100U", "CD127GeoMean_10U", "CD25GeoMean_100U", "CD25GeoMean_10U", "CD44GeoMean_100U", "CD44GeoMean_10U", "CD62LPercentage_100U", "CD62LPercentage_10U", "CXCR3GeoMean_100U", "CXCR3GeoMean_10U", "CXCR3Percentage_100U", "CXCR3Percentage_10U", "Lag3GeoMean_100U", "Lag3GeoMean_10U", "PD1GeoMean_100U", "PD1GeoMean_10U", "Tim3GeoMean_100U", "Tim3GeoMean_10U", "Tim3Percentage_100U", "Tim3Percentage_10U")
  
  new.df <- tsnf.df
  new.df$zScore <- unlist(z.score.df[key.col])
  new.df$pValue <- unlist(p.value.df[key.col])
  
  out.name <- paste(key.col, "_tsne-z-p.csv", sep="")
  write.csv(new.df, out.name, row.names = FALSE)
}

slt_anno <- function(df.x) {
  anno.vec <- character(nrow(df.x))
  anno.num <- 0
  for (x in 1:nrow(df.x)){
    if  ( (df.x$pValue[x] < 0.05) & (abs(df.x$zScore[x]) > 1) ){
      anno.vec[x] <- as.character(unlist(df.x$gene_name))[x]
      anno.num <- anno.num + 1
    }
    else {
      anno.vec[x] <- ""
    }
  }
  df.x$gene_name_anno <- anno.vec
  return(df.x)
}

cbs_plot <- function(in.file){
  in.df <- read.csv(in.file)
  in.df$nlog10p <- -log10(in.df$pValue + 0.00001)
  in.df <- slt_anno(in.df)
  
  # Save plot spreadsheet
  out.csv.name <- gsub(".csv", "_cbs.csv", in.file)
  write.csv(in.df, out.csv.name, row.names = FALSE)
  
  nlog10pval <- in.df$nlog10p
  sqrt_Z <- sqrt(abs(in.df$zScore)) * in.df$zScore / abs(in.df$zScore)
  cbs.plot <- ggplot(in.df, aes(x=in.df$t1, y=in.df$t2, color=sqrt_Z)) +
    geom_point(size=5, stroke=0, ) +
    scale_color_gradient2(midpoint=0, limits = c(-2,2), low="dodgerblue3", mid="white", high="firebrick3", space ="Lab") +
    scale_x_continuous(limits = c(-80, 80)) +
    scale_y_continuous(limits = c(-80, 80)) +
    labs(x ="t1", y = "t2") +
    theme(panel.background = element_rect(fill = "gray50"),
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "black"), 
          panel.grid.minor = element_line(size = 0.05, linetype = 'solid',colour = "black"),
          text = element_text(size=10))
  # Save without annotation
  out.name <- gsub(".csv", ".pdf", in.file)
  ggsave(out.name, plot = cbs.plot, device = "pdf", units = "cm", width = 22, height = 20)
  # Save with annotation
  cbs.plot <- cbs.plot +
    geom_text_repel(aes(label = in.df$gene_name_anno), size = 6)
  out.name <- gsub(".csv", "_anno.pdf", in.file)
  ggsave(out.name, plot = cbs.plot, device = "pdf", units = "cm", width = 22, height = 20)
}

########## Main ##########

###----- Build dataframes
#--- perplexity = 5
if (FALSE){
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/2_Amt_compile_cluster/0_bbplot_source"
  setwd(wk.dir)
  
  tsne <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/1_Amt_compile/Amt_normbycontrolZP_t-test.by.geneavg_z-score_tsne_per5.csv"
  z.score <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/1_Amt_compile/Amt_normbycontrolZP_t-test.by.geneavg_z-score.csv"
  p.value <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/1_Amt_compile/Amt_normbycontrolZP_t-test.by.genet-test_p-value.csv"
  
  colnames.vec <- c("CD103Percentage_100U", "CD103Percentage_10U", "CD127GeoMean_100U", "CD127GeoMean_10U", "CD25GeoMean_100U", "CD25GeoMean_10U", "CD44GeoMean_100U", "CD44GeoMean_10U", "CD62LPercentage_100U", "CD62LPercentage_10U", "CXCR3GeoMean_100U", "CXCR3GeoMean_10U", "CXCR3Percentage_100U", "CXCR3Percentage_10U", "Lag3GeoMean_100U", "Lag3GeoMean_10U", "PD1GeoMean_100U", "PD1GeoMean_10U", "Tim3GeoMean_100U", "Tim3GeoMean_10U", "Tim3Percentage_100U", "Tim3Percentage_10U")
  
  for (i in colnames.vec) {
    build_df(tsne, z.score, p.value, i)
  }
}

###----- Plot
#--- Tsne-per5
if (FALSE){
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/2_Amt_compile_cluster/0_bbplot_source"
  setwd(wk.dir)  
  files <- list.files(path = wk.dir, pattern = "tsne-z-p.csv", full.name=FALSE, recursive=FALSE)
  for (file.x in files){
    cbs_plot(file.x)
  }
}