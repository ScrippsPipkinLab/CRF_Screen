#################### UMAP NetworkPlot ####################
# Author: Huitian (Yolanda) Diao
# June 7th, 2019

#################### Libraries ####################
library(dplyr)
library(tidyverse)
library(ggrepel)

#################### Self defined functions ####################

### Find the coordinate of edges from coordinate dataframe
edgeCoor <- function(edgeTab, coorTab, setAlpha){
  idx1x <- numeric(nrow(edgeTab))
  idx1y <- numeric(nrow(edgeTab))
  idx2x <- numeric(nrow(edgeTab))
  idx2y <- numeric(nrow(edgeTab))
  for (i in c(1:nrow(edgeTab))){
    idx1.i <- edgeTab$idx1[i] + 1
    idx2.i <- edgeTab$idx2[i] + 1
    idx1.i.x <- unlist(slice(coorTab, idx1.i), use.names=FALSE)[2]
    idx1.i.y <- unlist(slice(coorTab, idx1.i), use.names=FALSE)[3]
    idx2.i.x <- unlist(slice(coorTab, idx2.i), use.names=FALSE)[2]
    idx2.i.y <- unlist(slice(coorTab, idx2.i), use.names=FALSE)[3]
    idx1x[i] <- as.numeric(idx1.i.x)
    idx1y[i] <- as.numeric(idx1.i.y)
    idx2x[i] <- as.numeric(idx2.i.x)
    idx2y[i] <- as.numeric(idx2.i.y)
  }
  edgeTab$idx1x <- idx1x
  edgeTab$idx1y <- idx1y
  edgeTab$idx2x <- idx2x
  edgeTab$idx2y <- idx2y
  edgeTab$tsp <- rep(setAlpha, nrow(edgeTab))
  return(edgeTab)
}

### Set transparency based on complex name
tspGene <- function(inList, listx, tspSet){
  outList <- numeric(length(inList))
  for (i in c(1:length(inList))){
    strx <- inList[i]
    if (strx %in% listx){
      outList[i] <- tspSet
    }
    else{
      outList[i] <- tspSet*0.65
    }
  }
  return(outList)
}

### Color top and bottom quarter
tbColor <- function(inList){
  orderList <- sort(inList)
  outList <- character(length(inList))
  lq <- orderList[ceiling(length(inList)/4)]
  hq <- orderList[length(inList) - floor(length(inList)/4)]
  for (i in c(1: length(inList))){
    if(inList[i] <= lq){
      outList[i] <- "bottom"
    } else if (inList[i] >= hq){
      outList[i] <- "top"
    } else{
      outList[i] <- "middle"
    }
  }
  return(outList)
}

### Select annonames
sltAnno <- function(gnList, cpxList, sltCpx){
  outList <- character(length(gnList))
  for (i in c(1: length(gnList))){
    if (cpxList[i] %in% sltCpx){
      outList[i] <- gnList[i]
    }
    else {
      outList[i] <- ""
    }
  }
  return(outList)
}

### Convert gene names
cvtGNames <- function(gnList){
  inList <- c("KDM8", "KAT8", "KAT7", "KAT6B", "KAT6A")
  outList <- c("JMJD5", "MYST1", "MYST2", "MYST4", "MYST3")
  newList <- character(length(gnList))
  for (i in c(1: length(gnList))){
    if (gnList[i] %in% inList){
      idxi <- match(gnList[i], inList)
      newList[i] <- outList[idxi]
    }
    else{
      newList[i] <- gnList[i]
    }
  }
  return(newList)
}

######################################## Main ########################################
scPlotAnno <- function(coorFile, edgeFile, edgeFilter, annoFile, highlightTerms, termNames, termCols, outName){
  coorFile <- "HGSCore_only-CRF_sq_dist_nn12-mdist0.5_coor.csv"
  edgeFile <- "HGSCore_only-CRF_sq_dist_nn12-mdist0.5_nn_dupr.csv"
  annoFile <- "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/3_CRM_protein_HGSScore/1_umap/HGSCore_only-CRF_anno.csv"
  edgeFilter <- 0.6
  highlightTerms <- c("BAF")
  termNames <- c("PRC1", "BAF", "SET-MLL", "PRC2", "FACT", 
                 "HDAC", "E2F6_COM-1", "NU4A", "L3MBTL", "ISWI", 
                 "HBO1", "MOZ-MORF", "CAF1", "SAGA", "SEC", "unknown")
  termCols <- c("aquamarine1", "indianred", "forestgreen", "goldenrod1", "gray0",
                "deeppink", "darkorange", "darkturquoise", "darkslateblue","darkseagreen4",
                "chartreuse1", "firebrick4", "burlywood", "purple", "seagreen3", "gray")
  outName <- "networkplot.anno.pdf"
  
  # Set transparency
  tspEdge <- 0.5
  tspDot <- 0.9
  outName2 <- gsub(".pdf", "_ref.csv", outName)
  
  # Read coordinates and edge files
  coor.tab <- read_csv(coorFile)
  edge.tab <- read_csv(edgeFile) %>% filter(val > edgeFilter)
  anno.tab <- read_csv(annoFile)
  
  # Annotate genes by complexes
  coor.tab <- coor.tab %>% left_join(anno.tab, by="gene_name") %>% mutate(complexNames = replace_na(complexNames, "unknown"))
  coor.tab$tsp <- tspGene(coor.tab$complexNames, highlightTerms, tspDot)
  coor.tab$complexNames <- factor(coor.tab$complexNames, levels=c(termNames))
  coor.tab$bd <- rep(0, nrow(coor.tab))
  coor.tab$annonames <- sltAnno(coor.tab$gene_name, coor.tab$complexNames, highlightTerms)
  
  # Annotate edge tibble, add coordiate and transparency
  edge.tab <- edgeCoor(edge.tab, coor.tab, tspEdge)
  
  # Plot
  sc.plot <- ggplot() +
    geom_segment(data=edge.tab, aes(x = idx1x, y = idx1y, xend = idx2x, yend = idx2y, alpha=tsp)) +
    geom_point(data=coor.tab, aes(x=x, y=y, color=complexNames, size=5, alpha=tsp, stroke=coor.tab$bd)) + 
    scale_color_manual(values=termCols) +
    geom_text_repel() +
    scale_x_continuous(limits = c(-2, 7)) +
    geom_text_repel(data=coor.tab, aes(x=x, y=y, label=coor.tab$annonames), segment.alpha=0.5, force=1, max.iter=50000, seed=123) +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())
  sc.plot
  ggsave(outName, sc.plot, device = 'pdf', width = 18, height = 12, units = "cm")
}

scPlotZScoreOverlay <- function(coorFile, edgeFile, edgeFilter, annoFile, highlighTerms, ZScoreFile, outName){
  coorFile <- "HGSCore_only-CRF_sq_dist_nn12-mdist0.5_coor.csv"
  edgeFile <- "HGSCore_only-CRF_sq_dist_nn12-mdist0.5_nn_dupr.csv"
  annoFile <- "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/3_CRM_protein_HGSScore/1_umap/HGSCore_only-CRF_anno.csv"
  ZScoreFile<- "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/5_zscore_div_sqrt_pval/all_z-score_div_sqrt-p_sqrt.csv"
  edgeFilter <- 0.6
  highlightTerms <- c("BAF")
  outName <- "networkplot.anno.pdf"
  
  # Set transparency
  tspEdge <- 0.5
  tspDot <- 0.8
  outName2 <- gsub(".pdf", "_ref.csv", outName)
  termCols <- c("tomato","deepskyblue", "gray20")
  
  # Read coordinates and edge files
  coor.tab <- read_csv(coorFile)
  edge.tab <- read_csv(edgeFile) %>% filter(val > edgeFilter)
  anno.tab <- read_csv(annoFile)
  z.tab <- read_csv(ZScoreFile) %>% mutate(gene_name = toupper(gene_name))
  
  # Mark top and bottom genes
  z.tab$Q4minusQ1_cor <- tbColor(z.tab$Q4minusQ1)
  z.tab$Q3minusOther_cor <- tbColor(z.tab$Q3minusOther)
  z.tab$InputMinusAvg_cor <- tbColor(z.tab$InputMinusAvg)
  
  # Annotate genes by complexes
  coor.tab <- coor.tab %>%
    left_join(anno.tab, by="gene_name")
  coor.tab$gene_name <- cvtGNames(coor.tab$gene_name)
  coor.tab <- coor.tab %>% 
    left_join(z.tab, by = "gene_name") %>% 
    mutate(Q4minusQ1 = replace_na(Q4minusQ1, 0)) %>%
    mutate(Q3minusOther = replace_na(Q3minusOther, 0)) %>%
    mutate(InputMinusAvg = replace_na(InputMinusAvg, 0)) %>%
    mutate(Q4minusQ1_cor = replace_na(Q4minusQ1_cor, "middle")) %>%
    mutate(Q3minusOther_cor = replace_na(Q3minusOther_cor, "middle")) %>%
    mutate(InputMinusAvg_cor = replace_na(InputMinusAvg_cor, "middle"))
  coor.tab$gene_name
  
  coor.tab$tsp <- tspGene(coor.tab$complexNames, highlightTerms, tspDot)
  coor.tab$bd <- rep(0, nrow(coor.tab))
  coor.tab$annonames <- sltAnno(coor.tab$gene_name, coor.tab$complexNames, highlightTerms)
  
  coor.tab$Q4minusQ1_cor <- factor(coor.tab$Q4minusQ1_cor, levels=c("top", "bottom", "middle"))
  coor.tab$Q3minusOther_cor <- factor(coor.tab$Q3minusOther_cor, levels=c("top", "bottom", "middle"))
  coor.tab$InputMinusAvg_cor <- factor(coor.tab$InputMinusAvg_cor, levels=c("top", "bottom", "middle"))
  
  # Annotate edge tibble, add coordiate and transparency
  edge.tab <- edgeCoor(edge.tab, coor.tab, tspEdge)
  
  outName <- "InputMinusAvg_color_noTXT.pdf"
  # Plot
  sc.plot <- ggplot() +
    geom_segment(data=edge.tab, aes(x = idx1x, y = idx1y, xend = idx2x, yend = idx2y, alpha=tsp)) +
    geom_point(data=coor.tab, aes(x=x, y=y, color=InputMinusAvg_cor, size=3, alpha=tsp, stroke=coor.tab$bd)) + 
    scale_color_manual(values=termCols) +
    #scale_x_continuous(limits = c(-2, 7)) +
    scale_alpha(range = c(0.1, 1)) +
    #geom_text_repel(data=coor.tab, aes(x=x, y=y, label=coor.tab$annonames), segment.alpha=0.5, force=1, max.iter=50000, seed=123) +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position='none',
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())
  #sc.plot
  ggsave(outName, sc.plot, device = 'pdf', width = 13, height = 12, units = "cm")
  
  outName <- "Q4minusQ1_color_noTXT.pdf"
  # Plot
  sc.plot <- ggplot() +
    geom_segment(data=edge.tab, aes(x = idx1x, y = idx1y, xend = idx2x, yend = idx2y, alpha=tsp)) +
    geom_point(data=coor.tab, aes(x=x, y=y, color=Q4minusQ1_cor, size=3, alpha=tsp, stroke=coor.tab$bd)) + 
    scale_color_manual(values=termCols) +
    #scale_x_continuous(limits = c(-2, 7)) +
    scale_alpha(range = c(0.1, 1)) +
    #geom_text_repel(data=coor.tab, aes(x=x, y=y, label=coor.tab$annonames), segment.alpha=0.5, force=1, max.iter=50000, seed=123) +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position='none',
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())
  #sc.plot
  ggsave(outName, sc.plot, device = 'pdf', width = 13, height = 12, units = "cm")

  outName <- "Q3minusOther_color_noTXT.pdf"
  # Plot
  sc.plot <- ggplot() +
    geom_segment(data=edge.tab, aes(x = idx1x, y = idx1y, xend = idx2x, yend = idx2y, alpha=tsp)) +
    geom_point(data=coor.tab, aes(x=x, y=y, color=Q3minusOther_cor, size=3, alpha=tsp, stroke=coor.tab$bd)) + 
    scale_color_manual(values=termCols) +
    #scale_x_continuous(limits = c(-2, 7)) +
    scale_alpha(range = c(0.1, 1)) +
    #geom_text_repel(data=coor.tab, aes(x=x, y=y, label=coor.tab$annonames), segment.alpha=0.5, force=1, max.iter=50000, seed=123) +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position='none',
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())
  #sc.plot
  ggsave(outName, sc.plot, device = 'pdf', width = 13, height = 12, units = "cm")   # With txt: 16, 12
  
}

#################### Plotting ####################
wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/3_CRM_protein_HGSScore/1_umap"
setwd(wk.dir)

### Generate all annotated reference plot
if(FALSE){
  coor.file <- "HGSCore_only-CRF_sq_dist_nn10-mdist0.4_coor.csv"
  coor.tab <- read_csv(coor.file)
  sc.plot <- ggplot(data=coor.tab, aes(x=x, y=y, size=8, stroke=0, label = coor.tab$gene_name)) +
    geom_point(alpha=0.3) + 
    geom_text_repel()
  ggsave("Anno_scplot.pdf", sc.plot, device = 'pdf', width = 35, height = 35, units = "cm")
}

### Genreate annotated plot for each complex
wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/3_CRM_protein_HGSScore/1_umap/HGSCore_only-CRF_sq_dist_nn12-mdist0.5"
setwd(wk.dir)

if (FALSE){
  count.file <- "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/2_GO_terms/complex_count_rank_use_true.txt"
  count.tb <- read_tsv(count.file)
  for (i in c(1:nrow(count.tb))){
    complexi <- count.tb$complex[i]
    if (count.tb$count[i] >= 3) {
      coor.file <- "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/3_CRM_protein_HGSScore/1_umap/HGSCore_only-CRF_sq_dist_nn10-mdist0.4/HGSCore_only-CRF_sq_dist_nn10-mdist0.4_coor.csv"
      edge.file <- "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/3_CRM_protein_HGSScore/1_umap/HGSCore_only-CRF_sq_dist_nn10-mdist0.4/HGSCore_only-CRF_sq_dist_nn10-mdist0.4_nn_dupr.csv"
      edge.filter <- 0.7
      terms.anno <- c(complexi)
      terms.names <- c("complex")
      terms.cols <- c("deepskyblue")
      out.name <- paste(complexi, ".pdf", sep="")
      out.name <- gsub("/", "", out.name)
      out.name <- gsub(" ", "_", out.name)
      scPlot(coor.file, edge.file, edge.filter, terms.anno, terms.names, terms.cols, out.name)
    }
  }
}


if(FALSE) {
  coor.file <- "HGSCore_only-CRF_sq_dist_nn10-mdist0.4_coor.csv"
  edge.file <- "HGSCore_only-CRF_sq_dist_nn10-mdist0.4_nn_dupr.csv"
  edge.filter <- 0.7
  terms.anno <- c("SWI/SNF complex", "NuRD complex", "Set1C/COMPASS complex", "MOZ/MORF histone acetyltransferase complex", 
                  "ESC/E(Z) complex", "PcG protein complex", "NURF complex","ELL-EAF complex",
                  "thiol-dependent ubiquitinyl hydrolase activity")
  terms.names <- c("SWISNF", "NuRD", "Set1cCompass", "MOZMORFHistoneAcetyltransferaseComplex", 
                   "ESCEZComplex", "PcGProteinComplex", "NURFComplex", "ELL-EAFComplex",
                   "Thiol-dependentUbiquitinylHydrolaseActivity")
  terms.cols <- c("firebrick1","deepskyblue", "limegreen", "orange", 
                  "purple", "aquamarine", "royalblue4", "yellow",
                  "palevioletred1")
  out.name <- "networkplot.label.pdf"
  
  scPlot(coor.file, edge.file, edge.filter, terms.anno, terms.names, terms.cols, out.name)
  
}




