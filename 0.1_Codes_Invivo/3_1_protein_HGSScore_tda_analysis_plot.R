#################### UMAP NetworkPlot ####################
# Author: Huitian (Yolanda) Diao
# June 7th, 2019

#################### Libraries ####################
library(dplyr)
library(tidyverse)
library(ggrepel)

#################### Self defined functions ####################

### Find genes of the target complex
# Reference file: generated from GO term annotation
findGoGenes <- function(term_name){
  # Read reference
  ref_go_file <- "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/2_GO_terms/CRM_complexes_count.csv"
  ref_go_tab <- read_csv(ref_go_file)
  ref_go_tab$gene_name <- toupper(ref_go_tab$gene_name)
  
  slt_cols <- c("gene_name",term_name)
  colNums <- match(slt_cols,names(ref_go_tab))
  ref_go_tab <- ref_go_tab %>% 
    select(colNums)
  ref_go_tab <- ref_go_tab %>%
    rename_at(vars(colnames(ref_go_tab)), ~ c("gene_name", "term_name")) %>%
    filter(term_name == "Yes")
  return(ref_go_tab)
}

### Annotate genes based on the complex they belong to
annoComplex <- function(all_genes, term_vector, term_newname, col_vector, setAlpha, outRefName){
  #term_vector <- c("SWI/SNF complex", "NuRD complex", "Set1C/COMPASS complex", "histone H3 acetylation", "ESC/E(Z) complex")
  #term_newname <- c("SWISNF", "NuRD", "Set1cCompass", "HistoneH3Acetylation", "ESCEZComplex")
  #col_vector <- c("firebrick1","deepskyblue", "limegreen", "orange", "purple")
  #all_genes <- c("PHF17")
  #setAlpha <- 1
  
  # Read reference
  ref_go_file <- "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/2_GO_terms/CRM_complexes_count.csv"
  ref_go_tab <- read_csv(ref_go_file)
  ref_go_tab$gene_name <- toupper(ref_go_tab$gene_name)
  
  gene_tb <- tibble(gene_name = all_genes)
  for (i in c(1: length(term_vector))){
    term <- term_vector[i]
    term_tb <- findGoGenes(term) 
    term_tb <- term_tb %>% rename(!!term_newname[i] := term_name)
    gene_tb <- gene_tb %>% left_join(term_tb)
  }
  
  gene_tsps <- numeric(length(all_genes))
  border <- numeric(length(all_genes))
  complex <- character(length(all_genes))
  for (i in c(1: nrow(gene_tb))){
    row_i <- unlist(slice(gene_tb, i), use.names=FALSE)[0:length(term_vector)+1]
    col_true <- 0
    for (j in c(1: length(term_vector))){
      if (!(row_i[1] %in% ref_go_tab$gene_name)){
        gene_tsps[i] <- setAlpha*0.85
        complex[i] <- "NoAnno"
        border[i] <- 1
      }
      else {
        if (!is.na(row_i[j+1])){
          col_true <- col_true + 1
          if (col_true == 1){
            gene_tsps[i] <- setAlpha
            complex[i] <- term_newname[j]
            border[i] <- 0
          }
          else{
            gene_tsps[i] <- setAlpha*0.85
            complex[i] <- "multi"
            border[i] <- 0
            msg <- paste("Multiple matches:", unlist(slice(gene_tb, i), use.names=FALSE)[1])
            print(msg)
          }
        }
        if (col_true == 0){
          gene_tsps[i] <- setAlpha*0.85
          complex[i] <- "other"
          border[i] <- 0
        }
      }
    }
  }
  return(list(gene_tsps, border, complex))
}

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

######################################## Main ########################################
scPlot <- function(coorFile, edgeFile, edgeFilter, termAnno, termNames, termCols, outName){
  #coorFile <- "HGSCore_only-CRF_sq_dist_nn12-mdist0.5_coor.csv"
  #edgeFile <- "HGSCore_only-CRF_sq_dist_nn12-mdist0.5_nn_dupr.csv"
  #edgeFilter <- 0.7
  #termAnno <- c("SWI/SNF complex", "NuRD complex", "Set1C/COMPASS complex", "histone H3 acetylation", 
  #                "ESC/E(Z) complex", "PcG protein complex", "PRC1 complex")
  #termNames <- c("SWISNF", "NuRD", "Set1cCompass", "HistoneH3Acetylation", "ESCEZComplex", "PcGProteinComplex", "PRC1Complex")
  #termCols <- c("firebrick1","deepskyblue", "limegreen", "orange", "purple", "aquamarine", "yellow")
  #outName <- "networkplot.new.pdf"
  
  # Set transparency
  tspEdge <- 0.5
  tspDot <- 1
  outName2 <- gsub(".pdf", "_ref.csv", outName)
  
  # Read coordinates and edge files
  coor.tab <- read_csv(coorFile)
  edge.tab <- read_csv(edgeFile) %>% filter(val > edgeFilter)
  
  # Annotate genes by complexes
  annos <- annoComplex(coor.tab$gene_name, termAnno, termNames, termCols, tspDot, outName2)
  coor.tab$complexNames <- factor(unlist(annos[3]), levels=c(termNames, "multi", "other", "NoAnno"))
  coor.tab$tsp <- unlist(annos[1])
  coor.tab$bd <- unlist(annos[2])

  # Annotate edge tibble, add coordiate and transparency
  edge.tab <- edgeCoor(edge.tab, coor.tab, tspEdge)
  
  ### sc.plot
  # Color setup
  if ("multi" %in% coor.tab$complexNames){
    termCols <- c(termCols, "black")
  }
  if ("other" %in% coor.tab$complexNames){
    termCols <- c(termCols, "grey")
  }
  if ("NoAnno" %in% coor.tab$complexNames){
    termCols <- c(termCols, "white")
  }
  
  # Plot
  sc.plot <- ggplot() +
    geom_segment(data=edge.tab, aes(x = idx1x, y = idx1y, xend = idx2x, yend = idx2y, alpha=tsp)) +
    geom_point(data=coor.tab, aes(x=x, y=y, color=complexNames, size=5, alpha=tsp, stroke=bd)) + 
    scale_color_manual(values=termCols) +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill = "gray55"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  sc.plot
  ggsave(outName, sc.plot, device = 'pdf', width = 25, height = 20, units = "cm")
}

#################### Plotting ####################
wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/3_CRM_protein_HGSScore/1_umap"
setwd(wk.dir)

### Generate all annotated reference plot
if(FALSE){
  coor.file <- "HGSCore_only-CRF_sq_dist_nn12-mdist0.5_coor.csv"
  coor.tab <- read_csv(coor.file)
  sc.plot <- ggplot(data=coor.tab, aes(x=x, y=y, size=8, stroke=0, label = coor.tab$gene_name)) +
    geom_point(alpha=0.3) + 
    geom_text_repel()
  ggsave("Anno_scplot.pdf", sc.plot, device = 'pdf', width = 40, height = 35, units = "cm")
}

### Genreate annotated plot for each complex
wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/3_CRM_protein_HGSScore/1_umap/HGSCore_only-CRF_sq_dist_nn10-mdist0.4"
setwd(wk.dir)

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


coor.file <- "HGSCore_only-CRF_sq_dist_nn12-mdist0.5_coor.csv"
edge.file <- "HGSCore_only-CRF_sq_dist_nn12-mdist0.5_nn_dupr.csv"
edge.filter <- 0.7
terms.anno <- c("SWI/SNF complex", "NuRD complex", "Set1C/COMPASS complex", "histone H3 acetylation", 
                "ESC/E(Z) complex", "PcG protein complex")
terms.names <- c("SWISNF", "NuRD", "Set1cCompass", "HistoneH3Acetylation", "ESCEZComplex", "PcGProteinComplex")
terms.cols <- c("firebrick1","deepskyblue", "limegreen", "orange", "purple", "aquamarine")
out.name <- "networkplot.new.pdf"

scPlot(coor.file, edge.file, edge.filter, terms.anno, terms.names, terms.cols, out.name)





