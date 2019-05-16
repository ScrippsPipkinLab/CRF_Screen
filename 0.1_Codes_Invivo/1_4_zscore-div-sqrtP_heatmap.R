########## Converted Data analysis ##########
# Author: Huitian (Yolanda) Diao
# May 15th, 2019
# For t-test with percentiles

########## Libraries ##########
library(BSDA)
library(dplyr)
library(tidyverse)
library(magrittr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

########## Self-defined functions ##########
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

in_vec <- function(refvec, vecx){
  out_vec <- numeric(length(vecx))
  for (x in c(1:length(vecx))){
    if (vecx[x] %in% refvec){
      out_vec[x] <- 1
    }
    else {
      out_vec[x] <- 0
    }
  }
  return(out_vec)
}

in_vec_name <- function(refvec, vecx){
  out_vec <- character(length(vecx))
  for (x in c(1:length(vecx))){
    if (vecx[x] %in% refvec){
      out_vec[x] <- vecx[x]
    }
    else {
      out_vec[x] <- ""
    }
  }
  return(out_vec)
}


########## Main ##########
###----- Calculate z-score divided by p-value
if (FALSE){
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/5_zscore_div_sqrt_pval"
  setwd(wk.dir)
  
  z.p.file <- "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/4_t-test_by_gene/all_z-score_p.csv"
  z.p.tb <- read_csv(z.p.file)
  colnames(z.p.tb)
  
  z.p.tb <- z.p.tb %>%
    mutate(InputMinusAvg_p_sqrt = sqrt(InputMinusAvg_p + 0.01)) %>%
    mutate(Q3minusOther_p_sqrt = sqrt(Q3minusOther_p + 0.01)) %>%
    mutate(Q4minusQ1_p_sqrt = sqrt(Q4minusQ1_p + 0.01)) %>%
    mutate(InputMinusAvg_z_p = InputMinusAvg_z/InputMinusAvg_p_sqrt) %>%
    mutate(Q3minusOther_z_p = Q3minusOther_z / Q3minusOther_p_sqrt) %>%
    mutate(Q4minusQ1_z_p = Q4minusQ1_z / Q4minusQ1_p_sqrt) %>%
    select(gene_name, InputMinusAvg_z_p, Q3minusOther_z_p, Q4minusQ1_z_p)
  
  colnames(z.p.tb)
  
  write_csv(z.p.tb, "all_z-score_div_sqrt-p.csv")
}


###----- Heatmap
if (FALSE) {
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/5_zscore_div_sqrt_pval"
  setwd(wk.dir)
  
  z.p.file <- "all_z-score_div_sqrt-p.csv"
  z.p.tb <- read_csv(z.p.file)
  z.p.tb <- z.p.tb %>% 
    mutate(Q4minusQ1 = abs(Q4minusQ1_z_p) / Q4minusQ1_z_p * sqrt(abs(Q4minusQ1_z_p))) %>%
    mutate(Q3minusOther = abs(Q3minusOther_z_p) / Q3minusOther_z_p * sqrt(abs(Q3minusOther_z_p))) %>%
    mutate(InputMinusAvg = abs(InputMinusAvg_z_p) / InputMinusAvg_z_p * sqrt(abs(InputMinusAvg_z_p))) %>%
    select(gene_name, Q4minusQ1, Q3minusOther, InputMinusAvg)
  write_csv(z.p.tb, "all_z-score_div_sqrt-p_sqrt.csv")
  z.p.tb <- z.p.tb %>% column_to_rownames("gene_name") 
  
  
  
  col.pal <- colorRampPalette(brewer.pal(n=11, name="RdBu"))
  col.use <- rev(col.pal(50))
  z.p.heatmap <- pheatmap(z.p.tb , color=col.use, cluster_cols = FALSE, show_rownames = FALSE) 
  
  save_pheatmap_pdf(z.p.heatmap, "zscore_div_sqrt-pval.pdf", 5, 45)
}



##########-------------------- Bar plot
if (FALSE) {
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/5_zscore_div_sqrt_pval"
  setwd(wk.dir)
  
  z.p.file <- "all_z-score_div_sqrt-p_sqrt.csv"
  z.p.tb <- read_csv(z.p.file)
  #z.p.tb <- z.p.tb %>% column_to_rownames("gene_name")
  
  #####---------- Genes to annotate
  anno.vec <- c("Tbx21", "Prdm1", "Id2", "Runx3", "Ncor1", "Tet2", "Mbd2", 
                "Ezh2", "Suv39h1", "Dnmt3a", "Kdm2b", "Rpa3", "Runx3", 
                "Ing2", "Ing3", "Ing4", "Ing5", "Bop1",
                "Cd4", "Cd14")
  
  
  #####---------- Q4 minus Q1
  out.name <- "Q4minusQ1.bar.pdf"
  # Rank order
  z.p.tb <- z.p.tb %>% arrange(Q4minusQ1)
  z.p.tb <- within(z.p.tb, z.p.tb$gene_name <- factor(z.p.tb$gene_name, levels=z.p.tb$gene_name))
  
  # Set color for top and bottom quarter
  col_panel <- c("steelblue1", "lightgrey", "indianred1")
  qt <- as.integer(floor(nrow(z.p.tb)/4))
  col.vec <- rep(col_panel[1], qt)
  col.vec <- c(col.vec, rep(col_panel[2], nrow(z.p.tb)-2*qt))
  col.vec <- c(col.vec, rep(col_panel[3], qt))
  z.p.tb$color_use <- col.vec
  
  # Select annotations
  z.p.tb <- z.p.tb %>% 
    mutate(pointsize = in_vec(anno.vec, gene_name)) %>%
    mutate(annoname = in_vec_name(anno.vec, gene_name))
  
  # Plot
  bar.plot <- ggplot(z.p.tb, aes(z.p.tb$gene_name, z.p.tb$Q4minusQ1, fill=col.vec)) +
    geom_col(alpha=0.7) +
    geom_point(size=z.p.tb$pointsize, stroke = 0) +
    scale_fill_manual(values=col_panel) +
    geom_text_repel(aes(label=annoname), force=5, min.segment.length=0) +
    coord_flip() +
    scale_y_continuous(position = "right", limits=c(-5.2, 5.2)) +
    geom_hline(yintercept=0, size=0.25) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white",colour = "white", size = 0.5, linetype = "solid"),
          axis.line.x = element_line(colour="black", size=0.5), axis.title.x = element_blank(), 
          axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.title.y = element_blank(),
          legend.position = "none")
  bar.plot

  ggsave(out.name, width=6, height=9, units="cm")
  
  #####---------- Q3 minus other
  out.name <- "Q3minusOther.bar.pdf"
  # Rank order
  z.p.tb <- z.p.tb %>% arrange(Q3minusOther)
  z.p.tb <- within(z.p.tb, z.p.tb$gene_name <- factor(z.p.tb$gene_name, levels=z.p.tb$gene_name))
  
  # Set color for top and bottom quarter
  col_panel <- c("steelblue1", "lightgrey", "indianred1")
  qt <- as.integer(floor(nrow(z.p.tb)/4))
  col.vec <- rep(col_panel[1], qt)
  col.vec <- c(col.vec, rep(col_panel[2], nrow(z.p.tb)-2*qt))
  col.vec <- c(col.vec, rep(col_panel[3], qt))
  z.p.tb$color_use <- col.vec
  
  # Select annotations
  z.p.tb <- z.p.tb %>% 
    mutate(pointsize = in_vec(anno.vec, gene_name)) %>%
    mutate(annoname = in_vec_name(anno.vec, gene_name))
  
  # Plot
  bar.plot <- ggplot(z.p.tb, aes(z.p.tb$gene_name, z.p.tb$Q3minusOther, fill=col.vec)) +
    geom_col(alpha=0.7) +
    geom_point(size=z.p.tb$pointsize, stroke = 0) +
    scale_fill_manual(values=col_panel) +
    geom_text_repel(aes(label=annoname), force=5, min.segment.length=0) +
    coord_flip() +
    scale_y_continuous(position = "right", limits=c(-5.2, 5.2)) +
    geom_hline(yintercept=0, size=0.25) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white",colour = "white", size = 0.5, linetype = "solid"),
          axis.line.x = element_line(colour="black", size=0.5), axis.title.x = element_blank(), 
          axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.title.y = element_blank(),
          legend.position = "none")
  bar.plot
  
  ggsave(out.name, width=6, height=9, units="cm")
  
  #####---------- Input v.s. output
  out.name <- "InputMinusAvg.bar.pdf"
  # Rank order
  z.p.tb <- z.p.tb %>% arrange(z.p.tb$InputMinusAvg)
  z.p.tb <- within(z.p.tb, z.p.tb$gene_name <- factor(z.p.tb$gene_name, levels=z.p.tb$gene_name))
  
  # Set color for top and bottom quarter
  col_panel <- c("steelblue1", "lightgrey", "indianred1")
  qt <- as.integer(floor(nrow(z.p.tb)/4))
  col.vec <- rep(col_panel[1], qt)
  col.vec <- c(col.vec, rep(col_panel[2], nrow(z.p.tb)-2*qt))
  col.vec <- c(col.vec, rep(col_panel[3], qt))
  z.p.tb$color_use <- col.vec
  
  # Select annotations
  z.p.tb <- z.p.tb %>% 
    mutate(pointsize = in_vec(anno.vec, gene_name)) %>%
    mutate(annoname = in_vec_name(anno.vec, gene_name))
  
  # Plot
  bar.plot <- ggplot(z.p.tb, aes(z.p.tb$gene_name, z.p.tb$InputMinusAvg, fill=col.vec)) +
    geom_col(alpha=0.7) +
    geom_point(size=z.p.tb$pointsize, stroke = 0) +
    scale_fill_manual(values=col_panel) +
    geom_text_repel(aes(label=annoname), force=5, min.segment.length=0) +
    coord_flip() +
    scale_y_continuous(position = "right", limits=c(-5.2, 5.2)) +
    geom_hline(yintercept=0, size=0.25) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white",colour = "white", size = 0.5, linetype = "solid"),
          axis.line.x = element_line(colour="black", size=0.5), axis.title.x = element_blank(), 
          axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.title.y = element_blank(),
          legend.position = "none")
  bar.plot
  
  ggsave(out.name, width=6, height=9, units="cm")
  
}














