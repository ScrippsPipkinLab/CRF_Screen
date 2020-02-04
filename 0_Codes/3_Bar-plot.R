######################################## Bar plot for in vivo screen results ########################################
# Author: Huitian (Yolanda) Diao
# Nov 14th, 2019

######################################## Libraries ########################################
library(dplyr)
library(enrichplot)
library(tidyverse)
library(ggrepel)

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
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
      out_vec[x] <- as.character(vecx[x])
    }
    else {
      out_vec[x] <- ""
    }
  }
  return(out_vec)
}


######################################## Main ########################################
##########-------------------- Bar plot
if (TRUE) {
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/1_Amt_compile"
  setwd(wk.dir)
  
  z.p.file <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/2_Amt_compile_cluster/0_bbplot_source/CD25GeoMean_100U_tsne-z-p_cbs.csv" 
  z.p.tb <- read_csv(z.p.file)
  #z.p.tb <- z.p.tb %>% column_to_rownames("gene_name")
  z.p.tb
  #####---------- Genes to annotate
  #anno.vec <- c("Tbx21", "Prdm1", "Id2", "Runx3", "Ncor1", "Tet2", "Mbd2", 
  #              "Ezh2", "Suv39h1", "Dnmt3a", "Kdm2b", "Rpa3", "Runx3", 
  #              "Ing2", "Ing3", "Ing4", "Ing5", "Bop1",
  #              "Cd4", "Cd14")
  #anno.vec <- c("Myst3", "Myst4", "Brpf1", "Ing5", "Ing4", "Ing3", 
  #              "Mll1", "Wdr5", "Rbbp5", "Ash2l", "Dpy30",
  #              "Cd4", "Runx3", "Tbx21",
  #              "Cxxc1", "Paf1")
  anno.vec <- c("ACTL6A", "ARID1A", "ARID1B", "ARID2", "BRD7", "BRD9", "PBRM1", "PHF10", 
                "SMARCA2", "SMARCA4", "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD1", "SMARCD2", 
                "SMARCD3", "SMARCE1", "Actl6a", "Actl6b", "Arid1a", "Arid1b", "Brd9", "Smarca2", 
                "Smarca4", "Smarcb1", "Smarcc1", "Smarcc2", "Smarcd1", "Smarcd2", "Smarcd3", "Smarce1",
                "CHD3", "CHD4", "CHD5", "HDAC1", "HDAC2", "KDM1A", 
                
                "MBD2", "MBD3", "MTA1", "MTA2", "MTA3", "RBBP4", "RBBP7", "SIN3A", "SIN3B",
                "Chd3", "Chd4", "Chd5", "Hdac1", "Hdac2", "Mbd3", "Mta1", "Mta2", 
                "Mta3", "Rbbp4", "Rbbp7", 
                
                "Cd4", "Runx3", "Tbx21", "Chd7")
  anno.vec <- tolower(anno.vec)
  anno.vec <- as.character(sapply(anno.vec, simpleCap))
  
  
  #####---------- CD25 100U
  out.name <- "CD25_100U.bar.pdf"
  
  # Rank order
  z.p.tb <- z.p.tb %>% filter(pValue <= 0.05) %>% arrange(zScore)
  #z.p.tb <- within(z.p.tb, z.p.tb$gene_name <- factor(z.p.tb$gene_name, levels=z.p.tb$gene_name))
  z.p.tb$gene_name <- factor(z.p.tb$gene_name, levels=z.p.tb$gene_name)
  
  # Set color for top and bottom quarter
  col_panel <- c( "deepskyblue", "tomato", "snow2")
  col.vec <- c()
  for (i in z.p.tb$zScore){
    if(i>0) {
      col.vec <- c(col.vec, col_panel[2])
    } else if (i<0) {
      col.vec <- c(col.vec, col_panel[1])
    } else {
      col.vec <- c(col.vec, col_panel[3])
    }
  }
  z.p.tb$color_use <- col.vec
  
  # Select annotations
  z.p.tb <- z.p.tb %>% 
    mutate(pointsize = in_vec(anno.vec, gene_name)) %>%
    mutate(annoname = in_vec_name(anno.vec, gene_name))
  
  # Plot
  bar.plot <- ggplot(z.p.tb, aes(z.p.tb$gene_name, z.p.tb$zScore, fill=col.vec)) +
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
  
  ggsave(out.name, width=6, height=6, units="cm")
}

