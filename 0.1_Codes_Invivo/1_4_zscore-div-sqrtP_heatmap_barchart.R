########## Converted Data analysis ##########
# Author: Huitian (Yolanda) Diao
# May 15th, 2019
# For t-test with percentiles
# Reference: https://bioconductor.statistik.tu-dortmund.de/packages/3.1/bioc/vignettes/ComplexHeatmap/inst/doc/ComplexHeatmap.html

########## Libraries ##########
library(dplyr)
library(tidyverse)

library(BSDA)
library(magrittr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)

#BiocManager::install("org.Mm.eg.db")
library("org.Mm.eg.db")
library(clusterProfiler)
library(ggplot2)


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

# Save dendrogram k-means clusters from ComplexHeatmap
# Reference: https://github.com/jokergoo/ComplexHeatmap/issues/136
k_means_save <- function(c_heatmap, c_tb, out_name) {
  #c_heatmap <- target.heatmap
  #c_tb <- target.tb
  #out_name <- "target_cluster.csv"
  r.dend <- row_dend(c_heatmap)
  rcl.list <- row_order(c_heatmap)
  #lapply(rcl.list, function(x) length(x))
  
  for (i in 1:length(row_order(c_heatmap))){
    if (i == 1) {
      clu <- t(t(row.names(c_tb[row_order(c_heatmap)[[i]],])))
      out <- cbind(clu, paste("cluster", i, sep=""))
      colnames(out) <- c("GeneID", "Cluster")
    } else {
      clu <- t(t(row.names(c_tb[row_order(c_heatmap)[[i]],])))
      clu <- cbind(clu, paste("cluster", i, sep=""))
      out <- rbind(out, clu)
    }
  }
  
  write.csv(out, file=out_name, row.names=FALSE)
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
  
  # Convert adjusted z-score by sqrt
  z.p.tb <- z.p.tb %>% 
    mutate(Q4minusQ1 = abs(Q4minusQ1_z_p) / Q4minusQ1_z_p * sqrt(abs(Q4minusQ1_z_p))) %>%
    mutate(Q3minusOther = abs(Q3minusOther_z_p) / Q3minusOther_z_p * sqrt(abs(Q3minusOther_z_p))) %>%
    mutate(InputMinusAvg = abs(InputMinusAvg_z_p) / InputMinusAvg_z_p * sqrt(abs(InputMinusAvg_z_p))) %>%
    select(gene_name, Q4minusQ1, Q3minusOther, InputMinusAvg)
  write_csv(z.p.tb, "all_z-score_div_sqrt-p_sqrt.csv")
  
  # Plot heatmap
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

##########-------------------- Seperated heatmap
if (FALSE) {
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/5_zscore_div_sqrt_pval"
  setwd(wk.dir)
  
  z.p.file <- "all_z-score_div_sqrt-p_sqrt.csv"
  z.p.tb <- read_csv(z.p.file)
  #z.p.tb <- z.p.tb %>% column_to_rownames("gene_name")
  
  colnames(z.p.tb)
  
  ###----- Find genes that are on top / bottom of list in different cases
  #--- Q4 v.s. Q1
  z.p.tb <- z.p.tb %>% arrange(Q4minusQ1)
  qt <- as.integer(floor(nrow(z.p.tb)/8))
  q4_q1_dn <- z.p.tb$gene_name[1: qt]
  q4_q1_up <- z.p.tb$gene_name[(nrow(z.p.tb)-qt+1): nrow(z.p.tb)]
  #--- Q3 v.s. Other
  z.p.tb <- z.p.tb %>% arrange(Q3minusOther)
  qt <- as.integer(floor(nrow(z.p.tb)/8))
  q3_o_dn <- z.p.tb$gene_name[1: qt]
  q3_o_up <- z.p.tb$gene_name[(nrow(z.p.tb)-qt+1): nrow(z.p.tb)]
  #--- Input v.s. Avg
  z.p.tb <- z.p.tb %>% arrange(InputMinusAvg)
  qt <- as.integer(floor(nrow(z.p.tb)/8))
  in_a_dn <- z.p.tb$gene_name[1: qt]
  in_a_up <- z.p.tb$gene_name[(nrow(z.p.tb)-qt+1): nrow(z.p.tb)]
  
  ###----- All genes that are on top / bottom in any case
  all.tgt <- c(q4_q1_dn, q4_q1_up, q3_o_dn, q3_o_up, in_a_dn, in_a_up)
  all.tgt <- unique(all.tgt)
  length(all.tgt)
  
  #####---------- Heatmaps
  #col.pal <- colorRampPalette(brewer.pal(n=11, name="RdBu"))
  #col.use <- rev(col.pal(50))
  ###----- Heatmap of all targets
  target.tb <- z.p.tb %>% filter(gene_name %in% all.tgt) %>% column_to_rownames("gene_name") 
  #target.heatmap <- pheatmap(target.tb , color=col.use, cluster_cols = FALSE, 
  #                           clustering_distance_rows="canberra") 
  
  set.seed(123)
  pdf("target_cluster_2.pdf", width=4, height=24)
  target.heatmap <-   Heatmap(target.tb, name="foo", km=6, cluster_columns=FALSE, 
                              clustering_distance_rows = "pearson")
  target.heatmap
  dev.off()
  
  ###----- Save dendrogram clusters
  k_means_save(target.heatmap, target.tb, "target_cluster.csv")
  write_csv(target.tb, "target.tb.csv")
  
}

##########-------------------- Pathway
if (FALSE){
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/5_zscore_div_sqrt_pval"
  setwd(wk.dir)
  
  if (FALSE){
    z.p.file <- "all_z-score_div_sqrt-p_sqrt.csv"
    z.p.tb <- read_csv(z.p.file)
    qt <- as.integer(floor(nrow(z.p.tb)/4))
    
    #####---------- Q3 minus other
    z.p.tb <- z.p.tb %>% arrange(Q3minusOther)
    z.p.tb <- within(z.p.tb, z.p.tb$gene_name <- factor(z.p.tb$gene_name, levels=z.p.tb$gene_name))
    
    q3_o_dn <- z.p.tb$gene_name[1:qt]
    q3_o_up <- z.p.tb$gene_name[(nrow(z.p.tb)-qt+1): nrow(z.p.tb)]
    
    #####---------- Q4 minus Q1
    z.p.tb <- z.p.tb %>% arrange(Q4minusQ1)
    z.p.tb <- within(z.p.tb, z.p.tb$gene_name <- factor(z.p.tb$gene_name, levels=z.p.tb$gene_name))
    
    q4_q1_dn <- z.p.tb$gene_name[1:qt]
    q4_q1_up <- z.p.tb$gene_name[(nrow(z.p.tb)-qt+1): nrow(z.p.tb)]
    
    #####---------- Input v.s. output
    z.p.tb <- z.p.tb %>% arrange(z.p.tb$InputMinusAvg)
    z.p.tb <- within(z.p.tb, z.p.tb$gene_name <- factor(z.p.tb$gene_name, levels=z.p.tb$gene_name))
    
    in_a_dn <- z.p.tb$gene_name[1:qt]
    in_a_up <- z.p.tb$gene_name[(nrow(z.p.tb)-qt+1): nrow(z.p.tb)]
    
    out.tb <- tibble(q3_o_dn=q3_o_dn, q3_o_up=q3_o_up, 
                     q4_q1_dn=q4_q1_dn, q4_q1_up=q4_q1_up, 
                     in_a_dn=in_a_dn, in_a_up=in_a_up)
    write_csv(out.tb, "topQuarter.csv")
  }
  
  in.df <- read.csv("topQuarter.csv")
  gn_list_names <- colnames(in.df)
  q3_o_dn <- in.df$q3_o_dn
  q3_o_up <- in.df$q3_o_up
  q4_q1_dn <- in.df$q4_q1_dn
  q4_q1_up <- in.df$q4_q1_up
  in_a_dn <- in.df$in_a_dn
  in_a_up <- in.df$in_a_up
  gn_list <- list(q3_o_dn, q3_o_up, q4_q1_dn, q4_q1_up, in_a_dn, in_a_up)
  
  for (x in c(1:6)){
    #x <- 1
    i <- paste(gn_list_names[x], sep="")
    genes.i <- as.character(unlist(gn_list[x]))
    
    genes.i.id <- select(org.Mm.eg.db, genes.i, c("ENTREZID"), "ALIAS")
    #genes.i.id$ENTREZID
    
    egoBP <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "none", pvalueCutoff = 0.05, readable = TRUE)
    egoCC <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "CC", pAdjustMethod = "none", pvalueCutoff = 0.05, readable = TRUE)
    egoMF <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "MF", pAdjustMethod = "none", pvalueCutoff = 0.05, readable = TRUE)
    
    # Dotplot visualization
    if (!is.null(egoBP)){
      pdf.name <- paste(i,"_BP_dotplot.pdf",sep="")
      csv.name <- paste(i,"_BP_dotplot.csv",sep="")
      write.csv(egoBP@result, file=csv.name, row.names=FALSE)
      egoBP.dotplot <- dotplot(egoBP, x="count", showCategory=25)
      ggsave(pdf.name, egoBP.dotplot, device = "pdf", width = 30, height = 20, units = "cm")  
      
    }
    if(!is.null(egoCC)){
      csv.name <- paste(i,"_CC_dotplot.csv",sep="")
      pdf.name <- paste(i,"_CC_dotplot.pdf",sep="")
      write.csv(egoCC@result, file=csv.name, row.names=FALSE)
      egoCC.dotplot <- dotplot(egoCC, x="count", showCategory=25)
      ggsave(pdf.name, egoCC.dotplot, device = "pdf", width = 30, height = 20, units = "cm")  
    }
    if(!is.null(egoMF)){
      csv.name <- paste(i,"_MF_dotplot.csv",sep="")
      pdf.name <- paste(i,"_MF_dotplot.pdf",sep="")
      write.csv(egoMF@result, file=csv.name, row.names=FALSE)
      egoMF.dotplot <- dotplot(egoMF, x="count", showCategory=25)
      ggsave(paste(i,"_MF_dotplot.pdf",sep=""), egoMF.dotplot, device = "pdf", width = 30, height = 20, units = "cm")  
    }
  }
  
}















