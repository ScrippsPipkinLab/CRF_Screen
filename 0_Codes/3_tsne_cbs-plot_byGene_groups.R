########## CBS-plot ##########
# Author: Huitian (Yolanda) Diao
# Dec 3rd, 2018
# For plottinng color-bar scatter plot

########## Libraries ##########
library(ggplot2)
library(ggrepel)
library("clusterProfiler")
#BiocManager::install("org.Mm.eg.db")
library("org.Mm.eg.db")

########## Self-defined functions ##########
build_df <- function(tsne.file, z.score.file, p.value.file, key.col){
  #key.col <- "CD127GeoMean_10U"
  tsnf.df <- read.csv(tsne.file)
  z.score.df <- read.csv(z.score.file)
  p.value.df <- read.csv(p.value.file)
  
  new.df <- tsnf.df
  new.df$zScore <- unlist(z.score.df[key.col])
  new.df$pValue <- unlist(p.value.df[key.col])
  
  out.name <- paste(key.col, "_tsnf-z-p.csv")
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

cbs_plot <- function(in.file, slt.list, out.appx){
  in.df <- read.csv(in.file)
  sub.df <- subset(in.df, in.df$gene_name %in% slt.list)
  
  # Save sub plot spreadsheet
  out.csv.name <- gsub(".csv", "", in.file)
  out.csv.name <- paste(out.csv.name, "_", out.appx, ".csv", sep="")
  write.csv(sub.df, out.csv.name, row.names = FALSE)
  
  
  cbs.plot <- ggplot(in.df, aes(x=in.df$t1, y=in.df$t2)) +
    geom_point(alpha = 0.6, shape = 19, color='black', size=2) +
    scale_x_continuous(limits = c(-100, 100)) +
    scale_y_continuous(limits = c(-100, 100)) +
    theme(panel.background = element_rect(fill = "white", colour = "white",size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "grey"), 
          panel.grid.minor = element_line(size = 0.15, linetype = 'solid',colour = "grey"),
          text = element_text(size=10))
  cbs.plot
  cbs.plot <- cbs.plot + 
    layer(data = sub.df, geom="point", position = "identity", 
          stat="identity", mapping=aes(x=sub.df$t1, y=sub.df$t2, color="red", alpha=0.8, size=2))
  #          geom_params = list(fill = "white", color = "steelblue")
    
  cbs.plot
  # Save without annotation
  out.name <- gsub(".csv","", in.file)
  out.name <- paste(out.name, "_", out.appx, ".pdf", sep="")
  out.name
  ggsave(out.name, plot = cbs.plot, device = "pdf", units = "cm", width = 25, height = 25)
}

##########---------- Main
in_file <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/compiled_Amt/per5/Amt_normbycontrolZP_t-test.by.geneavg_z-score_tsne_per5.csv"
in_df <- read.csv(in_file)

if(TRUE){
  ###----- Group1
  slt_df1 <- subset(in_df, in_df$t2>50)
  g1_list <- unlist(slt_df1$gene_name)
  
  
  ###----- Group2
  slt_df2 <- subset(in_df, (in_df$t1>50)&(in_df$t2<0))
  g2_list <- unlist(slt_df2$gene_name)

  
  ###----- Group3
  slt_df3 <- subset(in_df, in_df$t1 < -70)
  g3_list <- unlist(slt_df3$gene_name)

  
  ###----- Group4
  slt_df4 <- subset(in_df, (in_df$t1 < -25)&(in_df$t1 > -70) | ( ((in_df$t1 < 0)&(in_df$t1 > -50)) & ((in_df$t2>0)&(in_df$t2<50)) ) )
  g4_list <- unlist(slt_df4$gene_name)

  ###----- Group5
  slt_df5 <- subset(in_df, in_df$t2 < -30)
  g5_list <- unlist(slt_df5$gene_name)

  
  ###----- Group6
  slt_df6 <- subset(in_df, (in_df$t1 > -25)&(in_df$t1 <20)&(in_df$t2<0)&(in_df$t2> -30))
  g6_list <- unlist(slt_df6$gene_name)

  
  ###----- Group7
  slt_df7 <- subset(in_df, (in_df$t1 > 50)&(in_df$t2> 0))
  g7_list <- unlist(slt_df7$gene_name)

  
  ###----- Group8
  slt_df8 <- subset(in_df, (in_df$t1 > 30)&(in_df$t1< 50)&(in_df$t2<25))
  g8_list <- unlist(slt_df8$gene_name)

  
  ###----- Group9
  slt_df9 <- subset(in_df, (in_df$t1 < 30)&(in_df$t1 > 15)&(in_df$t2<25))
  g9_list <- unlist(slt_df9$gene_name)

  
  ###----- Group10
  slt_df10 <- subset(in_df, (in_df$t1 < 50)&(in_df$t1 > 0)&(in_df$t2<50)&(in_df$t2>25))
  g10_list <- unlist(slt_df10$gene_name)
}

if(FALSE){
  #cbs_plot(in_file, g1_list, "g1")
  #cbs_plot(in_file, g2_list, "g2")
  #cbs_plot(in_file, g3_list, "g3")
  #cbs_plot(in_file, g4_list, "g4")
  #cbs_plot(in_file, g5_list, "g5")
  #cbs_plot(in_file, g6_list, "g6")
  #cbs_plot(in_file, g7_list, "g7")
  #cbs_plot(in_file, g8_list, "g8")
  #cbs_plot(in_file, g9_list, "g9")
  #cbs_plot(in_file, g10_list, "g10")
}


###----- Plot all
if(FALSE){
  all.plot <- ggplot(in.df, aes(x=in_df$t1, y=in_df$t2)) +
    geom_point(alpha = 0.6, shape = 19, color='black', size=2) +
    scale_x_continuous(limits = c(-100, 100)) +
    scale_y_continuous(limits = c(-100, 100)) +
    theme(panel.background = element_rect(fill = "white", colour = "white",size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "grey"), 
          panel.grid.minor = element_line(size = 0.15, linetype = 'solid',colour = "grey"),
          text = element_text(size=10))
  
  
  group_colors <- c("dodgerblue1", "firebrick1","limegreen","orange","slateblue",
                    "turquoise1", "yellow","tan4","black", "palevioletred1")
  
  slt_df1$colors <- rep(group_colors[1], nrow(slt_df1))
  slt_df2$colors <- rep(group_colors[2], nrow(slt_df2))
  slt_df3$colors <- rep(group_colors[3], nrow(slt_df3))
  slt_df4$colors <- rep(group_colors[4], nrow(slt_df4))
  slt_df5$colors <- rep(group_colors[5], nrow(slt_df5))
  slt_df6$colors <- rep(group_colors[6], nrow(slt_df6))
  slt_df7$colors <- rep(group_colors[7], nrow(slt_df7))
  slt_df8$colors <- rep(group_colors[8], nrow(slt_df8))
  slt_df9$colors <- rep(group_colors[9], nrow(slt_df9))
  slt_df10$colors <- rep(group_colors[10], nrow(slt_df10))
  all_df <- rbind(slt_df1, slt_df2, slt_df3, slt_df4, slt_df5,
                  slt_df6, slt_df7, slt_df8, slt_df9, slt_df10)
  
  all_plot <- ggplot(all_df, aes(x=all_df$t1, y=all_df$t2, color=all_df$colors)) +
    geom_point(alpha = 0.6, shape = 19, size=2) +
    scale_x_continuous(limits = c(-100, 100)) +
    scale_y_continuous(limits = c(-100, 100)) +
    theme(panel.background = element_rect(fill = "white", colour = "white",size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "grey"), 
          panel.grid.minor = element_line(size = 0.15, linetype = 'solid',colour = "grey"),
          text = element_text(size=10))
  all_plot
  
  setwd("/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/compiled_Amt")
  ggsave("all_groups.pdf", plot = all_plot, device = "pdf", units = "cm", width = 25, height = 25)
  
  all_plot <- ggplot(all_df, aes(x=all_df$t1, y=all_df$t2, color=all_df$colors)) +
    geom_point(alpha = 0.6, shape = 19, size=2, color="black") +
    scale_x_continuous(limits = c(-100, 100)) +
    scale_y_continuous(limits = c(-100, 100)) +
    theme(panel.background = element_rect(fill = "white", colour = "white",size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "grey"), 
          panel.grid.minor = element_line(size = 0.15, linetype = 'solid',colour = "grey"),
          text = element_text(size=10))
  all_plot
  setwd("/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/compiled_Amt")
  ggsave("al.pdf", plot = all_plot, device = "pdf", units = "cm", width = 25, height = 25)
}


#####---------- Pathway analysis
group_lists <- list(g1_list, g2_list, g3_list, g4_list, g5_list,
                    g6_list, g7_list, g8_list, g9_list, g10_list)

### Go-term dot plot
if (FALSE) {
  for (i in c(1:length(group_lists))) {
    setwd("/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/compiled_Amt")
    genes.i <- as.character(unlist(group_lists[i]))
    genes.i.id <- select(org.Mm.eg.db, genes.i, c("ENTREZID"), "ALIAS")
    genes.i.id$ENTREZID

    egoBP <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
    egoCC <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
    egoMF <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
    
    # Dotplot visualization
    if (!is.null(egoBP)){
      egoBP.dotplot <- dotplot(egoBP, showCategory=25)
      ggsave(paste(i,"_BP_dotplot.pdf",sep=""), egoBP.dotplot, device = "pdf", width = 30, height = 20, units = "cm")      
    }
    if(!is.null(egoCC)){
      egoCC.dotplot <- dotplot(egoCC, showCategory=25)
      ggsave(paste(i,"_CC_dotplot.pdf",sep=""), egoCC.dotplot, device = "pdf", width = 30, height = 20, units = "cm")      
    }
    if(!is.null(egoMF)){
      egoMF.dotplot <- dotplot(egoMF, showCategory=25)
      ggsave(paste(i,"_MF_dotplot.pdf",sep=""), egoMF.dotplot, device = "pdf", width = 30, height = 20, units = "cm")      
    }
  }
}
