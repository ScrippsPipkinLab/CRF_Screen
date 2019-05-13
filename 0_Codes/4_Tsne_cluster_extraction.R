######################################## Tsne cluster extraction ########################################
# Author: Huitian (Yolanda) Diao
# April 22nd, 2019
# Reference: https://jmonlong.github.io/Hippocamplus/2018/02/13/tsne-and-clustering/

######################################## Imports ########################################
library(dplyr)

library(igraph)
library(FNN)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(rlang)
library(pheatmap)
library(RColorBrewer)

#BiocManager::install("org.Mm.eg.db")
library("org.Mm.eg.db")
library("clusterProfiler")

######################################## Self-defined fucntions ########################################
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


######################################## Main ########################################


#################### Extract clusters ####################
if (FALSE) {
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/2_Amt_compile_cluster/2_cluster_extraction"
  setwd(wk.dir)
  ###----- Read data
  tsne.file <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/1_Amt_compile/Amt_normbycontrolZP_t-test.by.geneavg_z-score_tsne_per5.csv"
  tsne.all <- read.csv(tsne.file)
  tsne.real <- subset(tsne.all, select=c("t1", "t2"))
  info.real <- tsne.all
  
  # Check tsne
  ggplot(tsne.real, aes(x = t1, y = t2)) + geom_point(alpha = 0.1) + theme_bw()
  
  # Baf complex
  baf.mem <- c("Actb", "Actl6a", "Actl6b", "Arid1a", "Arid1b", 
               "Bcl7a", "Bcl7b", "Bcl7c", "Bcl11a", "Bcl11b", 
               "Brd9", "Dpf1", "Dpf2", "Dpf3", "Phf10", "Smarca2", 
               "Smarca4", "Smarcb1", "Smarcc1", "Smarcc2", "Smarcd1", 
               "Smarcd2", "Smarcd3", "Smarce1", "Ss18", "Ss18l1")

  
  ###----- Clustering
  k <- 15
  knn.real <- get.knn(as.matrix(tsne.real), k = k)
  knn.real <- data.frame(from = rep(1:nrow(knn.real$nn.index), k), to = as.vector(knn.real$nn.index), weight = 1/(1 + as.vector(knn.real$nn.dist)))
  
  nw.real <- graph_from_data_frame(knn.real, directed = FALSE)
  nw.real <- simplify(nw.real)
  
  lc.real <- cluster_louvain(nw.real) #, gamma = 0.1
  info.real$louvain <- as.factor(membership(lc.real))
  lc.cent <- info.real %>% group_by(louvain) %>% select(t1, t2) %>% summarize_all(mean)
  
  label.plot <- ggplot(info.real, aes(x = t1, y = t2, colour = louvain, label = gene_name)) + 
    #geom_text() +
    geom_point(alpha = 0.3) + theme_bw() + 
    geom_label_repel(aes(label = louvain), data = lc.cent) + 
    #geom_text_repel(
    #  data = subset(info.real, info.real$gene_name %in% baf.mem),
    #  point.padding = NA, force = 2.5) +
    guides(colour = FALSE)
  label.plot
  
  ggsave("exp36-174_flt_avg-log10tpm_log2fc_CRF_tsne.pdf", width = 12, height = 10, units = "cm")
  
  ###----- Save
  write.csv(info.real, file = "info.real.csv",row.names=FALSE)
}

#################### Pathway ####################
wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/2_Amt_compile_cluster/2_cluster_extraction"
setwd(wk.dir)
info.real <- read.csv("info.real.csv")
for (x in c(1:13)){
  i <- paste("group", x, sep="")
  subset.x <- subset(info.real, info.real$louvain == x)
  genes.i <- as.character(unlist(subset.x$gene_name))
  genes.i <- unique(genes.i)
  
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


#################### Cluster heatmap ####################
### 
if (FALSE) {
  wk.dir <- "/Volumes/Yolanda/RNAseq_Compilation/CRM_analysis_2/1_GeneCluster/0.1_Cluster_heatmap"
  setwd(wk.dir)
  
  in.file <- "/Volumes/Yolanda/RNAseq_Compilation/CRM_analysis_2/1_GeneCluster/0_Cluster_extraction/info.real.csv"
  in.tb <- read_csv(in.file)
  heat.tb <- in.tb %>% select(gene_name, D5EE, D5TE, D8EE, D8TE, D8MP, NAV, EXH, MEM, louvain )
  
  #--- Construct log10tpm heatmap tibble
  tpm.names <- c("D5EE", "D5TE", "D8EE", "D8TE", "D8MP", "NAV", "EXH", "MEM")
  log10.names <- c("gene_name", "louvain")
  for (col.name in tpm.names) {
    new.col.name <- paste(col.name, "log", sep="_")
    log10.names <- c(log10.names, c(new.col.name))
    heat.tb <- heat.tb %>%
      mutate(!!sym(new.col.name) := log10(!!sym(col.name) + 1))
  }
  heat.log10.tb <- heat.tb %>% select(one_of(log10.names))
  # Write output
  write_csv(heat.log10.tb, "log10tpm.csv")
  
  #--- Calculate row Z-score
  if (FALSE) {
    log10.names <- log10.names[3:length(log10.names)]
    log10.names
    z.names <- c("gene_name", "louvain")
    for (col.name in log10.names) {
      new.col.name <- paste(col.name, "z", sep="_")
      z.names <- c(z.names, c(new.col.name))
      heat.log10.tb <- heat.log10.tb %>%
        rowwise() %>%
        mutate( !!sym(new.col.name) := ( !!sym(col.name) -  mean(c(D5EE_log,D5TE_log,D8EE_log,D8TE_log,D8MP_log,NAV_log,EXH_log,MEM_log)) ) /sd(c(D5EE_log,D5TE_log,D8EE_log,D8TE_log,D8MP_log,NAV_log,EXH_log,MEM_log)) )
    }
    heat.log10.z.tb <- heat.log10.tb %>% select(one_of(z.names))
    # Write output
    write_csv(heat.log10.z.tb, "log10tpm_row-zscore.csv")
  }
  
  #--- Calculate column Z-score
  if (FALSE) {
    heat.log10.tb <- heat.log10.tb %>%
      mutate(D5EE_log_z = scale(D5EE_log, center = TRUE, scale = TRUE)) %>%
      mutate(D5TE_log_z = scale(D5TE_log, center = TRUE, scale = TRUE)) %>%
      mutate(D8EE_log_z = scale(D8EE_log, center = TRUE, scale = TRUE)) %>%
      mutate(D8TE_log_z = scale(D8TE_log, center = TRUE, scale = TRUE)) %>%
      mutate(D8MP_log_z = scale(D8MP_log, center = TRUE, scale = TRUE)) %>%
      mutate(NAV_log_z = scale(NAV_log, center = TRUE, scale = TRUE)) %>%
      mutate(EXH_log_z = scale(EXH_log, center = TRUE, scale = TRUE)) %>%
      mutate(MEM_log_z = scale(MEM_log, center = TRUE, scale = TRUE)) %>%
      select(gene_name, louvain, D5EE_log_z, D5TE_log_z, D8EE_log_z, D8TE_log_z, D8MP_log_z, NAV_log_z, EXH_log_z, MEM_log_z)
    # Write output
    write_csv(heat.log10.tb, "log10tpm_col-zscore.csv")
  }

  
  #####----- Heatmap
  heat.log10.tb <- read_csv("log10tpm_col-zscore.csv")
  heat.log10.tb <- heat.log10.tb %>% column_to_rownames("gene_name")
  heat.log10.tb <- heat.log10.tb %>% select(D5EE_log_z, D5TE_log_z, D8EE_log_z, D8TE_log_z, D8MP_log_z, EXH_log_z, MEM_log_z, NAV_log_z)
  head(heat.log10.tb)
  
  col.pal <- colorRampPalette(brewer.pal(n=11, name="RdBu"), method="linear")
  
  z.heatmap <- pheatmap(heat.log10.tb, color=col.pal(50), cluster_cols = FALSE)
  
  save_pheatmap_pdf(z.heatmap, "CRF_col-zscore_heatmap.pdf", 10, 40)
  
  
  

}



