######################################## Correlate In vivo screen result with ScRNAseq ########################################
# Author: Huitian (Yolanda) Diao
# Aug. 6th, 2019

######################################## Libraries ########################################
library(dplyr)
library(tidyverse)

cvt_name <- function(vecx, in_names, out_names) {
  outvec <- c()
  for (i in c(1:length(vecx))){
    if (vecx[i] %in% in_names){
      idx <- match(vecx[i], in_names)
      outvec <- c(outvec, out_names[idx])
    } else {
      outvec <- c(outvec, vecx[i])
    }
  }
  return(outvec)
}

######################################## Main ########################################
wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/3_0_corr_Sc"
setwd(wk.dir)
###----- Plot sc expr on CRF layout
tsne.file <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/1_Amt_compile/Amt_normbycontrolZP_t-test.by.geneavg_z-score_tsne_per5.csv"
sc.expr.file <- "/Volumes/Yolanda/Exp334CD25KOSc/5_CR/1_PAGA_cluster_expr/all_norm_counts_named_nbPctl_z_avgLouvain_Z_CR.csv"
alt.name.file <- "/Volumes/Yolanda/CRF_Screen/Ref/CRF_alternative_gn.csv"

# Merge data tibbles
tsne.tb <- read_csv(tsne.file)
sc.expr.tb <- read_csv(sc.expr.file)
alt.name.tb <- read_csv(alt.name.file)
sc.expr.tb$gene_name <- cvt_name(sc.expr.tb$gene_name,  alt.name.tb$gene_name, alt.name.tb$Alternative)
plot.tb <- inner_join(tsne.tb, sc.expr.tb, by="gene_name")

# Change colnames
plot.colnames <- colnames(plot.tb)
new.colnames <- c()
for (i in plot.colnames) {
  if ((i != "gene_name") && (i != "t1") && (i != "t2")) {
    new.colnames <- c(new.colnames, paste("C", i, sep=""))
  } else(
    new.colnames <- c(new.colnames, i)
  )
}
colnames(plot.tb) <- new.colnames

# Plot
for (i in c(0:9)) {
  #i <- 0
  Ci <- paste("C", as.character(i), sep="")
  Ci_idx <- match(Ci, colnames(plot.tb))
  Ci_z <- as.numeric(unlist(plot.tb[Ci_idx]))
  Ci_new_z <- as.numeric(unlist(scale(Ci_z)))
  Ci_new_z_sqrt <- sqrt(abs(Ci_new_z)) * Ci_new_z / abs(Ci_new_z)
  scplot <- ggplot(data=plot.tb, aes(x=t1, y=t2, color=Ci_new_z_sqrt)) +
    geom_point(stroke=0, size=3, alpha=0.7) +
    scale_color_gradient2(low="dodgerblue3", mid="white",high="firebrick3", midpoint = 0, limits=c(-2,2), name=Ci) +
    theme(panel.background = element_rect(fill = "grey50"),
          panel.grid.major = element_line(size = 0.2, linetype = 'solid',colour = "black"), 
          panel.grid.minor = element_line(size = 0.1, linetype = 'solid',colour = "black"))
  out_name <- paste(Ci, "_CR_sqrt_z.pdf", sep="")
  ggsave(out_name, scplot, device="pdf", units="cm", width=12, height=10)
}

  

###----- Select based on CD25 and CD127
wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/3_0_corr_Sc"
setwd(wk.dir)

p.val.file <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/1_Amt_compile/Amt_normbycontrolZP_t-test.by.genet-test_p-value.csv"
z.file <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/1_Amt_compile/Amt_normbycontrolZP_t-test.by.geneavg_z-score.csv"
p.val.tb <- read_csv(p.val.file)
z.tb <- read_csv(z.file)

p.val.tb.cd25.cd127.sig <- p.val.tb %>% filter(CD25GeoMean_100U <= 0.05)  %>% filter(CD127GeoMean_100U <= 0.05)
z.tb.cd25.cd127 <- z.tb %>% filter(CD25GeoMean_100U < 0 ) %>% filter(CD127GeoMean_100U > 0)
intersect(z.tb.cd25.cd127$gene_name, p.val.tb.cd25.cd127.sig$gene_name)
