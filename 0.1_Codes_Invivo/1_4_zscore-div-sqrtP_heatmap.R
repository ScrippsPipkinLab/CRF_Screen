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

########## Self-defined functions ##########
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


########## Main ##########
###----- Calculate z-score divided by p-value
if (FALSE){
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/6_zscore_div_sqrt-pval"
  setwd(wk.dir)
  
  z.p.file <- "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/5_p-val_byGene/all_z-score_p.csv"
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
wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/6_zscore_div_sqrt-pval"
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
















