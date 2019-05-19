########## Converted Data analysis ##########
# Author: Huitian (Yolanda) Diao
# May 15th, 2019
# For t-test with percentiles
# Reference: https://rdrr.io/bioc/TCGAbiolinks/man/getAdjacencyBiogrid.html

########## Libraries ##########
library(dplyr)
library(tidyverse)
library(TCGAbiolinks)

########## Self-defined functions ##########
convertGN <- function(in.vec, cov.old, cov.new){
  for (x in c(1:length(in.vec))){
    if (in.vec[x] %in% cov.old){
      idx <- match(in.vec[x], cov.old)
      in.vec[x] <- cov.new[idx]
    }
  }
  return(in.vec)
}


########## Main ##########
wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein"
setwd(wk.dir)

###----- Read inputs
tgt.file <- "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/5_zscore_div_sqrt_pval/target.heatmap.clusters.csv"
tgt.tb <- read_csv(tgt.file)

crf.gn.alt.file <- "/Volumes/Yolanda/CRF_Screen/Ref/CRF_alternative_gn.csv"
crf.gn.alt.tb <- read_csv(crf.gn.alt.file)

biogrid.ref <- "/Volumes/Yolanda/CRF_Screen/Ref/BIOGRID.use_simp.txt"
biogrid.tb <- read_delim(biogrid.ref, "\t")

###----- Convert mouse gene names to all upper case (match with human)
tgt.tb <- tgt.tb %>% 
  mutate(gene_name = convertGN(gene_name, crf.gn.alt.tb$Alternative, crf.gn.alt.tb$gene_name)) %>%
  mutate(gene_name = toupper(gene_name))
biogrid.tb <- biogrid.tb %>% 
  mutate(gene_name = toupper(OFFICIAL_SYMBOL_A)) %>%
  mutate(gene_name2 = toupper(OFFICIAL_SYMBOL_B)) %>%
  select(gene_name, gene_name2)

###----- Select mouse / human genes from biogrid spreadsheet
biogrid.tb.tgt <- biogrid.tb %>%
  filter(gene_name %in% tgt.tb$gene_name) %>%
  filter(gene_name2 %in% tgt.tb$gene_name) %>%
  distinct() %>%
  left_join(tgt.tb, by="gene_name") %>%
  mutate(interation="pp")


write_csv(biogrid.tb.tgt, "CRF_target_biogrid.csv")


