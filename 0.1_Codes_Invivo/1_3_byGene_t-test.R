########## Converted Data analysis ##########
# Author: Huitian (Yolanda) Diao
# May 13th, 2019
# For t-test with percentiles

########## Libraries ##########
library(BSDA)
library(dplyr)
library(tidyverse)
library(magrittr)

########## Self-defined functions ##########
uniformCase <- function(vec.x){
  #vec.x <- c("CD4.1", "cd4.2", "RUNX3.2")
  new.vec <- character(length(vec.x))
  for (x in 1:length(vec.x)){
    new.vec[x] <- paste(toupper(substring(vec.x[x], 1, 1)), tolower(substring(vec.x[x], 2)), sep="", collapse=" ")
  }
  return(new.vec)
}


shRNAtoGene <- function(vec.x){
  #vec.x <- c("Smarce1.1", "Smarce1.2",  "Smarce1.3",  "Smarce1.4", "Msl3.1")
  new.vec <- character(length(vec.x))
  for (x in 1:length(vec.x)){
    new.vec[x] <- unlist(strsplit(vec.x[x], "\\."))[1]
  }
  return(new.vec)
}

tByGene <- function(file.x){
  #file.x <- "/Users/yolandatiao/Desktop/CRF_Screen/1_2_normtocontrol_ZP/AmtGFP/AmtGFP_CD25GeoMean_10U_shRNA_normbycontrolZP.csv"
  out.name <- gsub(".csv", "_t-test.by.gene.csv",file.x)
  x.df <- read.csv(file.x)
  shRNA.vec <- as.character(unlist(x.df$shRNA))
  gene.vec <- unique(shRNAtoGene(shRNA.vec))
  gene.vec <- sort(gene.vec)
  all.p.vec <- as.numeric(unlist(x.df$value))
  out.df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(out.df) <- c("gene_name", "avg_percentile", "t-test_p-value")
  for (gene in gene.vec){
    gene.df <- subset(x.df, grepl(gene, x.df$shRNA))
    if ( nrow(gene.df) > 1) {
      gene.p.vec <- as.numeric(unlist(gene.df$value))
      gene.t.test <- t.test(x = gene.p.vec, y = all.p.vec)
      new.row <- c(gene, mean(gene.p.vec), gene.t.test$p.value)
      out.df[nrow(out.df) + 1,] = new.row     
    }
  }
  write.csv(out.df, out.name, row.names = FALSE)
}

########## Main ##########
###----- Combine pools
if (FALSE) {
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/2_Gate_comparisons/byConstruct"
  setwd(wk.dir)
  p1_inputAvg <- read_csv("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/2_Gate_comparisons/byConstruct/P1-7_InputMinusAvg.csv")
  p1_q3_other <- read_csv("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/2_Gate_comparisons/byConstruct/P1-7_Q3minusOther.csv")
  p1_q4_q1 <- read_csv("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/2_Gate_comparisons/byConstruct/P1-7_Q4minusQ1.csv")
  p1_inputAvg <- p1_inputAvg %>% mutate(pool="p1")
  p1_q3_other <- p1_q3_other %>% mutate(pool="p1")
  p1_q4_q1 <- p1_q4_q1 %>% mutate(pool="p1")
  
  p2_inputAvg <- read_csv("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/2_Gate_comparisons/byConstruct/P8-14_InputMinusAvg.csv")
  p2_q3_other <- read_csv("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/2_Gate_comparisons/byConstruct/P8-14_Q3minusOther.csv")
  p2_q4_q1 <- read_csv("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/2_Gate_comparisons/byConstruct/P8-14_Q4minusQ1.csv")
  p2_inputAvg <- p2_inputAvg %>% mutate(pool="p2")
  p2_q3_other <- p2_q3_other %>% mutate(pool="p2")
  p2_q4_q1 <- p2_q4_q1 %>% mutate(pool="p2")
  
  p3_inputAvg <- read_csv("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/2_Gate_comparisons/byConstruct/P15-21_InputMinusAvg.csv")
  p3_q3_other <- read_csv("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/2_Gate_comparisons/byConstruct/P15-21_Q3minusOther.csv")
  p3_q4_q1 <- read_csv("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/2_Gate_comparisons/byConstruct/P15-21_Q4minusQ1.csv")
  p3_inputAvg <- p3_inputAvg %>% mutate(pool="p3")
  p3_q3_other <- p3_q3_other %>% mutate(pool="p3")
  p3_q4_q1 <- p3_q4_q1 %>% mutate(pool="p3")
  
  inputAvg <- bind_rows(p1_inputAvg, p2_inputAvg, p3_inputAvg)
  inputAvg <- inputAvg %>% select(shRNA, inputMinusAvg_nbPctg, pool) %>% set_colnames(c("shRNA", "value", "pool")) %>% mutate(shRNA = uniformCase(shRNA))
  q3_other <- bind_rows(p1_q3_other, p2_q3_other, p3_q3_other)
  q3_other <- q3_other %>% select(shRNA, q3minusOther_nbPctg, pool) %>% set_colnames(c("shRNA", "value", "pool")) %>% mutate(shRNA = uniformCase(shRNA))
  q4_q1 <- bind_rows(p1_q4_q1, p2_q4_q1, p3_q4_q1)
  q4_q1 <- q4_q1 %>% select(shRNA, q4minusq1_nbPctg, pool) %>% set_colnames(c("shRNA", "value", "pool")) %>% mutate(shRNA = uniformCase(shRNA))
  
  write_csv(inputAvg, "InputMinusAvg.csv")
  write_csv(q3_other, "Q3minusOther.csv")
  write_csv(q4_q1, "Q4minusQ1.csv")
}


###----- Calculate p-value of each gene
if (FALSE){
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/2_Gate_comparisons/byConstruct"
  setwd(wk.dir)
  tByGene("InputMinusAvg.csv")
  tByGene("Q3minusOther.csv")
  tByGene("Q4minusQ1.csv")
}


















