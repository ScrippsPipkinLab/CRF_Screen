########## Converted Data analysis ##########
# Author: Huitian (Yolanda) Diao
# Dec 3rd, 2018
# For t-test with percentiles

########## Libraries ##########
library(BSDA)

########## Self-defined functions ##########
simpName <- function(vec.x) {
  new.vec <- character(length(vec.x))
  for (x in 1:length(vec.x)){
    new.x <- gsub("_shRNA_normbycontrolZP_t-test.by.gene.csv", "", vec.x[x])
    new.x <- gsub("Amt_", "", new.x)
    new.x <- gsub("AmtGFP_", "", new.x)
    new.vec[x] <- new.x
  }
  return(new.vec)
}

simpNameOne <- function(str.x) {
  new.x <- gsub("_shRNA_normbycontrolZP_t-test.by.gene.csv", "", str.x)
  new.x <- gsub("Amt_", "", new.x)
  new.x <- gsub("AmtGFP_", "", new.x)
  return(new.x)
}

shRNAtoGene <- function(vec.x){
  #vec.x <- c("Smarce1.1", "Smarce1.2",  "Smarce1.3",  "Smarce1.4", "Msl3.1")
  new.vec <- numeric(length(vec.x))
  for (x in 1:length(vec.x)){
    new.vec[x] <- unlist(strsplit(vec.x[x], "\\."))[1]
  }
  return(new.vec)
}

tByGene <- function(file.x){
  #file.x <- "/Users/yolandatiao/Desktop/CRF_Screen/1_2_normtocontrol_ZP/AmtGFP/AmtGFP_CD25GeoMean_10U_shRNA_normbycontrolZP.csv"
  out.name <- gsub(".csv", "_t-test.by.gene.csv",file.x)
  x.df <- read.csv(file.x) 
  x.df <- subset(x.df, 
                 !(grepl("11", x.df$Position) & grepl("control", x.df$Position)) 
                 | !grepl("Cd4", x.df$shRNA) ) # Exclude weird shCd4 sample
  shRNA.vec <- as.character(unlist(x.df$shRNA))
  gene.vec <- unique(shRNAtoGene(shRNA.vec))
  gene.vec <- sort(gene.vec)
  all.p.vec <- as.numeric(unlist(x.df$precentile))
  out.df <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(out.df) <- c("gene_name", "avg_percentile", "avg_z-score", "t-test_p-value")
  for (gene in gene.vec){
    gene.df <- subset(x.df, grepl(gene, x.df$shRNA))
    if ( nrow(gene.df) > 1) {
      gene.p.vec <- as.numeric(unlist(gene.df$precentile))
      gene.z.vec <- as.numeric(unlist(gene.df$z))
      gene.t.test <- t.test(x = gene.p.vec, y = all.p.vec)
      new.row <- c(gene, mean(gene.p.vec), mean(gene.z.vec), gene.t.test$p.value)
      out.df[nrow(out.df) + 1,] = new.row     
    }
  }
  write.csv(out.df, out.name, row.names = FALSE)
}


compileData <- function(in.dir, out.base.name){
  setwd(in.dir)
  files <- list.files(path=in.dir, pattern="Amt_*", full.name=FALSE, recursive=FALSE)
  newheader <- simpName(files)
  newheader <- c("gene_name", newheader)
  df.1 <- read.csv(files[1])
  out.df.p <- data.frame(matrix(ncol = 0, nrow = nrow(df.1)))
  out.df.z <- data.frame(matrix(ncol = 0, nrow = nrow(df.1)))
  out.df.t <- data.frame(matrix(ncol = 0, nrow = nrow(df.1)))
  out.df.p$gene_name <- as.character(unlist(df.1$gene_name))
  out.df.z$gene_name <- as.character(unlist(df.1$gene_name))
  out.df.t$gene_name <- as.character(unlist(df.1$gene_name))
  
  for (file in files){
    print(file)
    df.x <- read.csv(file)
    df.p <- data.frame("gene_name" = df.x$gene_name)
    df.z <- data.frame("gene_name" = df.x$gene_name)
    df.t <- data.frame("gene_name" = df.x$gene_name)
    
    df.p$pval <- as.numeric(unlist(df.x$avg_percentile))
    df.z <- cbind(df.z, as.numeric(unlist(df.x$avg_z.score)))
    df.t$t <- as.numeric(unlist(df.x$t.test_p.value))
   
    colnames(df.p) <- c("gene_name", simpNameOne(file))
    colnames(df.z) <- c("gene_name", simpNameOne(file))
    colnames(df.t) <- c("gene_name", simpNameOne(file))
    
    out.df.p <- merge(out.df.p, df.p, by="gene_name")
    out.df.z <- merge(out.df.z, df.z, by="gene_name")
    out.df.t <- merge(out.df.t, df.t, by="gene_name")
  }
  
  colnames(out.df.p) <- newheader
  colnames(out.df.z) <- newheader
  colnames(out.df.t) <- newheader
  
  out.name.p <- paste(out.base.name, "_avg_percentile.csv", sep="")
  out.name.z <- paste(out.base.name, "avg_z-score.csv", sep="")
  out.name.t <- paste(out.base.name, "t-test_p-value.csv", sep="")
  
  write.csv(out.df.p, out.name.p, row.names = FALSE)
  write.csv(out.df.z, out.name.z, row.names = FALSE)
  write.csv(out.df.t, out.name.t, row.names = FALSE)
  
}


########## Main ##########
###----- Calculate p-value of each gene for different datasets
#--- Amt
if (FALSE){
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/1_2_normtocontrol_ZP/Amt"
  setwd(wk.dir)
  files <- list.files(path=wk.dir, pattern="*ZP.csv", full.name=FALSE, recursive=FALSE)
  for (file in files) {
    tByGene(file)
  }
}

#--- Amt GFP
if (FALSE){
  wk.dir <-"/Volumes/Yolanda/CRF_Screen/InVitro/1_2_normtocontrol_ZP/AmtGFP"
  setwd(wk.dir)
  files <- list.files(path=wk.dir, pattern="*ZP.csv", full.name=FALSE, recursive=FALSE)
  for (file in files) {
    tByGene(file)
  }
}

###----- Compile data
if (FALSE){
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/Sep_Amt"
  out.base.name <- "Amt_normbycontrolZP_t-test.by.gene"
  compileData(wk.dir, out.base.name)
  
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/2_0_t-test_by_gene/Sep_AmtGFP"
  out.base.name <- "AmtGFP_normbycontrolZP_t-test.by.gene"
  compileData(wk.dir, out.base.name)
}


