########## Data Conversion for CRF library in vitro screening results ##########
# Author: Huitian (Yolanda) Diao
# Dec 3rd, 2018
# For data normalization and staticial analysis

########## Libraries ##########

########## Self-defined functions ##########
dataStats <- function(file.x){
  df.x <- read.csv(file.x)
  num.x.vec <- as.numeric(unlist(df.x[,4]))
  num.x.vec <- na.omit(num.x.vec)
  return(c(mean(num.x.vec), sd(num.x.vec)))
}

simpName <- function(vec.x) {
  new.vec <- character(length(vec.x))
  for (x in 1:length(vec.x)){
    new.x <- gsub("_shRNA_normZP.csv", "", vec.x[x])
    new.x <- gsub("_shRNA_normbycontrolZP.csv", "", new.x)
    new.x <- gsub("Amt_", "", new.x)
    new.x <- gsub("AmtGFP_", "", new.x)
    new.vec[x] <- new.x
  }
  return(new.vec)
}

simpNameOne <- function(str.x) {
  new.x <- gsub("_shRNA_normZP.csv", "", str.x)
  new.x <- gsub("_shRNA_normbycontrolZP.csv", "", new.x)
  new.x <- gsub("Amt_", "", new.x)
  new.x <- gsub("AmtGFP_", "", new.x)
  return(new.x)
}

dataStatsControl <- function(file.x, cal.col){
  #file.x <- "Amt_CD25GeoMean_10U_shRNA.csv"
  df.x <- read.csv(file.x)
  df.c <- subset(df.x, grepl("control",df.x$Position))
  df.c1 <- subset(df.c, df.c$Pool==1)
  df.c2 <- subset(df.c, df.c$Pool==2)
  df.c3 <- subset(df.c, df.c$Pool==3)
  num.c1.vec <- as.numeric(unlist(df.c1[,cal.col]))
  num.c1.vec <- na.omit(num.c1.vec)
  num.c2.vec <- as.numeric(unlist(df.c2[,cal.col]))
  num.c2.vec <- na.omit(num.c2.vec)
  num.c3.vec <- as.numeric(unlist(df.c3[,cal.col]))
  num.c3.vec <- na.omit(num.c3.vec)
  return(c(mean(num.c1.vec), mean(num.c2.vec), mean(num.c3.vec),
           sd(num.c1.vec), sd(num.c2.vec), sd(num.c3.vec)))
}

dataStatsByPool <- function(file.x, cal.col){
  #file.x <- "Amt_CD25GeoMean_10U_shRNA.csv"
  df.x <- read.csv(file.x)
  df.c1 <- subset(df.x, df.x$Pool==1)
  df.c2 <- subset(df.x, df.x$Pool==2)
  df.c3 <- subset(df.x, df.x$Pool==3)
  num.c1.vec <- as.numeric(unlist(df.c1[,cal.col]))
  num.c1.vec <- na.omit(num.c1.vec)
  num.c2.vec <- as.numeric(unlist(df.c2[,cal.col]))
  num.c2.vec <- na.omit(num.c2.vec)
  num.c3.vec <- as.numeric(unlist(df.c3[,cal.col]))
  num.c3.vec <- na.omit(num.c3.vec)
  return(c(mean(num.c1.vec), mean(num.c2.vec), mean(num.c3.vec),
           sd(num.c1.vec), sd(num.c2.vec), sd(num.c3.vec)))
}

zpByPool <- function(file.x){
  #file.x <- "Amt_CD25GeoMean_10U_shRNA.csv"
  new.name <- gsub(".csv", "_normZP.csv", file.x)
  df.x <- read.csv(file.x)
  df.x <- na.omit(df.x)
  df.x <- subset(df.x, shRNA != "Empty")
  df.c1 <- subset(df.x, df.x$Pool==1)
  df.c2 <- subset(df.x, df.x$Pool==2)
  df.c3 <- subset(df.x, df.x$Pool==3)
  num.c1.vec <- as.numeric(unlist(df.c1[,4]))
  num.c2.vec <- as.numeric(unlist(df.c2[,4]))
  num.c3.vec <- as.numeric(unlist(df.c3[,4]))
  c1.p <- (pnorm(num.c1.vec, mean=mean(num.c1.vec), sd=sd(num.c1.vec)))*100
  c1.z <- as.numeric(unlist(scale(num.c1.vec)))
  df.c1$precentile <- c1.p
  df.c1$z <- c1.z
  head(df.c1)
  c2.p <- (pnorm(num.c2.vec, mean=mean(num.c2.vec), sd=sd(num.c2.vec)))*100
  c2.z <- as.numeric(unlist(scale(num.c2.vec)))
  df.c2$precentile <- c2.p
  df.c2$z <- c2.z
  head(df.c2)
  c3.p <- (pnorm(num.c3.vec, mean=mean(num.c3.vec), sd=sd(num.c3.vec)))*100
  c3.z <- as.numeric(unlist(scale(num.c3.vec)))
  df.c3$precentile <- c3.p
  df.c3$z <- c3.z
  head(df.c3)
  new.df <- rbind(df.c1, df.c2, df.c3)
  write.csv(new.df, new.name, row.names = FALSE)
}

zpByPooltoControl <- function(file.x){
  #file.x <- "Amt_CD25GeoMean_10U_shRNA.csv"
  new.name <- gsub(".csv", "_normbycontrolZP.csv", file.x)
  df.x <- read.csv(file.x)
  df.x <- na.omit(df.x)
  df.x <- subset(df.x, shRNA != "Empty")
  df.c <- subset(df.x, grepl("control",df.x$Position))
  df.c1 <- subset(df.c, df.c$Pool==1)
  df.c2 <- subset(df.c, df.c$Pool==2)
  df.c3 <- subset(df.c, df.c$Pool==3)
  num.c1.vec <- as.numeric(unlist(df.c1[,4]))
  num.c2.vec <- as.numeric(unlist(df.c2[,4]))
  num.c3.vec <- as.numeric(unlist(df.c3[,4]))
  
  df.g1 <- subset(df.x, df.x$Pool==1)
  df.g2 <- subset(df.x, df.x$Pool==2)
  df.g3 <- subset(df.x, df.x$Pool==3)
  num.g1.vec <- as.numeric(unlist(df.g1[,4]))
  num.g2.vec <- as.numeric(unlist(df.g2[,4]))
  num.g3.vec <- as.numeric(unlist(df.g3[,4]))
  g1.p <- (pnorm(num.g1.vec, mean=mean(num.c1.vec), sd=sd(num.c1.vec)))*100
  g1.z <- (num.g1.vec - mean(num.c1.vec))/sd(num.c1.vec)
  df.g1$precentile <- g1.p
  df.g1$z <- g1.z
  g2.p <- (pnorm(num.g2.vec, mean=mean(num.c2.vec), sd=sd(num.c2.vec)))*100
  g2.z <- (num.g2.vec - mean(num.c2.vec))/sd(num.c2.vec)
  df.g2$precentile <- g2.p
  df.g2$z <- g2.z
  g3.p <- (pnorm(num.g3.vec, mean=mean(num.c3.vec), sd=sd(num.c3.vec)))*100
  g3.z <- (num.g3.vec - mean(num.c3.vec))/sd(num.c3.vec)
  df.g3$precentile <- g3.p
  df.g3$z <- g3.z
  new.df <- rbind(df.g1, df.g2, df.g3)
  
  write.csv(new.df, new.name, row.names = FALSE)
}

compileData <- function(in.dir, out.base.name){
  setwd(in.dir)
  ref.file <- "/Volumes/Yolanda/CRF_Screen/InVitro/position_Ref.csv"
  ref.df <- read.csv(ref.file)

  files <- list.files(path=in.dir, pattern="Amt_*", full.name=FALSE, recursive=FALSE)
  newheader <- simpName(files)
  newheader
  newheader <- c("Position", "shRNA", newheader)
  
  newheader
  out.df.p <- ref.df
  out.df.z <- ref.df
  
  for (file in files){
    #file <- files[1]
    df.x <- read.csv(file)
    df.p <- data.frame("plate_well" = df.x$Position)
    df.z <- data.frame("plate_well" = df.x$Position)
    
    df.p <- cbind(df.p, as.numeric(unlist(df.x$precentile)))
    df.z <- cbind(df.z, as.numeric(unlist(df.x$z)))
    
    colnames(df.p) <- c("plate_well", simpNameOne(file))
    colnames(df.z) <- c("plate_well", simpNameOne(file))
    
    out.df.p <- merge(out.df.p, df.p, by="plate_well")
    out.df.z <- merge(out.df.z, df.z, by="plate_well")
  }
  
  colnames(out.df.p) <- newheader
  colnames(out.df.z) <- newheader
  
  out.name.p <- paste(out.base.name, "_percentile.csv", sep="")
  out.name.z <- paste(out.base.name, "_z-score.csv", sep="")
  
  write.csv(out.df.p, out.name.p, row.names = FALSE)
  write.csv(out.df.z, out.name.z, row.names = FALSE)
}

########## Main ##########
###----- Create a data stats spreadsheet
if(FALSE){
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/1_1_shRNAmatched"
  setwd(wk.dir)
  files <- list.files(path=wk.dir, pattern="*.csv", full.names=FALSE, recursive=FALSE)
  
  stats.df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(stats.df) <- c("dataSet", "mean", "sd")
  
  for (file in files){
    new.row <- c(file, dataStats(file))
    stats.df[nrow(stats.df) + 1,] = new.row
  }
  write.csv(stats.df, "dataStats.csv", row.names = FALSE)
}

###----- Create a data stats spreadsheet for the control plates
if(FALSE){
  wk.dir <-  "/Volumes/Yolanda/CRF_Screen/InVitro/1_1_shRNAmatched"
  setwd(wk.dir)
  files <- list.files(path=wk.dir, pattern="*.csv", full.names=FALSE, recursive=FALSE)
  
  stats.df <- data.frame(matrix(ncol = 7, nrow = 0))
  colnames(stats.df) <- c("dataSet", "mean_control1", "mean_control2", "mean_control3",
                          "sd_control1", "sd_control2", "sd_control3")
  for (file in files){
    new.row <- c(file, dataStatsControl(file))
    stats.df[nrow(stats.df) + 1,] = new.row
  }
  write.csv(stats.df, "dataStats_controlPlates.csv", row.names = FALSE)
}

###----- Create a data stats spreadsheet by pool
if(FALSE){
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/1_1_shRNAmatched"
  setwd(wk.dir)
  files <- list.files(path=wk.dir, pattern="*.csv", full.names=FALSE, recursive=FALSE)
  
  stats.df <- data.frame(matrix(ncol = 7, nrow = 0))
  colnames(stats.df) <- c("dataSet", "mean_pool1", "mean_pool2", "mean_pool3",
                          "sd_pool1", "sd_pool2", "sd_pool3")
  for (file in files){
    new.row <- c(file, dataStatsByPool(file))
    stats.df[nrow(stats.df) + 1,] = new.row
  }
  write.csv(stats.df, "dataStats_byPool.csv", row.names = FALSE)
}

###----- Calculate percentile and Z-score by pool
if(FALSE){
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/1_1_shRNAmatched/Amt"
  setwd(wk.dir)
  files <- list.files(path=wk.dir, pattern="*.csv", full.names=FALSE, recursive=FALSE)
  for (file in files){
    zpByPool(file)
  }
  
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/1_1_shRNAmatched/AmtGFP"
  setwd(wk.dir)
  files <- list.files(path=wk.dir, pattern="*.csv", full.names=FALSE, recursive=FALSE)
  for (file in files){
    zpByPool(file)
  }
}

###----- Calculate percentile and Z-score by control for each pool
if(FALSE){
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/1_1_shRNAmatched/Amt"
  setwd(wk.dir)
  files <- list.files(path=wk.dir, pattern=".csv", full.names=FALSE, recursive=FALSE)
  for (file in files){
    zpByPooltoControl(file)
  }
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/1_1_shRNAmatched/AmtGFP"
  setwd(wk.dir)
  files <- list.files(path=wk.dir, pattern=".csv", full.names=FALSE, recursive=FALSE)
  for (file in files){
    zpByPooltoControl(file)
  }
}

###----- Compile percentile and z-score
if (FALSE){
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/1_2_normtocontrol_ZP/Amt"
  base.name <- "Amt_normbycontrolZP"
  compileData(wk.dir, base.name)
  
  wk.dir <- "/Volumes/Yolanda/CRF_Screen/InVitro/1_2_normtocontrol_ZP/AmtGFP"
  base.name <- "AmtGFP_normbycontrolZP"
  compileData(wk.dir, base.name)
}

