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

names.genes.de <- c("PLCB1","MCL1","PRDX4","TTF2","TACC3", "PARP4","LSM1")
tmp.biogrid <- data.frame("Official.Symbol.Interactor.A" = names.genes.de,
                          "Official.Symbol.Interactor.B" = rev(names.genes.de))
net.biogrid.de <- getAdjacencyBiogrid(tmp.biogrid, names.genes.de)
net.biogrid.de