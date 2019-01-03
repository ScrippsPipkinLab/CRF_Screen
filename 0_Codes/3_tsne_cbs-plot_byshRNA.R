########## CBS-plot ##########
# Author: Huitian (Yolanda) Diao
# Dec 3rd, 2018
# For plottinng color-bar scatter plot

########## Libraries ##########
library(ggplot2)
library(ggrepel)

########## Self-defined functions ##########
shRNAtoGene <- function(str.x){
  str.x.vec <- strsplit(str.x, "\\.")
  return(unlist(str.x.vec)[1])
}

sltAnnoCols <- function(df.x, gene.vec){
  new.vec <- character(nrow(df.x))
  for (x in 1:nrow(df.x)){
    gene.x <- shRNAtoGene(as.character(unlist(df.x$shRNA))[x])
    if(gene.x %in% gene.vec){
      new.vec[x] <- gene.x
      #print(gene.x)
    }
    else{
      new.vec[x] <- "other"
    }
  }
  return(new.vec)
}

sltAnnoGenes <- function(df.x, gene.vec){
  new.vec <- character(nrow(df.x))
  for (x in 1:nrow(df.x)){
    gene.x <- shRNAtoGene(as.character(unlist(df.x$shRNA))[x])
    if(gene.x %in% gene.vec){
      new.vec[x] <- as.character(unlist(df.x$shRNA))[x]
      #print(gene.x)
    }
    else{
      new.vec[x] <- ""
    }
  }
  return(new.vec)
}



cbs_plot <- function(in.file, anno_genes, anno_cols){
  #in.file <- "Amt_normbycontrolZP_z-score_tsne_per10.csv"
  #anno_genes <- c("A1l")
  #anno_cols <- c("red", "gray80")
  in.df <- read.csv(in.file)
  in.df$col <- sltAnnoCols(in.df, anno_genes)
  in.df$anno <- sltAnnoGenes(in.df, anno_genes)
  
  cbs.plot <- ggplot(in.df, aes(x=in.df$t1, y=in.df$t2, color=in.df$col)) +
    geom_point(alpha = 0.6, shape = 19,  size = 2) +
    scale_x_continuous(limits = c(-100, 100)) +
    scale_y_continuous(limits = c(-100, 100)) +
    scale_color_manual(values = anno_cols) +
    theme(panel.background = element_rect(fill = "white", colour = "white",size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "grey"), 
          panel.grid.minor = element_line(size = 0.15, linetype = 'solid',colour = "grey"),
          text = element_text(size=10))
  
  out.name <- paste(anno_genes, collapse="-")
  out.name <- paste(out.name, in.file, sep="--")
  out.name <- gsub(".csv", ".png", out.name)
  
  ggsave(out.name, plot = cbs.plot, device = "png", units = "cm", width = 15, height = 12)
}

########## Main ##########

###----- Build dataframes
###----- Plot
crf.genes <- c("A1l", "A2l", "Actl6b", "Aff1", "Aff4", "Arid1a", "Arid1b", "Arid2", "Asxl1", "Asxl2", 
               "Asxl3", "Atad2", "Atad2b", "Atm", "Atr", "Atrx", "Aurkb", "Baz1a", "Baz1b", "Baz2a", 
               "Baz2b", "Bmi1", "Bop1", "Bptf", "Brd1", "Brd2", "Brd3", "Brd4", "Brd7", "Brd8", "Brd9", 
               "Brdt", "Brpf1", "Brpf3", "Brwd1", "Brwd3", "Carm1", "Cbx1", "Cbx2", "Cbx3", "Cbx4", 
               "Cbx5", "Cbx6", "Cbx7", "Cbx8", "Cd4", "Cecr2", "Chaf1a", "Chaf1b", "Chd1", "Chd1l", 
               "Chd2", "Chd3", "Chd4", "Chd5", "Chd6", "Chd7", "Chd8", "Chd9", "Clock", "Crebbp", 
               "Ctsl", "Cxxc1", "Cxxcl", "Dmap1", "Dnmt1", "Dnmt3a", "Dnmt3b", "Dnmt3l", "Dot1l", 
               "Dpy30", "Dub2a", "Eaf1", "Eed", "Ehmt1", "Ehmt2", "Elp3", "Elp4", "Empty", "Ep300", 
               "Ep400", "Epc1", "Epc2", "Ercc5", "Ezh1", "Ezh2", "Fbxl19", "Fbxo17", "Fbxo44", "Fbxw9", 
               "Fkbp1a", "Fkbp2", "Fkbp5", "Gadd45a", "Gadd45b", "Gtf2b", "Gtf2f1", "Gtf2h1", "Gtf3c4", 
               "H2afz", "Hat1", "Hcfc1", "Hdac1", "Hdac10", "Hdac11", "Hdac2", "Hdac3", "Hdac4", "Hdac5", 
               "Hdac6", "Hdac7", "Hdac8", "Hdac9", "Hells", "Hira", "Hltf", "Ing2", "Ing3", "Ing4", "Ing5", 
               "Ino80", "Iws1", "Jarid2", "Jhdm1d", "Jmjd1c", "Jmjd4", "Jmjd5", "Jmjd6", "Jmjd8", "Kat2a", 
               "Kat2b", "Kat5", "Kdm1a", "Kdm1b", "Kdm2a", "Kdm2b", "Kdm3a", "Kdm3b", "Kdm4a", "Kdm4b", "Kdm4c", 
               "Kdm4d", "Kdm5a", "Kdm5b", "Kdm5c", "Kdm5d", "Kdm6a", "Kdm6b", "Klf2", "L3mbtl1", "L3mbtl2", "L3mbtl3", 
               "L3mbtl4", "Mbd1", "Mbd2", "Mbd3", "Mbd4", "Mecom", "Mecp2", "Men1", "Mll1", "Mll2", "Mll3", "Mll5", 
               "Morf4l1", "Msl3", "Mta1", "Mta2", "Mta3", "Mtf2", "Myst1", "Myst2", "Myst3", "Myst4", "Nap1l1", 
               "Nap1l2", "Nap1l3", "Ncoa1", "Ncoa3", "Ncor1", "Ncor2", "Nsd1", "Orc1", "Padi1", "Padi2", "Padi4", 
               "Padi6", "Paf1", "Parp1", "Parp2", "Paxip1", "Pbrm1", "Pcgf1", "Pcgf2", "Pcgf5", "Pcgf6", "Pcmt1", 
               "Phc1", "Phc2", "Phc3", "Phf1", "Phf10", "Phf17", "Phf19", "Phf2", "Phf20", "Phf20l1", "Phf8", 
               "Phip", "Polr2b", "Ppargc1a", "Prdm1", "Prdm10", "Prdm11", "Prdm12", "Prdm13", "Prdm14", 
               "Prdm15", "Prdm16", "Prdm2", "Prdm4", "Prdm5", "Prdm6", "Prdm8", "Prdm9", "Prkaa1", "Prkaa2", 
               "Prkcd", "Prmt1", "Prmt2", "Prmt3", "Prmt5", "Prmt6", "Prmt7", "Prmt8", "Psip1", "Rbbp4", 
               "Rbbp5", "Rbbp7", "Ring1", "Rnf2", "Rnf20", "Rnf217", "Rnf40", "Rpa3", "Runx3", "Satb1", 
               "Scmh1", "Scml2", "Scml4", "Setd1a", "Setd2", "Setd3", "Setd4", "Setd5", "Setd7", "Setd8", 
               "Setdb1", "Setdb2", "Setmar", "Sfmbt1", "Sfmbt2", "Sin3a", "Sin3b", "Sirt1", "Sirt3", "Sirt4", 
               "Sirt5", "Sirt6", "Sirt7", "Smarca1", "Smarca2", "Smarca4", "Smarca5", "Smarcb1", "Smarcc1", 
               "Smarcc2", "Smarcd1", "Smarcd2", "Smarcd3", "Smarce1", "Smyd1", "Smyd2", "Smyd3", "Smyd4", 
               "Smyd5", "Sp100", "Sp110", "Sp140", "Ssrp1", "Supt16h", "Suv39h1", "Suv39h2", "Suv420h1", 
               "Suv420h2", "Suz12", "Taf1", "Taf3", "Tcea1", "Tdg", "Tet1", "Tet2", "Tet3", "Trim24", "Trim28", 
               "Trim33", "Trim66", "Ube2a", "Ube2b", "Ube2e1", "Ube2i", "Uhrf1", "Usp22", "Usp27x", "Usp51", 
               "Wbp7", "Wdr5", "Wdr82", "Whsc1", "Whsc1l1", "Zmynd11", "Zmynd8")
if (FALSE){
  wk.dir <- "/Users/yolandatiao/Desktop/CRF_Screen/1_2_normtocontrol_ZP/Compiled_Amt"
  setwd(wk.dir) 
  tsne.per10 <- "Amt_normbycontrolZP_z-score_tsne_per10.csv"
  
  for (gene in crf.genes) {
    gene.vec <- c(gene)
    col.vec <- c("red", "gray80")
    cbs_plot(tsne.per10, gene.vec, col.vec)
  }
  
  tsne.per5 <- "Amt_normbycontrolZP_z-score_tsne_per5.csv"
  
  for (gene in crf.genes) {
    gene.vec <- c(gene)
    col.vec <- c("red", "gray80")
    cbs_plot(tsne.per5, gene.vec, col.vec)
  }
  
}


