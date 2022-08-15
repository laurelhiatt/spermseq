setwd("/Users/quinlan/Documents/Quinlan-PhD/SpermSeq")
###wherever files are
list.files()
spermseq_csv = read.csv("GleesonSpermVariants.vcf", sep = '\t', header = TRUE)

spermseq_bed = spermseq_csv

spermseq_bed$start <- NA
spermseq_bed$end <- NA

spermseq_bed$start = (spermseq_csv$POS-1)
spermseq_bed$end = (spermseq_csv$POS)
library(dplyr)

spermseq_bed <- subset (spermseq_bed, select = -POS)
library(tidyverse)  
write_tsv(spermseq_bed, "GleesonSpermVariants.bed")