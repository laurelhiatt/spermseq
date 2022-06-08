#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

### First, load necessary libraries. I've put installation code that should work
### if that is necessary, not sure what to expect in the protected environment.

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
# this is the library to access the reference genome

#BiocManager::install("MutationalPatterns")
library(MutationalPatterns)
# this actually generates the mut count matrix

# install.packages("dplyr")
library(dplyr)
# this is to convert imported files to a list, which is necessary for input
# for the read_vcfs_as_granges
# there may be a "faster" way to do this but I think this is the most concise

# setwd('/Users/quinlan/Documents/Quinlan-PhD/SpermSeq/spermseq')
# set working direcory
# this should be whatever working directory you want to set
# not necessary if your arguments are full file paths, but it's an option

## Get a reference genome BSgenome object.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"

###sample_names <- read.csv('sample_names.txt', header = FALSE, sep = ',')
sample_names <- read.csv(args[1], header = FALSE, sep = ',')
# I would make this a csv file, ignore the txt
# you may need to play with whatever file type you want to put the sample_names
# as, but I can help if we need to change the separator.
# header = FALSE indicates no header, change to TRUE if there's a header

sample_names <- dplyr::pull(sample_names, V1)
# this converts the data frame of the imported file into a vector
# pulls column 1 (V1) of the sample_names df to a singular vector
# this is less code than importing the file directly into a list/vector
# feel free to play with it if you hate it but this should work

## We assemble a list of files we want to load. These files match the
## sample names defined above.
### vcf_files <- read.csv('vcf_file_names.txt', header = FALSE, sep = ',')
vcf_files <- read.csv(args[2], header = FALSE, sep = ',')

vcf_files <- dplyr::pull(vcf_files, V1)

## This function loads the files as GRanges objects.
## For backwards compatability reasons it only loads SNVs by default
#vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)

## Here I have snv as explicit because I think it's good to know: 
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome, type = "snv")

## Loading only indels can be done like this.
# vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome, type = "indel")

### Generating the final matrix!
mut_mat <- mut_matrix(vcfs, ref_genome)


###write.csv(mut_mat,"/Users/quinlan/u1264408/mutation_signatures/19610testMutationCountMatrix.csv", row.names = TRUE)
write.csv(mut_mat, args[3], row.names = TRUE)
#output the matrix wherever you want so that I can use it!
