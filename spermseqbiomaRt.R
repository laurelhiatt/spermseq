library(BiocManager)
library(biomaRt)

# http://127.0.0.1:28249/library/biomaRt/doc/accessing_ensembl.html#step1-identifying-the-database-you-need

###library(EnsDb.Hsapiens.v79) This is the option below; recommended on stackoverflow as most accurate
### If interested, separate syntax is necessary.
getwd()
setwd("/Users/quinlan/Documents/Quinlan-PhD/SpermSeq/spermseq") ###or wherever you want to be
list.files() ###make sure the file you want is where you want it

data <- read.csv("appris_gtex_gtf_pfam_uniprot_3dhotspot_chang.txt", stringsAsFactors = TRUE, header = TRUE, sep='\t')
# import our data 

#ensembl.genes <- data$gene_id ##filter by whatever row has your ensembl info

keeps <- c("chr","start", "stop", "gene_name")
df = data[keeps]
# i want to sort by the loci and gene_name so i'm subsetting the data to that

#ensembl.transcript <- data$transcript_id # 


martthatworks <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

listEnsembl()

datasets <- listDatasets(ensembl)
# head(datasets)

ensembl <- useEnsembl(biomart="snps", dataset="hsapiens_snp")
# i try several of these... 
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_SNP", 
                      dataset = "hsapiens_snp", 
                      mirror = "useast")

listMarts()
mart <- useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")


filters <- listFilters(ensembl)
# here we can get the filters that form our input 

library(dplyr)
df <- df %>%
  rename(chr_name = chr,
         end = stop,
         ensemble_gene = gene_name)

attributes <- listAttributes(ensembl)
# here we can get the attributes we're interested in

snpmart = useEnsembl(biomart = "snp", dataset="hsapiens_snp")
snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
## these are examples from stackoverflow to try and make things work 

snp_mart = useMart(biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp")

data_list <- getBM(attributes = c('chr_name', 'chrom_start', 'chrom_end', 'ensembl_gene_name', 'pmid'),
                   filters = c('chr_name', 'start', 'end', 'ensembl_gene'),
                   values = list(df=c("chr_name", "start", "end", "ensembl_gene"), 
                   mart = ensembl))

data_list <- getBM(attributes = c('chr_name', 'chrom_start', 'chrom_end', 'ensembl_gene_name', 'pmid'),
                   filters = c('chr_name', 'start', 'end', 'ensembl_gene'),
                   values = list(df=c("chr_name", "start", "end", "ensembl_gene"), 
                                 mart = snpmart))

data_list <- getBM(attributes = c('chr_name', 'chrom_start', 'chrom_end', 'ensembl_gene_name', 'pmid'),
                   filters = c('chr_name', 'start', 'end', 'ensembl_gene'),
                   values = list(df=c("chr_name", "start", "end", "ensembl_gene"), 
                                 mart = snp_mart))

data_list <- getBM(attributes = c('chr_name', 'chrom_start', 'chrom_end', 'ensembl_gene_name', 'pmid'),
                   filters ='ensembl_gene',
                   values = list(df$ensemble_gene, 
                                 mart = mart))

?useMart

write.table(G_list, file='nameoffile.tsv', quote=FALSE, index = FALSE, sep='\t')

