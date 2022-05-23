setwd("/Users/quinlan/Documents/Quinlan-PhD/SpermSeq")
###wherever files are
list.files()
spermseq_csv = read.csv("shared_blood_sperm_mutation_analysis.csv", sep = ',', header = TRUE)

spermseq_bed = spermseq_csv
spermseq_bed$start = (spermseq_csv$start-1)
spermseq_bed$end = (spermseq_csv$end-1)


  
write_tsv(spermseq_bed, path = "shared_blood_sperm_mutation_analysis.bed")