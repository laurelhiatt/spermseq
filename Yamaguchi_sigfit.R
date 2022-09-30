library(usethis) 
usethis::edit_r_environ()


# https://htmlpreview.github.io/?https://github.com/HirofumiNakaoka/endometrium_natcommun_2021/blob/main/code_source_data_04-05.html

#install.packages('barplot3d')
library(barplot3d)
library(ggplot2)
library(rgl)
#devtools::install_github("kgori/sigfit", build_vignettes = TRUE,
 #                        build_opts = c("--no-resave-data", "--no-manual"))

library(sigfit)

data("cosmic_signatures_v3.2")

# here are the Ceph trinucleotide counts by generation quantified
secondceph <- read.csv('/Users/quinlan/Documents/Quinlan-PhD/SpermSeq/spermseq/Ceph2ndGenTriCollapsedPaternal.csv', header = FALSE, sep = ' ')
thirdceph <- read.csv('/Users/quinlan/Documents/Quinlan-PhD/SpermSeq/spermseq/Ceph3rdGenTriCollapsedPaternal.csv', header = FALSE, sep = ' ')

#do a switcheroo
secondcephswitched <- t(as.data.frame(secondceph[1]))
thirdcephswitched <- t(as.data.frame(thirdceph[1]))

#get column names
colnames(secondcephswitched) <- colnames(cosmic_signatures_v3.2)
colnames(thirdcephswitched) <- colnames(cosmic_signatures_v3.2)

#create combined file
combinedcephswitched <- secondcephswitched + thirdcephswitched

#get proportions
secondcephswitched<- secondcephswitched/rowSums(secondcephswitched)
thirdcephswitched<- thirdcephswitched/rowSums(thirdcephswitched)
combinedcephswitched<-combinedcephswitched/rowSums(combinedcephswitched)


blood_count <- read.csv('/Users/quinlan/Documents/Quinlan-PhD/SpermSeq/spermseq/BloodMutationCountMatrix.csv')
sperm_count <- read.csv('/Users/quinlan/Documents/Quinlan-PhD/SpermSeq/spermseq/SpermMutationCountMatrix.csv')

bloodswitched <- t(as.data.frame(blood_count[-1]))
spermswitched <- t(as.data.frame(sperm_count[-1]))
colnames(bloodswitched) <- colnames(cosmic_signatures_v3.2)
colnames(spermswitched) <- colnames(cosmic_signatures_v3.2)

plot_spectrum(secondcephswitched)
plot_spectrum(thirdcephswitched)
plot_spectrum(combinedcephswitched)


### FOR COSMIC STUFF!!
mcmc_samples_fit <- fit_extract_signatures(
  counts=spermswitched,
  signatures=secondcephswitched,
  num_extra_sigs = 1,
  iter=50000,
  warmup=25000,
  chains=1,
  seed=1896)

?fit_extract_signatures

extr_sigs <- retrieve_pars(mcmc_samples_fit, "signatures")
plot_spectrum(extr_sigs)

extr_sigs <- retrieve_pars(mcmc_samples, "signatures")
plot_spectrum(cosmic_signatures_v2[7, ], pdf_path = "COSMIC_Sig7.pdf", name="COSMIC sig. 7")
plot_spectrum(extr_sigs, pdf_path = "Extracted_Sigs.pdf")

### this is from before, fitting to cosmic
mcmc_samples_fit <- fit_signatures(
  counts=spermswitched,
  signatures=cosmic_signatures_v3.2,
  iter=50000,
  warmup=25000,
  chains=4,
  seed=1896)

exposures <- retrieve_pars(
  mcmc_samples_fit,
  par="exposures",
  hpd_prob=0.90)

data <- do.call(rbind.data.frame, exposures)

meandata <- head(data,24)

library(reshape2)
library(tibble)
df <- tibble::rownames_to_column(meandata, "sample")

meltsperm <- melt(df, "sample")

# for the fit_extract
library(cowplot)
ggplot(meltsperm, aes(sample, value,  fill = variable)) + 
  geom_bar(stat="identity", show_guide=FALSE)  + 
  theme_cowplot(12) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.box.background = element_rect(fill='#fffafa'),
        legend.background = element_blank()) + 
  xlab('Sperm Signatures') + ylab('Proportion')


ggplot(meltblood, aes(variable, value,  fill = variable)) + 
  facet_wrap(~ sample, scales = 'free') + geom_bar(stat="identity", show_guide=FALSE)  + 
  theme_cowplot(12) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.box.background = element_rect(fill='#fffafa'),
        legend.background = element_blank()) + 
  xlab('Blood Signatures') + ylab('Proportion')
