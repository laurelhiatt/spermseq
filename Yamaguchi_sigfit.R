library(usethis) 
usethis::edit_r_environ()


# https://htmlpreview.github.io/?https://github.com/HirofumiNakaoka/endometrium_natcommun_2021/blob/main/code_source_data_04-05.html

install.packages('barplot3d')
library(barplot3d)
library(ggplot2)
library(rgl)
#devtools::install_github("kgori/sigfit", build_vignettes = TRUE,
 #                        build_opts = c("--no-resave-data", "--no-manual"))

library(sigfit)

data("cosmic_signatures_v3.2")

blood_count <- read.csv('/Users/quinlan/Documents/Quinlan-PhD/SpermSeq/BloodMutationCountMatrix.csv')
sperm_count <- read.csv('/Users/quinlan/Documents/Quinlan-PhD/SpermSeq/SpermMutationCountMatrix.csv')


data("counts_21breast")

bloodswitched <- t(as.data.frame(blood_count[-1]))
spermswitched <- t(as.data.frame(sperm_count[-1]))
colnames(bloodswitched) <- colnames(counts_21breast)
colnames(spermswitched) <- colnames(counts_21breast)


mcmc_samples_fit <- fit_signatures(
  counts=bloodswitched,
  signatures=cosmic_signatures_v3.2,
  iter=50000,
  warmup=25000,
  chains=4,
  seed=1896)

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

library(cowplot)
ggplot(meltsperm, aes(variable, value,  fill = variable)) + 
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
