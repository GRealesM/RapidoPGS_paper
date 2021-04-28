# Preprocessing datasets
# 2021/02/22

# This script will serve to process datasets (Asthma, RA, T1D-noHLA, T1D, T2D, BRCA, PRCA, MDD, CAD, BMI and Height).
# Also, we'll apply QC recommended by Florian Privé in previous versions.


library(data.table)
library(bigsnpr)
library(ggplot2)
setDTthreads(8)

## Helper function for quant traits

sdY.est <- function(vbeta, maf, n) {
  warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
  oneover <- 1/vbeta
  nvx <- 2 * n * maf * (1-maf)
  m <- lm(nvx ~ oneover - 1)
  cf <- coef(m)[['oneover']]
  if(cf < 0)
    stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
  return(sqrt(cf))
}

## QC variants
# For consistency across methods, we'll filter our summary statistic datasets by HapMap3 variants.
hm3 <- readRDS(url("https://ndownloader.figshare.com/files/25503788"))
hm3 <- hm3[,1:4]
# Additionally, we'll create another set of datasets, filtered by QC'd (ie. maf > 0.01 and hwe p-value > 1e-10) 1000 Genomes variants, for RápidoPGS comparison across different SNP sets.
varlist <- fread("../datasets/RapidoPGS_panel_variants.txt", header=FALSE)
varlist <- varlist[,V1]


## Asthma dataset

asthma <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/DemenaisF_29273806_GCST006862/harmonised/29273806-GCST006862-EFO_0000270-build37.f.tsv.gz",
	select = c(1,  3:8),
	col.names = c("CHR", "BP", "REF", "ALT", "BETA", "SE", "P"))# European dataset Demenais et al., 2018

asthma$n_eff <- 4 / (1 / 19954 + 1 / 107715)

# We need the frequencies
fwrite(asthma, "../datasets/Asthma-NOFREQS.tsv.gz", sep="\t")

# Let's compute the frequencies using an external tool and the RápidoPGS panel
system("Compute_freqs.sh -f ../datasets/Asthma-NOFREQS.tsv.gz -p CEU")

# Recover the file
asthma <- fread("../datasets/Asthma-withfreqs.tsv.gz")

## Apply QC recommended by LDpred2 authors.

sd_ss <- with(asthma, 2 / sqrt(n_eff * SE^2))
sd_val <- with(asthma, sqrt(2*ALT_FREQ*(1-ALT_FREQ)))

is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05 
table(is_bad)
# is_bad
#   FALSE    TRUE 
# 1987447    5685 

qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

ggsave("../figures/sd-approx-Asthma.png", width = 10, height = 7)

asthma <- asthma[!is_bad,]

asthma[,SNPID:=paste(CHR,BP,sep=":")]
# Filter by 1kG QC variants
invarlist <- asthma$SNPID %in% varlist
asthma.1kg <- asthma[invarlist,]

# Save file with 1kG variants
fwrite(asthma.1kg, "../datasets/Asthma-1kg.tsv.gz", sep="\t")

# Filter (and align) to hm3 variants manifest
setnames(asthma, c("CHR","BP","REF","ALT", "BETA"), c("chr","pos","a0","a1","beta"))
asthma.hm3 <- as.data.table(snp_match(asthma, hm3, strand_flip = FALSE))
#1,987,447 variants to be matched.
#0 ambiguous SNPs have been removed.
#932,300 variants have been matched; 0 were flipped and 0 were reversed.
setnames(asthma.hm3,c("chr","pos","a0","a1","beta"), c("CHR","BP","REF","ALT", "BETA"))
asthma.hm3[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]

# Save file with HapMap3 variants
fwrite(asthma.hm3, "../datasets/Asthma-hm3.tsv.gz", sep="\t")

rm(asthma, asthma.1kg, asthma.hm3, sd_ss, sd_val, is_bad, invarlist)


## RA (Rheumatoid Arthritis) dataset

# Now I chose to analyse the European-only meta-analysis instead. (RA_1 in our repo).
ra <- fread("http://plaza.umin.ac.jp/~yokada/datasource/files/GWASMetaResults/RA_GWASmeta_European_v2.txt.gz", 
	select= c(2:6,9),
	col.names = c("CHR", "BP", "ALT", "REF", "OR", "P"))
ra[,BETA:=log(OR)][,chi2:= qchisq(P, df = 1, lower.tail = FALSE)][,SE:= ifelse(chi2 > 1e-4, abs(BETA) / sqrt(chi2), NA)][, n_eff:= 4 / (1 / 14361 + 1 / 43923)] # Updated numbers here too


## This file didn't have MAF available, so we'll need to compute them using an external tool.
fwrite(ra, "../datasets/RA-NOFREQS.tsv.gz", sep="\t")
system("Compute_freqs.sh -f ../datasets/RA-NOFREQS.tsv.gz -p CEU")

## Recover the file
ra <- fread("../datasets/RA-withfreqs.tsv.gz")

## Apply QC recommended by LDpred2 authors.

sd_ss <- with(ra, 2 / sqrt(n_eff * SE^2))
sd_val <- with(ra, sqrt(2*ALT_FREQ*(1-ALT_FREQ)))

is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05 
table(is_bad)
# is_bad
#   FALSE    TRUE 
# 5443957 2617608 

qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

ggsave("../figures/sd-approx-RA.png", width = 10, height = 7)

ra <- ra[!is_bad,]

ra[,SNPID:=paste(CHR,BP,sep=":")]
ra[,c("OR","chi2"):=NULL]

# Filter by 1000G QC variants
invarlist <- ra$SNPID %in% varlist
ra.1kg <- ra[invarlist,]

# Write 1000G qc-filtered file
fwrite(ra.1kg, "../datasets/RA-1kg.tsv.gz", sep ="\t")

# Filter (and align) to hm3 variants manifest
setnames(ra, c("CHR","BP","REF","ALT", "BETA"), c("chr","pos","a0","a1","beta"))
ra.hm3 <- as.data.table(snp_match(ra, hm3, strand_flip=FALSE))
#5,443,957 variants to be matched.
#0 ambiguous SNPs have been removed.
#781,276 variants have been matched; 0 were flipped and 379,961 were reversed.
setnames(ra.hm3,c("chr","pos","a0","a1","beta"), c("CHR","BP","REF","ALT", "BETA"))
ra.hm3[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]

# Save file with HapMap3 variants
fwrite(ra.hm3, "../datasets/RA-hm3.tsv.gz", sep="\t")

# little test
#hm3f <- readRDS(url("https://ndownloader.figshare.com/files/25503788"))
#test <- as.data.table(snp_match(ra, hm3f))
#plot(test$ALT_FREQ, test$af_UKBB) 
# I wanted to test correlation between CEU allele freqs and UKBB. Despite the fact that half of SNPs are flipped (thus giving me a nice X plot), the correlation is quite high.
# Let's see what happens if we consider the MAF
#test$ALT_FREQ  <- ifelse(test$ALT_FREQ > 0.5, 1 - test$ALT_FREQ, test$ALT_FREQ)
#test$af_UKBB  <- ifelse(test$af_UKBB > 0.5, 1 - test$af_UKBB, test$af_UKBB)
#plot(test$ALT_FREQ, test$af_UKBB) 

rm(ra, ra.hm3, ra.1kg, sd_ss, sd_val, is_bad, invarlist)


## T1D (Type 1 diabetes) dataset

## We had access to the original Cooper dataset, which comes in one file.
t1d <- fread("../../../02-Processed/T1D_Cooper_doi101101120022_1-hg38.tsv.gz",
		select=c(2:6,16,17,19,28),
		col.names=c("CHR", "BP","REF", "ALT","ALT_FREQ","BETA","SE","P","qc.check"))

# NOTE: This file comes with allele frequencies, which we'll use in this case for QC.

## We apply QC check as in LDpred2 paper, including removing alleles with qc.check different than PASS, and those with ALT_FREQ < 1 / sqrt(n_eff)
t1d <- t1d[qc.check == "PASS",][,n_eff:=4 / (1 / 5913 + 1 / 8828)][ALT_FREQ > 1/sqrt(n_eff),][,qc.check:=NULL]
t1d <- na.omit(t1d)

## Apply QC recommended by LDpred2 authors.

sd_ss <- with(t1d, 2 / sqrt(n_eff * SE^2))
sd_val <- with(t1d, sqrt(2*ALT_FREQ*(1-ALT_FREQ)))

is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05 
table(is_bad)
# is_bad
#   FALSE    TRUE 
# 4182567    5023 

qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")
ggsave("../figures/sd-approx-T1D.png", width = 10, height = 7)

t1d <- t1d[!is_bad,]
t1d[,SNPID:=paste(CHR,BP,sep=":")]

# Filter by 1000G QC variants
invarlist <- t1d$SNPID %in% varlist
t1d.1kg <- t1d[invarlist,]

# Write 1000G qc-filtered file
fwrite(t1d.1kg, "../datasets/T1D-1kg.tsv.gz", sep ="\t")

# Filter (and align) to hm3 variants manifest
setnames(t1d, c("CHR","BP","REF","ALT", "BETA"), c("chr","pos","a0","a1","beta"))
t1d.hm3 <- as.data.table(snp_match(t1d, hm3, strand_flip=FALSE,  match.min.prop=0.4))
#4,182,567 variants to be matched.
#487,732 variants have been matched; 0 were flipped and 0 were reversed.
setnames(t1d.hm3,c("chr","pos","a0","a1","beta"), c("CHR","BP","REF","ALT", "BETA"))
t1d.hm3[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]

# Save file with HapMap3 variants
fwrite(t1d.hm3, "../datasets/T1D-hm3.tsv.gz", sep="\t")

rm(t1d, t1d.hm3, t1d.1kg, sd_ss,sd_val, is_bad, invarlist)


## T2D (Type 2 diabetes) dataset

t2d <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/ScottRA_28566273_GCST004773/METAANALYSIS_DIAGRAM_SE1.txt",
	col.names=c("SNPID","ALT","REF","BETA","SE","P","N"))
t2d[,c("CHR","BP"):=tstrsplit(SNPID, ":", fixed = TRUE)]
t2d[,n_eff:=4 / (1 / 26676 + 1 / 132532)][,CHR:=as.numeric(CHR)][,BP:=as.numeric(BP)][,N:=NULL]

# This file doesn't have ALT_FREQ, so we'll compute it using an external tool
fwrite(t2d, "../datasets/T2D-NOFREQS.tsv.gz", sep="\t")
system("Compute_freqs.sh -f ../datasets/T2D-NOFREQS.tsv.gz -p CEU")

# Let's recover the file
t2d <- fread("../datasets/T2D-withfreqs.tsv.gz")

## Apply QC recommended by LDpred2 authors.

sd_ss <- with(t2d, 2 / sqrt(n_eff * SE^2))
sd_val <- with(t2d, sqrt(2*ALT_FREQ*(1-ALT_FREQ)))

is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05 
table(is_bad)
# is_bad
#   FALSE    TRUE 
# 6550140 1514988 

qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")
ggsave("../figures/sd-approx-T2D.png", width = 10, height = 7)

t2d <- t2d[!is_bad,]
# Reordering
t2d <- t2d[,c("SNPID","CHR","BP", "REF","ALT","BETA","SE","P","n_eff")]

# Filter by 1000G QC variants
invarlist <- t2d$SNPID %in% varlist
t2d.1kg <- t2d[invarlist,]

# Write 1000G qc-filtered file
fwrite(t2d.1kg, "../datasets/T2D-1kg.tsv.gz", sep ="\t")

# Filter (and align) to hm3 variants manifest
setnames(t2d, c("CHR","BP","REF","ALT", "BETA"), c("chr","pos","a0","a1","beta"))
t2d.hm3 <- as.data.table(snp_match(t2d, hm3, strand_flip=FALSE,  match.min.prop=0.4))
#6,550,140 variants to be matched.
#1,046,689 variants have been matched; 0 were flipped and 508,688 were reversed.
setnames(t2d.hm3,c("chr","pos","a0","a1","beta"), c("CHR","BP","REF","ALT", "BETA"))
t2d.hm3[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]

# Save file with HapMap3 variants
fwrite(t2d.hm3, "../datasets/T2D-hm3.tsv.gz", sep="\t")

rm(t2d, t2d.hm3, t2d.1kg, sd_ss,sd_val, is_bad, invarlist)


## BRCA (Breast Cancer) dataset

brca <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
	      select = c("chr", "position_b37", "a0", "a1", "bcac_onco_icogs_gwas_eaf_controls", "bcac_onco_icogs_gwas_beta", "bcac_onco_icogs_gwas_se","bcac_onco_icogs_gwas_P1df"),
             col.names = c("CHR", "BP", "REF", "ALT", "ALT_FREQ",
                                 "BETA", "SE", "P"),
	     na.strings="NULL")
# NOTE: This dataset includes a meta-analysis from multiple European and East Asian patients and controls.
# It includes ALT_FREQ in the controls, but it also has imputed SNPs in source GWAS.

brca[,n_eff:=4 / (1 / 137045 + 1 / 119078)]

## Apply QC recommended by LDpred2 authors.

## NOTE: The following line wasn't present in the LDpred2 paper code, but I include it here for consistency.
brca <- brca[pmin(ALT_FREQ, 1-ALT_FREQ) > (1 / sqrt(n_eff)),]


sd_ss <- with(brca, 2 / sqrt(n_eff * SE^2))
sd_val <- with(brca, sqrt(2*ALT_FREQ*(1-ALT_FREQ)))

is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05 
table(is_bad)
# is_bad
#    FALSE     TRUE 
# 10234809  1557549 

qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

ggsave("../figures/sd-approx-BRCA.png", width = 10, height = 7)

brca <- brca[!is_bad,]

# Filter by 1000G QC variants
# This dataset didn't have SNPID, so let's create it
brca[,SNPID:=paste(CHR,BP,sep=":")]
invarlist <- brca$SNPID %in% varlist
brca.1kg <- brca[invarlist,]

# Write 1000G qc-filtered file
fwrite(brca.1kg, "../datasets/BRCA-1kg.tsv.gz", sep ="\t")

# Filter (and align) to hm3 variants manifest
setnames(brca, c("CHR","BP","REF","ALT", "BETA"), c("chr","pos","a0","a1","beta"))
brca.hm3 <- as.data.table(snp_match(brca, hm3, strand_flip=FALSE))
#10,234,809 variants to be matched.
#1,053,580 variants have been matched; 0 were flipped and 0 were reversed.
setnames(brca.hm3,c("chr","pos","a0","a1","beta"), c("CHR","BP","REF","ALT", "BETA"))
brca.hm3[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]

# Save file with HapMap3 variants
fwrite(brca.hm3, "../datasets/BRCA-hm3.tsv.gz", sep="\t")

rm(brca, brca.hm3, brca.1kg, sd_ss,sd_val, is_bad, invarlist)


## PRCA (Prostate Cancer) dataset

prca <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/SchumacherFR_29892016_GCST006085/meta_v3_onco_euro_overall_ChrAll_1_release.txt", 
	       select = c("Chr", "position", "Allele1", "Allele2", "Effect", "StdErr",
             "Pvalue", "Freq1"),
  col.names = c("CHR", "BP", "ALT", "REF", "BETA", "SE", "P", "ALT_FREQ"))
prca[,c("REF","ALT"):=list(toupper(REF),toupper(ALT))][,n_eff:= 4 / (1 / 79194 + 1 / 61112)]
# NOTE: I updated the numbers here. After checking the paper, they slightly differed from what Privé et al., used

## Apply QC recommended by LDpred2 authors.

prca <- prca[pmin(ALT_FREQ, 1-ALT_FREQ) > (1 / sqrt(n_eff)),]

sd_ss <- with(prca, 2 / sqrt(n_eff * SE^2))
sd_val <- with(prca, sqrt(2*ALT_FREQ*(1-ALT_FREQ)))

is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05 
table(is_bad)
# is_bad
#    FALSE     TRUE 
# 10554634  3056442 

qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  #coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

ggsave("../figures/sd-approx-PRCA.png", width = 10, height = 7)

prca <- prca[!is_bad,]

# Filter by 1000G QC variants
# This dataset didn't have SNPID, so let's create it
prca[,SNPID:=paste(CHR,BP,sep=":")]
invarlist <- prca$SNPID %in% varlist
prca.1kg <- prca[invarlist,]

# Write 1000G qc-filtered file
fwrite(prca.1kg, "../datasets/PRCA-1kg.tsv.gz", sep ="\t")

# Filter (and align) to hm3 variants manifest
setnames(prca, c("CHR","BP","REF","ALT", "BETA"), c("chr","pos","a0","a1","beta"))
prca.hm3 <- as.data.table(snp_match(prca, hm3, strand_flip=FALSE))
#10,234,809 variants to be matched.
#1,053,580 variants have been matched; 0 were flipped and 0 were reversed.
setnames(prca.hm3,c("chr","pos","a0","a1","beta"), c("CHR","BP","REF","ALT", "BETA"))
prca.hm3[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]

# Save file with HapMap3 variants
fwrite(prca.hm3, "../datasets/PRCA-hm3.tsv.gz", sep="\t")

rm(prca, prca.hm3, prca.1kg, sd_ss,sd_val, is_bad, invarlist)


## MDD (Major depression) dataset

mdd <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/WrayNR_29700475_GCST005839/MDD2018_ex23andMe.gz",
	     fill = TRUE,
                   select = c("CHR", "BP", "A1", "A2", "OR", "SE", "P", "Nca", "Nco"),
                   col.names = c("CHR", "BP", "ALT", "REF", "OR", "SE", "P", "N1", "N0"))
mdd  <- na.omit(mdd)
mdd[,BETA:=log(OR)][,OR:=NULL][,chr:=as.numeric(CHR)][,BP:=as.numeric(BP)][,n_eff:=4 / (1 / N1 + 1 / N0)]
mdd  <- mdd[n_eff > (0.5 * max(n_eff)),]

## NOTE: This file didn't have MAF available, so we'll need to compute them using an external tool.
fwrite(mdd, "../datasets/MDD-NOFREQS.tsv.gz", sep="\t")
system("Compute_freqs.sh -f ../datasets/MDD-NOFREQS.tsv.gz -p CEU")
# Recover the file
mdd <- fread("../datasets/MDD-withfreqs.tsv.gz")

## Apply QC recommended by LDpred2 authors.

sd_ss <- with(mdd, 2 / sqrt(n_eff * SE^2))
sd_val <- with(mdd, sqrt(2*ALT_FREQ*(1-ALT_FREQ)))

is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05 
table(is_bad)
# is_bad
#   FALSE    TRUE
# 6617035   29400

qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")
ggsave("../figures/sd-approx-MDD.png", width = 10, height = 7)

mdd <- mdd[!is_bad,]
mdd[,chr:=NULL]

# Filter by 1000G QC variants
# This dataset didn't have SNPID, so let's create it
mdd[,SNPID:=paste(CHR,BP,sep=":")]
invarlist <- mdd$SNPID %in% varlist
mdd.1kg <- mdd[invarlist,]

# Write 1000G qc-filtered file
fwrite(mdd.1kg, "../datasets/MDD-1kg.tsv.gz", sep ="\t")

# Filter (and align) to hm3 variants manifest
setnames(mdd, c("CHR","BP","REF","ALT", "BETA"), c("chr","pos","a0","a1","beta"))
mdd.hm3 <- as.data.table(snp_match(mdd, hm3, strand_flip=FALSE))
# 6,617,035 variants to be matched.
# 1,043,246 variants have been matched; 0 were flipped and 506,944 were reversed.
setnames(mdd.hm3,c("chr","pos","a0","a1","beta"), c("CHR","BP","REF","ALT", "BETA"))
mdd.hm3[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]

# Save file with HapMap3 variants
fwrite(mdd.hm3, "../datasets/MDD-hm3.tsv.gz", sep="\t")

rm(mdd, mdd.hm3, mdd.1kg, sd_ss,sd_val, is_bad, invarlist)



## CAD (Coronary artery disease) dataset 

cad <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/NikpayM_26343387_GCST003116/cad.add.160614.website.txt",
	select = c("chr", "bp_hg19", "noneffect_allele",
                              "effect_allele", "beta", "se_dgc", "p_dgc"),
        col.names = c("CHR", "BP", "REF", "ALT", "BETA", "SE", "P"))
cad[,n_eff:=4 / (1 / 60801 + 1 / 123504)]

## NOTE: This file didn't have MAF available, so we'll need to compute them using an external tool.
fwrite(cad, "../datasets/CAD-NOFREQS.tsv.gz", sep="\t")
system("Compute_freqs.sh -f ../datasets/CAD-NOFREQS.tsv.gz -p CEU")
# Recover the file
cad <- fread("../datasets/CAD-withfreqs.tsv.gz")

## Apply QC recommended by LDpred2 authors.

sd_ss <- with(cad, 2 / sqrt(n_eff * SE^2))
sd_val <- with(cad, sqrt(2*ALT_FREQ*(1-ALT_FREQ)))

is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05 
table(is_bad)
# is_bad
#   FALSE    TRUE 
# 6278405  670834 

qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")
ggsave("../figures/sd-approx-CAD.png", width = 10, height = 7)

cad <- cad[!is_bad,]

# Filter by 1000G QC variants
# This dataset didn't have SNPID, so let's create it
cad[,SNPID:=paste(CHR,BP,sep=":")]
invarlist <- cad$SNPID %in% varlist
cad.1kg <- cad[invarlist,]

# Write 1000G qc-filtered file
fwrite(cad.1kg, "../datasets/CAD-1kg.tsv.gz", sep ="\t")

# Filter (and align) to hm3 variants manifest
setnames(cad, c("CHR","BP","REF","ALT", "BETA"), c("chr","pos","a0","a1","beta"))
cad.hm3 <- as.data.table(snp_match(cad, hm3, strand_flip=FALSE))
# 6,278,405 variants to be matched.
# 1,044,311 variants have been matched; 0 were flipped and 725,285 were reversed.
setnames(cad.hm3,c("chr","pos","a0","a1","beta"), c("CHR","BP","REF","ALT", "BETA"))
cad.hm3[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]

# Save file with HapMap3 variants
fwrite(cad.hm3, "../datasets/CAD-hm3.tsv.gz", sep="\t")

rm(cad, cad.hm3, cad.1kg, sd_ss,sd_val, is_bad, invarlist)



## BMI (Body mass Index) dataset

bmi <- fread("https://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz",
	col.names=c("SNPID", "ALT", "REF", "ALT_FREQ", "BETA", "SE", "P","N"))

# Let's compute n_eff
bmi[,sdy:=sdY.est(SE^2,ALT_FREQ, N)][,n_eff:=4 * N /sdy^2]

## NOTE: This dataset doesn't have coordinates, we'll need to compute them from panel, using 00a-Missing_coordinates. Then we'll reimport it and replace SNPID
# Rscript ../../../GWAS_tools/00a-Missing_coordinates/Fetch_coordinates.R BMI_Locke_25673413_1-nofiltNOCOORDS.tsv.gz
fwrite(bmi, "../datasets/BMI-NOCOORDS.tsv.gz", sep="\t")

# Retrieve coordinates 
system("Rscript --vanilla Fetch_coordinates.R ../datasets/BMI-NOCOORDS.tsv.gz")

bmi <- fread("../datasets/BMI-withcoords.tsv.gz")
bmi[,SNPID:=paste(CHR,BP, sep = ":")]
bmi <- bmi[,c("SNPID","CHR","BP", "REF","ALT","ALT_FREQ","BETA","SE","P","N", "n_eff", "sdy")]
bmi <- na.omit(bmi)

## Apply QC recommended by LDpred2 authors.
# In this particular case, we'll apply a bit different QC, since it's a continuous trait.
# I computed sdy already, so I'll filter SNPs following advice found at https://github.com/privefl/bigsnpr/issues/168
sd_ss <- with(bmi, sdy / sqrt(n_eff * SE^2))
sd_val <- with(bmi, sqrt(2*ALT_FREQ*(1-ALT_FREQ)))

is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05 
table(is_bad)
# is_bad
#   FALSE    TRUE 
# 1878199   97509 

qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")
ggsave("../figures/sd-approx-BMI.png", width = 10, height = 7)

bmi <- bmi[!is_bad,]

# Filter by 1000G QC variants
invarlist <- bmi$SNPID %in% varlist
bmi.1kg <- bmi[invarlist,]

fwrite(bmi.1kg, "../datasets/BMI-1kg.tsv.gz", sep="\t")

# Filter (and align) to hm3 variants manifest
setnames(bmi, c("CHR","BP","REF","ALT", "BETA"), c("chr","pos","a0","a1","beta"))
bmi.hm3 <- as.data.table(snp_match(bmi, hm3, strand_flip=FALSE))
# 1,878,199 variants to be matched.
# Some duplicates were removed.  
# 763,023 variants have been matched; 0 were flipped and 376,261 were reversed.
setnames(bmi.hm3,c("chr","pos","a0","a1","beta"), c("CHR","BP","REF","ALT", "BETA"))
bmi.hm3[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]

# Save file with HapMap3 variants
fwrite(bmi.hm3, "../datasets/BMI-hm3.tsv.gz", sep="\t")


rm(bmi, bmi.hm3, bmi.1kg, invarlist)



## Height dataset

height <- fread("https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz",
	col.names=c("SNPID", "ALT", "REF", "ALT_FREQ", "BETA", "SE", "P","N"))

### Let's compute n_eff
height[,sdy:=sdY.est(SE^2,ALT_FREQ, N)][,n_eff:=4 * N /sdy^2]

## NOTE: This dataset doesn't have coordinates, we'll need to compute them from panel, using 00a-Missing_coordinates. Then we'll reimport it and replace SNPID
fwrite(height, "../datasets/Height-NOCOORDS.tsv.gz", sep="\t")

# Retrieve coordinates 
system("Rscript --vanilla Fetch_coordinates.R ../datasets/Height-NOCOORDS.tsv.gz")


height <- fread("../datasets/Height-withcoords.tsv.gz")
height[,SNPID:=paste(CHR,BP, sep = ":")]
height <- height[,c("SNPID","CHR","BP", "REF","ALT","ALT_FREQ","BETA","SE","P","N", "n_eff", "sdy")]
height <- na.omit(height)

## Apply QC recommended by LDpred2 authors.
# In this particular case, we'll apply a bit different QC, since it's a continuous trait.
# I computed sdy already, so I'll filter SNPs following advice found at https://github.com/privefl/bigsnpr/issues/168
sd_ss <- with(height, sdy / sqrt(n_eff * SE^2))
sd_val <- with(height, sqrt(2*ALT_FREQ*(1-ALT_FREQ)))

is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05 
table(is_bad)
# is_bad
#   FALSE    TRUE 
# 1443721  567305 

qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")
ggsave("../figures/sd-approx-HEIGHT.png", width = 10, height = 7)

height <- height[!is_bad,]

# Filter by 1000G QC variants
invarlist <- height$SNPID %in% varlist
height.1kg <- height[invarlist,]

fwrite(height.1kg, "../datasets/Height-1kg.tsv.gz", sep="\t")

# Filter (and align) to hm3 variants manifest
setnames(height, c("CHR","BP","REF","ALT", "BETA"), c("chr","pos","a0","a1","beta"))
height.hm3 <- as.data.table(snp_match(height, hm3, strand_flip=FALSE))
# 1,443,721 variants to be matched.
# 581,721 variants have been matched; 0 were flipped and 282,777 were reversed.

setnames(height.hm3,c("chr","pos","a0","a1","beta"), c("CHR","BP","REF","ALT", "BETA"))
height.hm3[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]

# Save file with HapMap3 variants
fwrite(height.hm3, "../datasets/Height-hm3.tsv.gz", sep="\t")


rm(height, height.hm3, height.1kg, invarlist)



# RA (International)

# RA (European) failed for LDpred2-auto, so we'll fall back to the transethnic dataset (RA_Okada_3) we used in previous versions (and that LDpred2 authors used in their own paper).
raint <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/OkadaY_24390342_GCST002318/RA_GWASmeta_TransEthnic_v2.txt.gz", select = 2:7,
                   col.names = c("CHR", "BP", "ALT", "REF", "OR", "P"))
raint[,BETA:=log(OR)][,chi2:= qchisq(P, df = 1, lower.tail = FALSE)][,SE:= ifelse(chi2 > 1e-4, abs(BETA) / sqrt(chi2), NA)][, n_eff:= 4 / (1 / 19234 + 1 / 61565)] # Updated numbers here too


## This file didn't have MAF available, so we'll need to compute them using an external tool.
fwrite(raint, "../datasets/RAint-NOFREQS.tsv.gz", sep="\t")
system("Compute_freqs.sh -f ../datasets/RAint-NOFREQS.tsv.gz -p CEU")

## Recover the file
raint <- fread("../datasets/RAint-withfreqs.tsv.gz")

## Apply QC recommended by LDpred2 authors.

sd_ss <- with(raint, 2 / sqrt(n_eff * SE^2))
sd_val <- with(raint, sqrt(2*ALT_FREQ*(1-ALT_FREQ)))

is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05 
table(is_bad)
# is_bad
#   FALSE    TRUE 
# 4992972 3183373 

qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

ggsave("../figures/sd-approx-RAint.png", width = 10, height = 7)

raint <- raint[!is_bad,]

raint[,SNPID:=paste(CHR,BP,sep=":")]
raint[,c("OR","chi2"):=NULL]

# Filter by 1000G QC variants
invarlist <- raint$SNPID %in% varlist
raint.1kg <- raint[invarlist,]

# Write 1000G qc-filtered file
fwrite(raint.1kg, "../datasets/RAint-1kg.tsv.gz", sep ="\t")

# Filter (and align) to hm3 variants manifest
setnames(raint, c("CHR","BP","REF","ALT", "BETA"), c("chr","pos","a0","a1","beta"))
raint.hm3 <- as.data.table(snp_match(raint, hm3, strand_flip=FALSE))
# 4,992,972 variants to be matched.
# 728,597 variants have been matched; 0 were flipped and 354,315 were reversed.
setnames(raint.hm3,c("chr","pos","a0","a1","beta"), c("CHR","BP","REF","ALT", "BETA"))
raint.hm3[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]

# Save file with HapMap3 variants
fwrite(raint.hm3, "../datasets/RAint-hm3.tsv.gz", sep="\t")

rm(raint, raint.hm3, raint.1kg, sd_ss, sd_val, is_bad, invarlist)



