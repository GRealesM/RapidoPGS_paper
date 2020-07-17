# This script will serve to process new datasets (BRCA, PRCA, MDD, CAD, RA and T1D-HLA, and their corresponding UKBB_Neale counterparts for analysis) 
# and as a template for future clean code.
# We had previously downloaded and processed some of the datasets (Asthma, T1D, and T2D, from another project, so the code to download + preprocess those 
# may be a bit different albeit equivalent to that showed here.
# Bear in mind that there is an important step prior to PGS generation that must be performed, as many of these datasets lack allele frequencies.
# We computed them using 1000Genomes Phase III CEU population, using other scripts. Datasets with no allele frequencies are saved as "NOFREQS", and the frequency computing step is due.
# Regarding the missing UKBB Neale datasets, they were processed exactly the same way. I will include the rest of the code for completeness.


library(data.table)
library(bigsnpr)

## We'll apply a filter to our datasets, using the SNPs in sAUC reference panel
## Allegedly 1000 genomes panel I, although not specified in the paper.
## This panel contains 16M SNPs, and 379 individuals of European origin
sauc.map  <- fread("../references/reference_sAUC.txt", col.names=c("SNPID","chr","pos", "a1","a0"))

## BRCA (Breast Cancer) dataset

brca <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
	      select = c("chr", "position_b37", "a0", "a1", "bcac_onco_icogs_gwas_eaf_controls", "bcac_onco_icogs_gwas_beta", "bcac_onco_icogs_gwas_se","bcac_onco_icogs_gwas_P1df"),
             col.names = c("chr", "pos", "a0", "a1", "freq",
                                 "beta", "beta_se", "p"),
	     na.strings="NULL")

brca.aligned <- as.data.table(snp_match(brca,sauc.map))
#11,792,542 variants to be matched.
#1,633,946 ambiguous SNPs have been removed.
#7,978,580 variants have been matched; 0 were flipped and 1,328,824 were reversed.
brca.aligned[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]
names(brca.aligned) <- c("CHR19","BP19", "REF","ALT","ALT_FREQ","BETA","SE","P","SNPID")
fwrite(brca.aligned, "../datasets/BRCA_Michailidou_29059683_1-sAUCfilt.tsv.gz", sep="\t")


## PRCA (Prostate Cancer) dataset

prca <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/SchumacherFR_29892016_GCST006085/meta_v3_onco_euro_overall_ChrAll_1_release.txt", 
	       select = c("Chr", "position", "Allele1", "Allele2", "Effect", "StdErr",
             "Pvalue", "Freq1"),
  col.names = c("chr", "pos", "a1", "a0", "beta", "beta_se", "p", "freq"))
prca[,c("a0","a1"):=list(toupper(a0),toupper(a1))]

prca.aligned <- as.data.table(snp_match(prca,sauc.map))
#20,370,946 variants to be matched.
#2,866,342 ambiguous SNPs have been removed.
#10,115,631 variants have been matched; 497 were flipped and 4,276,822 were reversed.

prca.aligned[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]
names(prca.aligned) <- c("CHR19","BP19", "REF","ALT","BETA","SE","P","ALT_FREQ","SNPID")
fwrite(prca.aligned, "../datasets/PRCA_Schumacher_29892016_1-sAUCfilt.tsv.gz", sep="\t")



#### MDD (Major depression) dataset #####

mdd <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/WrayNR_29700475_GCST005839/MDD2018_ex23andMe.gz",
	     fill = TRUE,
                   select = c("CHR", "BP", "A1", "A2", "OR", "SE", "P", "Nca", "Nco"),
                   col.names = c("chr", "pos", "a1", "a0", "or", "beta_se", "p", "N1", "N0"))
mdd[,beta:=log(or)]
mdd[,or:=NULL]
mdd[,chr:=as.numeric(chr)]

mdd.aligned <- as.data.table(snp_match(mdd,sauc.map))
#13,554,550 variants to be matched.
#1,868,797 ambiguous SNPs have been removed. 
#9,758,197 variants have been matched; 552 were flipped and 4,132,132 were reversed
mdd.aligned[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]
names(mdd.aligned) <- c("CHR19","BP19", "REF","ALT","SE","P","N1","N0","BETA","SNPID")
# NOTE: This file didn't have MAF available, so we'll need to compute them using an external tool.
fwrite(mdd.aligned, "../datasets/MDD_Wray_29700475_1-sAUCfiltNOFREQS.tsv.gz", sep="\t")


### CAD (Coronary artery disease) ###

cad <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/NikpayM_26343387_GCST003116/cad.add.160614.website.txt",
	select = c("chr", "bp_hg19", "noneffect_allele",
                              "effect_allele", "beta", "se_dgc", "p_dgc"),
        col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se", "p"))
cad.aligned  <- as.data.table(snp_match(cad,sauc.map))
#9,455,778 variants to be matched.
#1,332,725 ambiguous SNPs have been removed.
#7,285,220 variants have been matched; 0 were flipped and 6,505,655 were reversed
cad.aligned[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]
names(cad.aligned) <- c("CHR19","BP19", "REF","ALT","BETA","SE","P","SNPID")
# NOTE: This file didn't have MAF available, so we'll need to compute them using an external tool.
fwrite(cad.aligned, "../datasets/CAD_Nikpay_26343387_1-sAUCfiltNOFREQS.tsv.gz", sep="\t")


## T1D No HLA

# We already processed this dataset in a previous step (not depicted here - yet), so we can simply load it
t1d <- fread("../datasets/T1D_Cooper_1-sAUCfilt.tsv.gz")

hlaloci <- which(t1d$CHR19 == 6 & t1d$BP19 > 20000000 & t1d$BP19 < 40000000)
t1d.nohla <- t1d[-hlaloci,]
fwrite(t1d.nohla, "../datasets/T1DnoHLA_Cooper_1-sAUCfilt.tsv.gz",sep="\t")



## RA - The right dataset
# I was using RA_Okada_1, while the right dataset is actually RA_Okada_3, let's correct that. Since I downloaded RA_Okada_3 from harmonised files at GWAS catalog, and it didn't contain original, hg19 coordinates, I will download the original, unharmonised one

ra <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/OkadaY_24390342_GCST002318/RA_GWASmeta_TransEthnic_v2.txt.gz", 
	select=1:7,
	col.names = c("SNPID","chr", "pos", "a1", "a0", "or", "p"))
ra[,beta:=log(or)][,chi2:= qchisq(p, df = 1, lower.tail = FALSE)][,beta_se:= ifelse(chi2 > 1e-4, abs(beta) / sqrt(chi2), NA)]

ra.aligned  <- as.data.table(snp_match(ra,sauc.map))
# 9,739,303 variants to be matched.
# 1,498,996 ambiguous SNPs have been removed.
# 7,384,567 variants have been matched; 0 were flipped and 3,243,960 were reversed.

ra.aligned[,c("_NUM_ID_.ss","_NUM_ID_"):=NULL]
ra.aligned[,c("SNPID.ss","or","chi2"):=NULL]
names(ra.aligned) <- c("CHR19","BP19", "REF","ALT","P","BETA","SE","SNPID")
# NOTE: This file didn't have MAF available, so we'll need to compute them using an external tool.
fwrite(ra.aligned, "../datasets/RA_Okada_24390342_3-sAUCfiltNOFREQS.tsv.gz", sep="\t")





###### UKBB validation datasets ########

## BRCA UKBB ##

# NOTE: This is only one of several possible phenotypes related with cancer. This corresponds to C50: diagnosis of malignant neoplasm of breast (females)
download.file("https://www.dropbox.com/s/oleri4jlu77e6oh/C50.gwas.imputed_v3.female.tsv.bgz?dl=1", destfile="../datasets/BRCA_Neale_UKBB_1.tsv.gz") 
brca.ukbb  <- fread("../datasets/BRCA_Neale_UKBB_1.tsv.gz")

# According to the manifest in https://www.dropbox.com/s/0a87i4y049vje12/phenotypes.female.v2.tsv.bgz?dl=0 
# N0 = 185928
# N1 = 8246 


brca.ukbb[,c("CHR19","BP19","REF","ALT"):=tstrsplit(variant, ":")][, pid:=paste(CHR19,BP19, sep=":")]
brca.ukbb <- brca.ukbb[CHR19 != "X",]
brca.ukbb <- brca.ukbb[,c("pid", "ALT", "minor_AF", "beta", "pval")]
setnames(brca.ukbb, old=c("minor_AF","beta","pval"), new=c("MAF", "BETA", "P"))
brca.ukbb <- na.omit(brca.ukbb)
fwrite(brca.ukbb, "../datasets/BRCA_Neale_UKBB_1-forvalidation.tsv.gz", sep="\t")

## PRCA UKBB ##

# NOTE: This is only one of several possible phenotypes related with cancer. This corresponds to C61: diagnosis of malignant neoplasm of prostate (males)
download.file("https://www.dropbox.com/s/x5qh2uwbslcf84k/C61.gwas.imputed_v3.male.tsv.bgz?dl=1", destfile="../datasets/PRCA_Neale_UKBB_1.tsv.gz") 
prca.ukbb  <- fread("../datasets/PRCA_Neale_UKBB_1.tsv.gz")

# According to the manifest in https://www.dropbox.com/s/k8h7j85awav0lrt/phenotypes.male.v2.tsv.bgz?dl=0
# N0 = 162678	
# N1 = 4342

prca.ukbb[,c("CHR19","BP19","REF","ALT"):=tstrsplit(variant, ":")][, pid:=paste(CHR19,BP19, sep=":")]
prca.ukbb <- prca.ukbb[CHR19 != "X",]
prca.ukbb <- prca.ukbb[,c("pid", "ALT", "minor_AF", "beta", "pval")]
setnames(prca.ukbb, old=c("minor_AF","beta","pval"), new=c("MAF", "BETA", "P"))
fwrite(prca.ukbb, "../datasets/PRCA_Neale_UKBB_1-forvalidation.tsv.gz", sep="\t")


## MDD UKBB ##

# NOTE: This is only one of several possible phenotypes related with cancer. This corresponds to 20002_1286 (both sexes).
download.file("https://www.dropbox.com/s/iti2uh5tvhkphg5/20002_1286.gwas.imputed_v3.both_sexes.tsv.bgz?dl=1", destfile="../datasets/MDD_Neale_UKBB_1.tsv.gz") 
mdd.ukbb  <- fread("../datasets/MDD_Neale_UKBB_1.tsv.gz")

# According to the manifest in https://www.dropbox.com/s/j1h9e0lblmukiko/phenotypes.both_sexes.v2.tsv.bgz?dl=0
# N0 = 340493	
# N1 = 20648

mdd.ukbb[,c("CHR19","BP19","REF","ALT"):=tstrsplit(variant, ":")][, pid:=paste(CHR19,BP19, sep=":")]
mdd.ukbb <- mdd.ukbb[CHR19 != "X",]
mdd.ukbb <- mdd.ukbb[,c("pid", "ALT", "minor_AF", "beta", "pval")]
setnames(mdd.ukbb, old=c("minor_AF","beta","pval"), new=c("MAF", "BETA", "P"))
fwrite(mdd.ukbb, "../datasets/MDD_Neale_UKBB_1-forvalidation.tsv.gz", sep="\t")

## CAD UKBB ##

# NOTE: This is only one of several possible phenotypes related with coronary heart disease. This corresponds to I9CHD (both sexes). 
download.file("https://www.dropbox.com/s/yxyk4t6df6i4bgg/I9_CHD.gwas.imputed_v3.both_sexes.tsv.bgz?dl=1", destfile="../datasets/CAD_Neale_UKBB_1.tsv.gz") 
cad.ukbb  <- fread("../datasets/CAD_Neale_UKBB_1.tsv.gz")

# According to the manifest in https://www.dropbox.com/s/j1h9e0lblmukiko/phenotypes.both_sexes.v2.tsv.bgz?dl=0
# N0 = 351037	
# N1 = 10157
cad.ukbb[,c("CHR19","BP19","REF","ALT"):=tstrsplit(variant, ":")][, pid:=paste(CHR19,BP19, sep=":")]
cad.ukbb <- cad.ukbb[CHR19 != "X",]
cad.ukbb <- cad.ukbb[,c("pid", "ALT", "minor_AF", "beta", "pval")]
setnames(cad.ukbb, old=c("minor_AF","beta","pval"), new=c("MAF", "BETA", "P"))
fwrite(cad.ukbb, "../datasets/CAD_Neale_UKBB_1-forvalidation.tsv.gz", sep="\t")






