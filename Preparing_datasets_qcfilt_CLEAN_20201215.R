# Preprocessing datasets
# This script will serve to process datasets used in the paper (Asthma, RA, T1D, T2D, BRCA, PRCA, MDD, CAD, BMI and Height).
# Also, we'll apply QC recommended by Florian Privé in the LDpred2 preprint
# 2020/12/15


library(data.table)
library(bigsnpr)

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
# We'll prefilter all datasets by the variants contained in the RápidoPGS-mult panel that we created, so RápidoPGS-single will use the same variants as RápidoPGS-mult
# This won't affect LDpred2, since that panel experienced the same QC filt (ie. maf > 0.01 and hwe p-value > 1e-10) plus filter by HapMap3 variants, so its more strigent than RápidoPGS panel.
# Datasets will be filtered by LDpred2 at runtime.
varlist <- fread("../datasets/RapidoPGS_panel_variants.txt", header=FALSE)
varlist <- varlist[,V1]

## Asthma dataset

 download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/DemenaisF_29273806_GCST006862/TAGC_meta-analyses_results_for_asthma_risk.zip", destfile = "sumstats_asthma.zip")
 unzip("sumstats_asthma.zip")
# We already downloaded this dataset
asthma <- fread("TAGC meta-analyses results for asthma risk/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv",
  select = c(1, 3:5, 18:20),
  col.names = c("CHR", "BP", "REF", "ALT", "BETA", "SE", "P")
)
asthma$n_eff <- 4 / (1 / 19954 + 1 / 107715)


## Apply QC by removing variants with sd_ss lower than 0.1. Since we don't have a validation dataset, we can't compute sd_val,
# so we won't compute it

sd_ss <- with(asthma, 1 / sqrt(n_eff / 4 * SE^2))

is_bad <- sd_ss < 0.1

asthma <- asthma[!is_bad,]

# Filter by QC variants
asthma[,SNPID:=paste(CHR,BP,sep=":")]
invarlist <- asthma$SNPID %in% varlist
asthma <- asthma[invarlist,]

fwrite(asthma, "../datasets/AST_Demenais_29273806_1-qcfiltNOFREQS.tsv.gz", sep="\t")

####### This step is no longer required, so dataset can be saved as qcfilt.tsv.gz straight away ######
# Let's compute the frequencies using an external tool and the RápidoPGS panel
system("Compute_freqs.sh -f ../datasets/AST_Demenais_29273806_1-qcfiltNOFREQS.tsv.gz -p CEU")
file.rename("../datasets/AST_Demenais_29273806_1-withfreqs.tsv.gz","../datasets/AST_Demenais_29273806_1-qcfilt.tsv.gz")
######################################################################################################

rm(asthma, sd_ss,is_bad, invarlist)


## RA (Rheumatoid Arthritis) dataset

ra <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/OkadaY_24390342_GCST002318/RA_GWASmeta_TransEthnic_v2.txt.gz", 
	select=1:7,
	col.names = c("SNPID","CHR", "BP", "ALT", "REF", "OR", "P"))
ra[,BETA:=log(OR)][,chi2:= qchisq(P, df = 1, lower.tail = FALSE)][,SE:= ifelse(chi2 > 1e-4, abs(BETA) / sqrt(chi2), NA)]

ra$n_eff <- 4 / (1 / 19234 + 1 / 61565)

### Apply QC by removing variants with sd_ss lower than 0.1. Since we don't have a validation dataset, we can't compute sd_val,
## so we won't compute it

sd_ss <- with(ra, 1 / sqrt(n_eff / 4 * SE^2))

is_bad <- sd_ss < 0.1

ra <- ra[!is_bad,]


# Filter by QC variants
ra[,SNPID:=paste(CHR,BP,sep=":")]
invarlist <- ra$SNPID %in% varlist
ra <- ra[invarlist,]

ra[,c("OR","chi2"):=NULL]

fwrite(ra, "../datasets/RA_Okada_24390342_3-qcfiltNOFREQS.tsv.gz", sep="\t")

####### This step is no longer required, so dataset can be saved as qcfilt.tsv.gz straight away ######
# Let's compute the frequencies using an external tool and the RápidoPGS panel
system("Compute_freqs.sh -f ../datasets/RA_Okada_24390342_3-qcfiltNOFREQS.tsv.gz -p CEU")
file.rename("../datasets/RA_Okada_24390342_3-withfreqs.tsv.gz","../datasets/RA_Okada_24390342_3-qcfilt.tsv.gz")
######################################################################################################

rm(ra, sd_ss,is_bad, invarlist)


## T1D (Type 1 diabetes) dataset

## We had access to the original Cooper dataset, which comes in one file.
t1d <- fread("../../../02-Processed/T1D_Cooper_doi101101120022_1-hg38.tsv.gz",
		select=c(2:6,16,17,19,28),
		col.names=c("CHR", "BP","REF", "ALT","ALT_FREQ","BETA","SE","P","qc.check"))

## We apply QC check as in LDpred2 paper, including removing alleles with qc.check different than PASS, and those with ALT_FREQ < 1 / sqrt(n_eff)
t1d <- t1d[qc.check == "PASS",][,n_eff:=4 / (1 / 5913 + 1 / 8828)][ALT_FREQ > 1/sqrt(n_eff),][,qc.check:=NULL]
t1d <- na.omit(t1d)

### Apply QC by removing variants with sd_ss lower than 0.1. Since we don't have a validation dataset, we can't compute sd_val,
## so we won't compute it

sd_ss <- with(t1d, 1 / sqrt(n_eff / 4 * SE^2))

is_bad <- sd_ss < 0.1

t1d <- t1d[!is_bad,]

# Filter by QC variants
t1d[,SNPID:=paste(CHR,BP,sep=":")]
invarlist <- t1d$SNPID %in% varlist
t1d <- t1d[invarlist,]

fwrite(t1d,"../datasets/T1D_Cooper_doi101101120022_1-qcfilt.tsv.gz", sep="\t")
rm(t1d, sd_ss,is_bad, invarlist)


## T2D (Type 2 diabetes) dataset

t2d <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/ScottRA_28566273_GCST004773/METAANALYSIS_DIAGRAM_SE1.txt")
t2d[,c("CHR","BP"):=tstrsplit(`Chr:Position`, ":", fixed = TRUE)]
setnames(t2d, c("Chr:Position", "Allele1", "Allele2", "Effect", "StdErr", "P-value", "TotalSampleSize"), c("SNPID","ALT","REF","BETA","SE","P","N"))

t2d[,n_eff:=4 / (1 / 26676 + 1 / 132532)][,CHR:=as.numeric(CHR)][,BP:=as.numeric(BP)][,N:=NULL]

sd_ss <- with(t2d, 1 / sqrt(n_eff / 4 * SE^2))

is_bad <- sd_ss < 0.1

t2d <- t2d[!is_bad,]
# Reordering
t2d <- t2d[,c("SNPID","CHR","BP", "REF","ALT","BETA","SE","P","n_eff")]

# Filter by QC variants
invarlist <- t2d$SNPID %in% varlist
t2d <- t2d[invarlist,]

fwrite(t2d, "../datasets/T2D_Scott_28566273_1-qcfiltNOFREQS.tsv.gz", sep="\t")

####### This step is no longer required, so dataset can be saved as qcfilt.tsv.gz straight away ######
# Let's compute the frequencies using an external tool and the RápidoPGS panel
system("Compute_freqs.sh -f ../datasets/T2D_Scott_28566273_1-qcfiltNOFREQS.tsv.gz -p CEU")
file.rename("../datasets/T2D_Scott_28566273_1-withfreqs.tsv.gz","../datasets/T2D_Scott_28566273_1-qcfilt.tsv.gz")
######################################################################################################

rm(t2d, sd_ss,is_bad, invarlist)


## BRCA (Breast Cancer) dataset

brca <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
	      select = c("chr", "position_b37", "a0", "a1", "bcac_onco_icogs_gwas_eaf_controls", "bcac_onco_icogs_gwas_beta", "bcac_onco_icogs_gwas_se","bcac_onco_icogs_gwas_P1df"),
             col.names = c("CHR", "BP", "REF", "ALT", "ALT_FREQ",
                                 "BETA", "SE", "P"),
	     na.strings="NULL")

brca[,n_eff:=4 / (1 / 137045 + 1 / 119078)]

sd_ss <- with(brca, 1 / sqrt(n_eff / 4 * SE^2))

is_bad <- sd_ss < 0.1

brca <- brca[!is_bad,]

# Filter by QC variants
brca[,SNPID:=paste(CHR,BP, sep = ":")]
invarlist <- brca$SNPID %in% varlist
brca <- brca[invarlist,]


fwrite(brca, "../datasets/BRCA_Michailidou_29059683_1-qcfilt.tsv.gz", sep="\t")
rm(brca,  sd_ss,is_bad, invarlist) 


## PRCA (Prostate Cancer) dataset

prca <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/SchumacherFR_29892016_GCST006085/meta_v3_onco_euro_overall_ChrAll_1_release.txt", 
	       select = c("Chr", "position", "Allele1", "Allele2", "Effect", "StdErr",
             "Pvalue", "Freq1"),
  col.names = c("CHR", "BP", "ALT", "REF", "BETA", "SE", "P", "ALT_FREQ"))
prca[,c("REF","ALT"):=list(toupper(REF),toupper(ALT))][,n_eff:= 4 / (1 / 79148 + 1 / 61106)]
## This extra step is applied by Privé to some datasets with maf (eg. T1D),
## but not to others (eg. BRCA). To get our runs as accurate as possible,
## we'll apply it here as was in the LDpred2 paper
prca <- prca[pmin(ALT_FREQ, 1-ALT_FREQ) > (1 / sqrt(n_eff)),]

sd_ss <- with(prca, 1 / sqrt(n_eff / 4 * SE^2))

is_bad <- sd_ss < 0.1

prca <- prca[!is_bad,]

# Filter by QC variants
prca[,SNPID:=paste(CHR,BP, sep = ":")]
invarlist <- prca$SNPID %in% varlist
prca <- prca[invarlist,]

fwrite(prca, "../datasets/PRCA_Schumacher_29892016_1-qcfilt.tsv.gz", sep="\t")
rm(prca, sd_ss,is_bad, invarlist) 


## MDD (Major depression) dataset

mdd <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/WrayNR_29700475_GCST005839/MDD2018_ex23andMe.gz",
	     fill = TRUE,
                   select = c("CHR", "BP", "A1", "A2", "OR", "SE", "P", "Nca", "Nco"),
                   col.names = c("CHR", "BP", "ALT", "REF", "OR", "SE", "P", "N1", "N0"))
mdd  <- na.omit(mdd)
mdd[,BETA:=log(OR)][,OR:=NULL][,chr:=as.numeric(CHR)][,BP:=as.numeric(BP)][,n_eff:=4 / (1 / N1 + 1 / N0)]
mdd  <- mdd[n_eff > (0.5 * max(n_eff)),]

sd_ss <- with(mdd, 1 / sqrt(n_eff / 4 * SE^2))

is_bad <- sd_ss < 0.1

mdd <- mdd[!is_bad,]
# Filter by QC variants
mdd[,SNPID:=paste(CHR,BP, sep = ":")]
invarlist <- mdd$SNPID %in% varlist
mdd <- mdd[invarlist,]

fwrite(mdd, "../datasets/MDD_Wray_29700475_1-qcfiltNOFREQS.tsv.gz", sep="\t")

####### This step is no longer required, so dataset can be saved as qcfilt.tsv.gz straight away ######
# Let's compute the frequencies using an external tool and the RápidoPGS panel
system("Compute_freqs.sh -f ../datasets/MDD_Wray_29700475_1-qcfiltNOFREQS.tsv.gz -p CEU")
file.rename("../datasets/MDD_Wray_29700475_1-withfreqs.tsv.gz","../datasets/MDD_Wray_29700475_1-qcfilt.tsv.gz")
######################################################################################################
rm(mdd, sd_ss,is_bad, invarlist)

## CAD (Coronary artery disease) dataset 

cad <- fread("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/NikpayM_26343387_GCST003116/cad.add.160614.website.txt",
	select = c("chr", "bp_hg19", "noneffect_allele",
                              "effect_allele", "beta", "se_dgc", "p_dgc"),
        col.names = c("CHR", "BP", "REF", "ALT", "BETA", "SE", "P"))
cad[,n_eff:=4 / (1 / 60801 + 1 / 123504)]

sd_ss <- with(cad, 1 / sqrt(n_eff / 4 * SE^2))

is_bad <- sd_ss < 0.1

cad <- cad[!is_bad,]

# Filter by QC variants
cad[,SNPID:=paste(CHR,BP, sep = ":")]
invarlist <- cad$SNPID %in% varlist
cad <- cad[invarlist,]

fwrite(cad, "../datasets/CAD_Nikpay_26343387_1-qcfiltNOFREQS.tsv.gz", sep="\t")

####### This step is no longer required, so dataset can be saved as qcfilt.tsv.gz straight away ######
# Let's compute the frequencies using an external tool and the RápidoPGS panel
system("Compute_freqs.sh -f ../datasets/CAD_Nikpay_26343387_1-qcfiltNOFREQS.tsv.gz -p CEU")
file.rename("../datasets/CAD_Nikpay_26343387_1-withfreqs.tsv.gz","../datasets/CAD_Nikpay_26343387_1-qcfilt.tsv.gz")
######################################################################################################
rm(cad, sd_ss,is_bad, invarlist)



## BMI (Body mass Index) dataset

bmi <- fread("https://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz",
	col.names=c("SNPID", "ALT", "REF", "ALT_FREQ", "BETA", "SE", "P","N"))

# Let's compute n_eff
bmi[,sdy:=sdY.est(SE^2,ALT_FREQ, N)][,n_eff:=4 * N /sdy^2]

## NOTE: This dataset doesn't have coordinates, we'll need to compute them from panel, using 00a-Missing_coordinates. Then we'll reimport it and replace SNPID
# Rscript ../../../GWAS_tools/00a-Missing_coordinates/Fetch_coordinates.R BMI_Locke_25673413_1-nofiltNOCOORDS.tsv.gz
fwrite(bmi, "../datasets/BMI_Locke_25673413_1-nofiltNOCOORDS.tsv.gz", sep="\t")

# Retrieve coordinates 
system("Rscript --vanilla Fetch_coordinates.R ../datasets/BMI_Locke_25673413_1-nofiltNOCOORDS.tsv.gz")

bmi <- fread("../datasets/BMI_Locke_25673413_1-withcoords.tsv.gz")

# Filter by QC variants
bmi[,SNPID:=paste(CHR,BP, sep = ":")]
invarlist <- bmi$SNPID %in% varlist
bmi <- bmi[invarlist,]
bmi[,sdy:=NULL]
bmi <- bmi[,c("SNPID","CHR","BP", "REF","ALT","ALT_FREQ","BETA","SE","P","N", "n_eff")]
fwrite(bmi, "../datasets/BMI_Locke_25673413_1-qcfilt.tsv.gz", sep="\t")



## Height dataset

height <- fread("https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz",
	col.names=c("SNPID", "ALT", "REF", "ALT_FREQ", "BETA", "SE", "P","N"))

### Let's compute n_eff
height[,sdy:=sdY.est(SE^2,ALT_FREQ, N)][,n_eff:=4 * N /sdy^2]

fwrite(height, "../datasets/HEIGHT_Wood_25282103_1-nofiltNOCOORDS.tsv.gz", sep="\t")

# Retrieve coordinates 
system("Rscript --vanilla Fetch_coordinates.R ../datasets/HEIGHT_Wood_25282103_1-nofiltNOCOORDS.tsv.gz")


height <- fread("../datasets/HEIGHT_Wood_25282103_1-withcoords.tsv.gz")

# Filter by QC variants
height[,SNPID:=paste(CHR,BP,sep=":")]
invarlist <- height$SNPID %in% varlist
height <- height[invarlist,]
height[,sdy:=NULL]
height <- height[,c("SNPID","CHR","BP", "REF","ALT","ALT_FREQ","BETA","SE","P","N", "n_eff")]
fwrite(height, "../datasets/HEIGHT_Wood_25282103_1-qcfilt.tsv.gz", sep="\t")

