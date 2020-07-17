# This script is intended to filter BMI and Height datasets by the hapmap3 qc'd SNPs by Martin.
# These datasets come without hg coords, so in normal circumstances I'd pull them from 1000genomes,
# but since we'll need to filter them by those qc SNPs anyway, we'll take the coords (from them).

library(data.table)


bmi <- fread("BMI_Locke_25673413_1.tsv.gz")
height <- fread("HEIGHT_Wood_25282103_1.tsv.gz")

bmiref <- fread("reference_hapmap3_BMI_qcSNPs.txt")                           
heightref <- fread("reference_hapmap3_Height_qcSNPs.txt")

names(bmi) <- c("SNPID","ALT","REF","ALT_FREQ","BETA","SE","P", "N")
names(height) <- c("SNPID","ALT","REF","ALT_FREQ","BETA","SE","P", "N")


bmicoords <- merge(bmi, bmiref[, .(SNPID,CHR19,BP19)]) 
heightcoords <- merge(height, heightref[, .(SNPID,CHR19,BP19)]) 


fwrite(bmicoords, "BMI_Locke_25673413_1-qcfilt.tsv.gz", sep="\t")
fwrite(heightcoords, "HEIGHT_Wood_25282103_1-qcfilt.tsv.gz", sep="\t")
