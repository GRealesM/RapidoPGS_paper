## Generating a Genome-wide reference panel and aligning our summary statistics dataset to it

# This code will use 1000 Genomes Phase III reference panel to align our sumstats, thus (hopefully) improving results at coloc.
# Note that the first step will take some space (~50GB), so please make sure you have enough space.

##################################
### LOADING REQUIRED LIBRARIES ###
##################################

library(data.table)
library(bigsnpr)

system_verbose <- function(..., verbose) {
	system(..., ignore.stdout = !verbose, ignore.stderr = !verbose)
}

# Remove annoying timeout limit in download.file
timeout <- getOption('timeout')
options(timeout=10000)

##################################
###      DOWNLOADING PANEL     ###
##################################

# In this step we'll download 1000 Genomes Phase III in vcf format, and we'll subsequently convert it to Plink format. 

dir.create("../references/RapidoPGS-ref")
plink <- download_plink("../references/RapidoPGS-ref")

# Downloading 1000Genomes phase III
# These are (very) big files, so they'll take a while to download. Time for a cuppa.
# Here we'll consider autosomes only, but we could add X chromosomes, too.
for(chr in c(1:22)){
	download.file(paste0("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"),
              destfile = paste0("../references/RapidoPGS-ref/chr",chr,".vcf.gz"), mode = "wb")
}

# X and Y chromosomes have a bit different link, in case you'd like to use it
#download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz", destfile = "../references/RapidoPGS-ref/chr23.vcf.gz", mode = "wb")
#download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz", destfile = "../references/RapidoPGS-ref/chr24.vcf.gz", mode = "wb")

# Transform them to plink format. 
for(chr in 1:22){
	system_verbose(paste0(plink, " --vcf ../references/RapidoPGS-ref/chr", chr, ".vcf.gz --make-bed --out ../references/RapidoPGS-ref/chr", chr), verbose = TRUE)
}
# We don't need our hard-earned vcf files anymore, so we can delete them
unlink("../references/RapidoPGS-ref/*vcf.gz")


##############################
#### APPLYING QC TO FILES ####
##############################

ped  <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped")
#ped <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
pop <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv")

for(chr in 1:22){

    bed <- snp_plinkQC(plink, prefix.in = paste0("../references/RapidoPGS-ref/chr",chr),
                       prefix.out = paste0("../references/RapidoPGS-ref/chr",chr,"_QC"),
                       geno = 0, maf = 0.01, hwe = 1e-10)
    rds <- snp_readBed(bed)
    snp <- snp_attach(paste0("../references/RapidoPGS-ref/chr",chr,"_QC.rds"))  
    G <- snp$genotypes
    fam <- snp$fam
    fam <- dplyr::left_join(fam, ped, by = c("sample.ID" = "Individual ID"))
    fam <- dplyr::left_join(fam, pop, by = c("Population" = "Population Code"))
    
    snp$fam$family.ID <- paste(fam$`Super Population`, fam$Population, sep = "_")
    snp$fam$paternal.ID <- snp$fam$maternal.ID <- 0L

    maf <- snp_MAF(G)    
    bed <- snp_writeBed(snp, paste0("../references/RapidoPGS-ref/chr",chr,"_tmp.bed"),
                          ind.col = which(maf > 0.01))
    rds <- snp_readBed(bed)
    snp <- snp_attach(rds)

    unlink(paste0("../references/RapidoPGS-ref/chr", chr,"\\.*")) # Remove original plink files
    snp_writeBed(snp, paste0("../references/RapidoPGS-ref/chr",chr,".bed"))
    unlink("../references/RapidoPGS-ref/*_QC.*")
    unlink("../references/RapidoPGS-ref/*_tmp.*")
    message("Done!")
}
  unlink("../references/RapidoPGS-ref/plink")
  message("Done! Now you can find your reference panel ready at ../references/RapidoPGS-ref/.")

################################################################################
## EXTRA STEP: Generate a list of variants in the panel to filter datasets by ##
################################################################################

system("for i in {1..22}; do cat ../references/RapidoPGS-ref/chr$i.bim >> ../references/RapidoPGS-ref/variant_list.txt; done")
variant_list <- fread("../references/RapidoPGS-ref/variant_list.txt")
variant_list[,pid:=paste(V1,V4, sep=":")] 
variant_list <- variant_list[,pid]
fwrite(list(variant_list), "../datasets/RapidoPGS_panel_variants.txt", col.names=F)



