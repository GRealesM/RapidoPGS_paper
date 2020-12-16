# Getting public data from 1000genomes project.
# This script will try to follow the instructions to download, process, and understand
# public data in bigsnpr vignette (https://github.com/privefl/bigsnpr/blob/master/data-raw/public-data.R)

# NOTE: In this panel construction, we won't remove related individuals.


##################################
### LOADING REQUIRED LIBRARIES ###
##################################

library(bigsnpr)
library(data.table)

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

# We download the panel
dir.create("../references/LDpred2-ref/")
plink <- download_plink("../references/LDpred2-ref/")
# And the SNPlist
# NOTE: I believe there's a typo, and pos column is actually hg19. I don't know which build the other columns correspond to, but they don't seem to match
info <- readRDS(url("https://github.com/privefl/bigsnpr/raw/master/data-raw/hm3_variants.rds"))


# Downloading 1000Genomes phase III
# These are (very) big files, so they'll take a while to download. Time for a cuppa.
# Here we'll consider autosomes only, but we could add X chromosomes, too.
for(chr in c(1:22)){
	download.file(paste0("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"),
              destfile = paste0("../references/LDpred2-ref/chr",chr,".vcf.gz"), mode = "wb")
}

# X and Y chromosomes have a bit different link
#download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz", destfile = "../references/LDpred2-ref/chr23.vcf.gz", mode = "wb")
#download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz", destfile = "../references/LDpred2-ref/chr24.vcf.gz", mode = "wb")

# Transform them to plink format. 
for(chr in 1:22){
	system_verbose(paste0(plink, " --vcf ../references/LDpred2-ref/chr", chr, ".vcf.gz --make-bed --out ../references/LDpred2-ref/chr", chr), verbose = TRUE)
}
# We don't need our hard-earned vcf files anymore, so we can delete them
unlink("../references/LDpred2-ref/*vcf.gz")

ped  <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped")
#ped <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
pop <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv")

for(chr in 1:22){

    bed <- snp_plinkQC(plink, prefix.in = paste0("../references/LDpred2-ref/chr",chr),
                       prefix.out = paste0("../references/LDpred2-ref/chr",chr,"_QC"),
                       geno = 0, maf = 0.01, hwe = 1e-10)
    rds <- snp_readBed(bed)
    snp <- snp_attach(paste0("../references/LDpred2-ref/chr",chr,"_QC.rds"))  
    G <- snp$genotypes
    fam <- snp$fam
    fam <- dplyr::left_join(fam, ped, by = c("sample.ID" = "Individual ID"))
    fam <- dplyr::left_join(fam, pop, by = c("Population" = "Population Code"))
    
    snp$fam$family.ID <- paste(fam$`Super Population`, fam$Population, sep = "_")
    snp$fam$paternal.ID <- snp$fam$maternal.ID <- 0L
    
    # Only in HapMap3list, special filter for LDpred2
    info.chr  <- info[info$chr == chr,] 

    maf <- snp_MAF(G)    
    bed <- snp_writeBed(snp, paste0("../references/LDpred2-ref/chr",chr,"_tmp.bed"),
                          ind.col = which(maf > 0.01 & snp$map$physical.pos %in% info.chr$pos)) # Filter LDpred2 panel by HapMap3 variants
    rds <- snp_readBed(bed)
    snp <- snp_attach(rds)
    unlink(paste0("../references/LDpred2-ref/chr", chr,"\\.*"))
    snp_writeBed(snp, paste0("../references/LDpred2-ref/chr",chr,".bed"))
    unlink("../references/LDpred2-ref/*_QC.*")
    unlink("../references/LDpred2-ref/*_tmp.*")
    message("Done!")
}
  unlink("../references/LDpred2-ref/plink")
  message("Done! Now you can find your reference panel ready at ../references/LDpred2-ref/.")

