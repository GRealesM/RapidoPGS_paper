# Preprocessing sumstats for SBayesR.
# 2021/03/17

# This script will prepare sumstats datasets for SBayesR

library(data.table)
setDTthreads(8)

traitlist <- c("Asthma","RA","T1D","T2D","BRCA","PRCA","CAD","MDD","BMI","Height", "RAint")
# Define the variable for GWAS size, corresponding to traits:
# Asthma, RA, T1D, T2D, BRCA, PRCA, CAD, MDD, BMI, Height, RAint
# Note: MDD, BMI, and Height have per-SNP N included in the dataset, so we'll use those (see below).
gwasn <- c(127669, 58284, 14741, 159208, 256123, 140306, 184305, NA, NA, NA ,80799)

panel <- readRDS("../references/UKBB_mat/map.rds")
panel$SNPID <- paste(panel$chr,panel$pos, sep=":")
panel <- panel[,c("SNPID","rsid")]

for(trait in traitlist){

datasetpath <- paste0("../datasets/",trait, "-hm3.tsv.gz")
ds <- fread(datasetpath)

# This approach will require to have ALT_FREQ aligned to ALT - which is not required in other approaches. However, since we aligned datasets to HM3 dataset, and aligner doesn't reverse ALT_FREQ, we'll quickly get allele frequencies for them
ds[,ALT_FREQ:=NULL]
fwrite(ds, paste0("../datasets/",trait,"_SBayesR-nofreqs.tsv.gz"), sep="\t")
system(paste0("Compute_freqs.sh -f ../datasets/",trait,"_SBayesR-nofreqs.tsv.gz -p CEU"))

ds <- fread(paste0("../datasets/",trait,"_SBayesR-withfreqs.tsv.gz"))
ds[,SNPID:=paste(CHR,BP,sep=":")]
ds <- merge(ds,panel)
# Add N
if(trait == "MDD"){
	ds[,N:=N0+N1]
	ds <- ds[ N > quantile(N, 0.1),] # Remove SNPs in the lower 10% quantile, as recommended by SBayesR
} else if(trait == "BMI" || trait == "Height"){
       ds <- ds[ N > quantile(N, 0.1),]

} else{
       ds[,N:=gwasn[match(trait,traitlist)]]
}

	

# Select columns for SBayesR input
ds  <- ds[, .(rsid, ALT, REF, ALT_FREQ, BETA, SE, P, N)]
names(ds) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N" )

fwrite(ds, paste0("../datasets/",trait,"_SBayesR.ma"), sep = " ")

}
