#### Running PRScs-auto
## 2021/03/17

## This R script will prepare the datasets for PRScs and run PRScs-auto on then

library(data.table)

# Here we're using an array job to run the models, so we'll create an array in bash containing
# the indexes to select the traits.
# And we'll use those to generate the models
args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args) # Array index will be a 1-11 integer.

traitlist <- c("Asthma","RA","T1D","T2D","BRCA","PRCA","CAD","MDD","BMI","Height", "RAint")
trait <- traitlist[args] # This way we transform the index number into the real trait name

dataset <- paste0(trait, "-hm3.tsv.gz")
ds <- fread(paste0("../datasets/",dataset))
ds[,SNPID:=paste(CHR,BP, sep=":")]

panel <- fread("../references/PRScs/snpinfo_1kg_hm3")
panel[,SNPID:=paste(CHR,BP, sep=":")]
panel <- panel[,.(SNP,SNPID)]

ds <- merge(ds, panel)

# Create bim file
bim <- copy(ds)
bim[,geno:=0]
bim  <- bim[,.(CHR, SNP,  geno, BP, ALT, REF)]
fwrite(bim, paste0("../datasets/dataset_",args, ".bim"), col.names= FALSE, sep="\t")

# Prepare dataset for PRScs
ds <- ds[,.(SNP, ALT, REF, BETA, P)]
setnames(ds, c("ALT","REF"), c("A1","A2"))
fwrite(ds, paste0("../datasets/", args,"_PRScs.tsv"), sep="\t")
