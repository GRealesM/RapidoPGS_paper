# This script will create RápidoPGS models for 8 case-control and 2 quantitative trait datasets using single and multiple causal variant approach in rapidoPGS.
# 2021/03/17
# Note that the functions below are experimental, and will be latter to RapidoPGS R package

# Let's load the libraries
# remotes::install.packages("GRealesM/RapidoPGS") # v2.1.0 commit fd8a09c05e4b1518068b0d4537008d36cee34c02
library(RapidoPGS)
#library(dats.table)
#library(GenomicRanges)
## remotes::install_github("stephenslab/susieR") # v0.10.1
#library(susieR)
#library(bigsnpr)
## remotes::install_github("chr1swallace/coloc", ref="susie") #4.0-6
#library(coloc)

###############################################
#### GENERATING THE PGS MODELS   ##############
###############################################

# Here we're using an array job to run the models, so we'll create an array in bash containing
# the indexes to select the traits.
# And we'll use those to generate the models
args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args) # Array index will be a 1-10 integer.

traitlist <- c("Asthma","RA","T1D","T2D","BRCA","PRCA","CAD","MDD","BMI","Height", "RAint")
args <- traitlist[args] # This way we transform the index number into the real trait name

dataset <- paste0(args, "-hm3.tsv.gz")
ncores <- 15 # Use all cores available (15 in our case)
setDTthreads(15) # Give data.table 15 cores too

ds <- fread(paste0("../datasets/",dataset))
N_LD=362320 # N individuals to generate UKBB LD matrices, provided by Privé et al.
LDmatrices <- "../references/UKBB_mat"



# Check if the datasets are case-controls or quantitative, and adapt the parameters to it
if(args == "BMI" || args == "Height"){
 trait.type <- "quant"
 N <- "N" 
 # We set the sd priors for quantitative traits, for both RápidoPGS single and multi.
 if(args == "BMI"){
   priors  <- 0.1272537 
 } else if(args == "Height"){
   priors  <- 0.1415274
 }
} else {
  trait.type <- "cc"
  N  <-  NULL
  priors <- 0.2 # For RápidoPGS single and multiple (set prior)
}


###############################################
##### RÁPIDOPGS-SINGLE MODEL CREATION #########
###############################################

wakefield.model  <- rapidopgs_single(ds, N=N, trait=trait.type, sd.prior=priors) # We use RapidoPGS current function, which is equivalent to rapidopgs_single
fwrite(wakefield.model, paste0("../models/",args, "_hm3_RapidoPGSsingle.full.model"), sep="\t")
# By default, cc will have sd.prior == 0.2 and quantitative traits will have previously estimated sd.prior (see above).

#####################################################
####### RÁPIDOPGS-MULTI (SET PRIOR) MODEL CREATION ##
#####################################################

# Here we apply RápidoPGS multi, with alpha.block = 1e-3 and alpha.snp=0.1

multi1e301  <- rapidopgs_multi(ds, trait=trait.type, LDmatrices=LDmatrices, N_LD = N_LD, N=N, ncores=ncores, alpha.block = 1e-3, alpha.snp = 0.1, sd.prior=priors)
fwrite(multi1e301, paste0("../models/",args, "_hm3_RapidoPGSmulti_1e301_setprior.full.model"), sep="\t")

# Here we apply RápidoPGS multi, with alpha.block = 1e-4 and alpha.snp=0.01

multi1e4001  <- rapidopgs_multi(ds, trait=trait.type, LDmatrices=LDmatrices, N_LD = N_LD, N=N, ncores=ncores, alpha.block = 1e-4, alpha.snp = 0.01, sd.prior=priors)
fwrite(multi1e4001, paste0("../models/",args, "_hm3_RapidoPGSmulti_1e4001_setprior.full.model"), sep="\t")

#######################################################
######## RÁPIDOPGS-MULTI (AUTO PRIOR) MODEL CREATION ##
#######################################################

# Here we apply RápidoPGS multi to case-control datasets, with alpha.block = 1e-3 and alpha.snp=0.1

multi1e301  <- rapidopgs_multi(ds, trait=trait.type, LDmatrices=LDmatrices, N_LD = N_LD, N=N, ncores=ncores, alpha.block = 1e-3, alpha.snp = 0.1)
fwrite(multi1e301, paste0("../models/",args, "_hm3_RapidoPGSmulti_1e301_autoprior.full.model"), sep="\t")

# Here we apply RápidoPGS multi, with alpha.block = 1e-4 and alpha.snp=0.01

multi1e4001  <- rapidopgs_multi(ds, trait=trait.type, LDmatrices=LDmatrices, N_LD = N_LD, N=N, ncores=ncores, alpha.block = 1e-4, alpha.snp = 0.01)
fwrite(multi1e4001, paste0("../models/",args, "_hm3_RapidoPGSmulti_1e4001_autoprior.full.model"), sep="\t")

######## DEBUGGING ##########

#library(data.table)
#library(GenomicRanges)
## remotes::install_github("stephenslab/susieR")
#library(susieR)
#library(bigsnpr)
## remotes::install_github("chr1swallace/coloc", ref="susie")
#library(coloc)
#
#data <- fread("../datasets/Asthma-hm3.tsv.gz")
#trait = "cc"
#N_LD=362320
#chrs=3
#block=325
#build="hg19"
#LDmatrices="../references/UKBB_mat"
#alpha.snp=0.1
#sd.prior=0.2
#
## Preparing dataset 
#  if(!"data.table" %in% class(data))
#	  data <- as.data.table(data)
#  if(!trait %in% c("cc", "quant")) stop("Please, specify your study type, choose case-control ('cc') or quantitative ('quant').")
#  if(length(trait) != 1) stop("Please select only one study type")
#  if(trait == "quant" && is.null(N)) stop("N (sample size) is required for quantitative traits, please provide them, either as an integer or as column name containing it.")
#  if(trait == "quant"){
#  	mincol <- c("CHR","BP", "REF","ALT","BETA", "SE", "ALT_FREQ")
#	  if(!all(mincol %in% names(data)))
#    	stop("All minimum columns should be present in the summary statistics dataset. Please check, missing columns: ",
#         paste(setdiff(mincol,names(data)), collapse=", "))  
#  	ds <- copy(data) # avoid modifying input data.table
#	ds  <- na.omit(ds, cols=c("CHR","BP", "BETA", "SE", "ALT_FREQ")) # Remove NA in relevant columns
#  } else{
#  	## Here's a list of columns that the dataset must have, otherwise it will fail
# 	mincol <- c("CHR","BP", "REF","ALT","BETA", "SE")
#  	if(!all(mincol %in% names(data)))
#  	  stop("All minimum columns should be present in the summary statistics dataset. Please check, missing columns: ",
#  	       paste(setdiff(mincol,names(data)), collapse=", "))
#  	ds <- copy(data) # avoid modifying input data.table
#  	ds  <- na.omit(ds, cols=c("CHR","BP", "BETA", "SE")) # Remove NA in relevant columns
#	ds <- ds[SE != 0,] # Remove SNPs with exactly zero SE (sign of problems in the data)
#  }	
#  if(!"SNPID" %in% names(ds)){
#	ds[,SNPID:=paste(CHR,BP, sep=":")] # SNPID is not strictly required to be provided. If it's not, we create it using CHR:BP
#	       }
#   # Assign ld.blocks, in case they werent computed yet
#   if(!"ld.block" %in% names(ds)){
#   	if(build == "hg19"){ 
#   		blranges <- RapidoPGS::EUR_ld.blocks
#   	}else if(build == "hg38"){ 
#   		blranges <- RapidoPGS::EUR_ld.blocks38
#   	}else{ 
#   		stop("RapidoPGS only accepts hg19 or hg38 at the moment, please check.")
#   	}		
#   message("Assigning LD blocks...")
#   snpranges <- GRanges(seqnames=paste("chr",ds$CHR, sep=""), ranges=IRanges(start=ds$BP, end=ds$BP, names=ds$SNPID), strand="*")
#   ds[,ld.block:=findOverlaps(snpranges, blranges, select='last')]
#   message("Done!")
#   }
#  ds <- ds[!is.na(ld.block),] # Remove SNPs without assigned block
#
#  # We don't need all columns, so we can filter them by name
#  ds <- ds[,c("CHR","BP","SNPID","REF","ALT","BETA","SE","P","n_eff", "ld.block")]
#  # Dear snp_match require specific names, so let's abide
#  names(ds)[c(1:2,4:7)]  <- c("chr","pos","a0","a1","beta", "beta_se")
#
#    map <- as.data.table(readRDS(paste(LDmatrices, "map.rds", sep="/")))
#    map[,SNPID:=paste(chr,pos, sep=":")] 
#    # Check if ds is aligned to map, and align if not
#    ds.snp <- paste(ds$chr, ds$pos, ds$a0, ds$a1, sep=":")
#    map.snp <- paste(map$chr, map$pos, map$a0, map$a1, sep=":") 
#
#    if(!all(ds.snp %in% map.snp)){
#       ds <- snp_match(ds, map)
#    }
#
#
#    LD.chr <- readRDS(paste0(LDmatrices, "/LD_chr",chrs, ".rds"))
#    ds.chr  <- ds[chr == chrs,]
#    map.chr <- map[chr == chrs,]
#      
#      idxbar <- which(unique(ds.chr$ld.block) == block)
#      snp.block <- ds.chr[ds.chr$ld.block == block,]
#      
#      # Remove SNPs with P above a certain threshold
#      snp.block <- snp.block[snp.block$P < alpha.snp,]
#      
#
#      # Match ids of resulting SNPs with those in map manifest (and hence, LD matrix)
#      # Remove duplicates
#      if(all(!duplicated(snp.block$SNPID))){
#	snp.block <- snp.block[!duplicated(snp.block$SNPID),]
#	 }
#      
#      snp.idx <- match(snp.block$SNPID, map.chr$SNPID) 
#      LD.block <- as.matrix(LD.chr[snp.idx,snp.idx])
#      dimnames(LD.block) <- list(snp.block$SNPID, snp.block$SNPID)
#      snp.block <- snp.block[,c("chr","pos", "a0","a1","SNPID","beta","beta_se")]
##      names(snp.block) <- c("CHR","BP","REF","ALT","SNPID","BETA","SE")
#     # Prepare dataset for runsusie
#     names(snp.block) <- c("CHR","BP","REF","ALT","SNPID","BETA","SE")
#     susie.ds <- list(snp=snp.block$SNPID, beta=snp.block$BETA, varbeta=snp.block$SE^2, LD=LD.block, type=trait)
#     
#     if(is.null(sd.prior)){
#             prior_est = TRUE
#	     prior_var = NULL
#     } else{
#             prior_est = FALSE
#	     prior_var = sd.prior^2
#     }
#
#     ppi_susie <- runsusie(susie.ds,nref=N_LD, prior_variance=prior_var, estimate_prior_variance=prior_est, check_R=FALSE)
#     ppi_susie <- ppi_susie$pip[1:(length(ppi_susie$pip)-1)]
#     snp.block$pip_susie <- ppi_susie
#
#
#     saveRDS(snp.block, "../datasets/ds_block325.rds")
#     saveRDS(LD.block, "../datasets/LDmat_block325.rds")

######## WORKSHOP AREA ########
#
## Libraries
#library(GenomicRanges)
#
## Test dataset 
#data <- fread("../datasets/Asthma-hm3.tsv.gz")
#
## LD matrices path
#LDmatrices <- "../references/UKBB_mat"
#
#chrs=1
#block=1
#alpha.block=1e-4
#alpha.snp=0.01
#sd.prior=0.2
#pi_i = 1e-04
#iterations = 100
#ncores=12
#N_LD=362320
#
#model <- rapidopgs_multi(data, reference=NULL, LDmatrices="../references/UKBB_mat", N_LD = 362320, ancestry="EUR", pi_i = 1e-04, iterations = 100, ncores=8, alpha.block=1e-4, alpha.snp=0.01, sd.prior=0.2)
#
#
#
###############################

