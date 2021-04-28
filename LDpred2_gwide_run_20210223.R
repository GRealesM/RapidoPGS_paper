# LDpred2-auto genomewide
# 2021/02/24
# This script will attempt to run LDpred2 in a genowe-wide fashion, rather than by chr

# Load libraries
library(data.table)
library(bigsnpr)
setDTthreads(0)
set.seed(1)


ldpred2.auto.gwide <- function(sumstats, iterations=3000, initial_p=seq_log(1e-4,0.9,15), cores){

	# Load map, which will help us trace indices in sumstats to those in matrices
	map  <- as.data.table(readRDS("../references/UKBB_mat/map.rds"))
	map[,SNPID:=paste(chr,pos, sep=":")]

	# Create temporary directory to store matrices
        tmp <- tempfile(tmpdir = "tmp-data")
        on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
	
	# Prepare sumstats
	setnames(sumstats, c("BETA", "SE"), c("beta", "beta_se")) 
#	if(!"SNPID" %in% names(sumstats))
		sumstats[,SNPID:=paste(CHR,BP,sep=":")]
	# Check if dataset is aligned to map, and align if not
	fullid.map <- paste(map$chr, map$pos, map$a0, map$a1, sep =":")
	fullid.sumstats <- paste(sumstats$CHR, sumstats$BP, sumstats$REF, sumstats$ALT, sep=":")
	if(!all(fullid.sumstats %in% fullid.map))
		stop("Sumstats not aligned to map. Please check.")
	# General index map-sumstats
	map.idx <- match(fullid.map, fullid.sumstats)
	

	# Load matrices and match SNPs
	message("Preparing dataset and LD matrix...")

	for (chrs in 1:22){
	    message("Preparing LD matrix for chr ", chrs)	  
	    sumstats.chr <- sumstats[CHR == chrs,]
	    map.chr <- map[chr == chrs,]
	    idx.chr  <- match(sumstats.chr$SNPID, map.chr$SNPID)

	    # Load matrix
	    corr0 <- readRDS(paste0("../references/UKBB_mat/LD_chr",chrs,".rds"))[idx.chr, idx.chr]
	    
	    # Add up dataset and matrices
	    if (chrs == 1) {
		      df_beta <- sumstats.chr[, c("beta", "beta_se", "n_eff")]
      		      ld <- Matrix::colSums(corr0^2)
		      corr <- as_SFBM(corr0, tmp)
	    } else {
   		      df_beta <- rbind(df_beta, sumstats.chr[, c("beta", "beta_se", "n_eff")])
      		      ld <- c(ld, Matrix::colSums(corr0^2))
      		      corr$add_columns(corr0, nrow(corr))
    		   }
  	}

    
      ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff, blocks = NULL))
      h2_est <- ldsc[["h2"]]
      message("Computing LDpred2 auto, with initial p = ",paste(initial_p, collapse=", "), ", initial h2 = ", h2_est,", and ",iterations," iterations (+1000 burn-in)")
     
      multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                                 vec_p_init = initial_p, num_iter=iterations,
                                 ncores = cores)
      beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
      sumstats$weight <- rowMeans(beta_auto)
      setnames(sumstats, c("beta","beta_se"), c("BETA","SE"))
      return(sumstats)

}



#######################

# Here we're using an array job to run the models, so we'll create an array in bash containing
# the indexes to select the traits.
# And we'll use those to generate the models
args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args) # Array index will be a 1-10 integer.

traitlist <- c("Asthma","RA","T1D","T2D","BRCA","PRCA","CAD","MDD","BMI","Height")
args <- traitlist[args] # This way we transform the index number into the real trait name

dataset <- paste0(args, "-hm3.tsv.gz")
cores <- bigparallelr::nb_cores() # Use all cores available (15 in our case)

ds <- fread(paste0("../datasets/",dataset))

ldpred2.model <- ldpred2.auto.gwide(ds, iterations=3000, initial_p=seq_log(1e-4,0.9,15), cores=cores)
fwrite(ldpred2.model, paste0("../models/",args, "_hm3_LDpred2auto_gwide.full.model"), sep="\t")

######################


#### WORKSHOP #######
# Load libraries
#library(data.table)
#library(bigsnpr)
#set.seed(1)
#
#chrs=22
#cores <- 8
#
#sumstats <- fread("../datasets/Asthma-hm3.tsv.gz")
#
#map  <- as.data.table(readRDS("../references/UKBB_mat/map.rds"))
#map[,SNPID:=paste(chr,pos, sep=":")]
## Create temporary directory to store matrices
#tmp <- tempfile(tmpdir = "tmp-data")
## on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
#	
#	# Prepare sumstats
#	setnames(sumstats, c("BETA", "SE"), c("beta", "beta_se")) 
#	if(!"SNPID" %in% names(sumstats))
#		sumstats[,SNPID:=paste(CHR,BP,sep=":")]
#	# Check if dataset is aligned to map, and align if not
#	fullid.map <- paste(map$chr, map$pos, map$a0, map$a1, sep =":")
#	fullid.sumstats <- paste(sumstats$CHR, sumstats$BP, sumstats$REF, sumstats$ALT, sep=":")
#	if(!all(fullid.sumstats %in% fullid.map))
#		stop("Sumstats not aligned to map. Please check.")
#	# General index map-sumstats
#	map.idx <- match(fullid.map, fullid.sumstats)
#	    
#	sumstats.chr <- sumstats[CHR == chrs,]
#        map.chr <- map[map$chr == chrs,]
#        idx.chr  <- match(sumstats.chr$SNPID, map.chr$SNPID)
#
#	    # Load matrix
#       corr0 <- readRDS(paste0("../references/UKBB_mat/LD_chr",chrs,".rds"))[idx.chr, idx.chr]
# 	    
#	    # Add up dataset and matrices
#	    if (chr == 1) {
#		      df_beta <- sumstats.chr[, c("beta", "beta_se", "n_eff")]
#      		      ld <- Matrix::colSums(corr0^2)
#		      corr <- as_SFBM(corr0, tmp)
#	    } else {
#   		      df_beta <- rbind(df_beta, sumstats.chr[, c("beta", "beta_se", "n_eff")])
#      		      ld <- c(ld, Matrix::colSums(corr0^2))
#      		      corr$add_columns(corr0, nrow(corr))
#    		   }
#  	}
#
    

