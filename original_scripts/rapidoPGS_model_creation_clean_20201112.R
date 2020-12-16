# This script will create RápidoPGS models for 8 case-control and 2 quantitative trait datasets using single and multiple causal variant approach in rapidoPGS.
# 2020/09/21
# NOTE: This is a clean version of the original code used to run the models for illustrative purposes, although they should be equivalent.
# NOTE:This script will assign sd.prior=0.2 to rapidopgs_multi. Auto prior versions are commented
# Note that the functions below are experimental, which will be latter to RapidoPGS R package



# Let's load the libraries
# install.packages("RapidoPGS")
library(RapidoPGS)
library(data.table)
library(GenomicRanges)
# remotes::install_github("stephenslab/susieR")
library(susieR)
library(bigsnpr)

#############################################
########   Relevant functions   #############
#############################################

# These functions are the same as in the package, so let's see if they work ok - they should!

# ##' Posterior inclusion probabilities under a multiple causal variant model
# ##'
# ##' A wrapper around susie_RSS from the susieR package. See
# ##' \url{https://stephenslab.github.io/susieR} for the package and Wang et al.
# ##' JRSSB 2020 for the paper describing the model
# ##' \url{https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12388}
# ##'
# ##' NB: intended as an internal function, called from other functions in the
# ##' RapidoPGS package.
# ##' @param d data.table with columns SNPID, BETA, SE (other columns allowed, d
# ##'   will be ignored)
# ##' @param LD LD matrix with dimnames matching "SNPID" column
# ##' @param nref number of individuals used to create the LD matrix
# ##' @param pi_i prior probability that a given variant is causal. Default 1e
# ##' @param sd.prior Standard deviation prior of the trait, if NULL (default)  will be estimated
# ##' @param max_it maximum number of iterations for susie. The default of 100
# ##'   should be plenty, but increase if there is not convergence
# ##' @inheritParams rapidopgs_multi
# ##' @importFrom susieR susie_rss
# ##' @return returns pip for each SNP, which we will use as a snp weight in
# ##'   generating the PGS (=pip * BETA)
# ##' @examples
# ##' data(michailidou) # load example data
# ##' d=michailidou[hm_chrom==3 & abs(hm_pos-27303612) < 1e+5] # focus on a window of association
# ##' setnames(d, old = c("hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele",
# ##'   "hm_effect_allele", "hm_beta", "hm_effect_allele_frequency",
# ##'   "standard_error", "p_value"), new=c("SNPID","CHR", "BP","REF", "ALT",
# ##'   "BETA", "ALT_FREQ", "SE", "P")) # rename
# ##' LD=(1 - abs(outer(d$ALT_FREQ,d$ALT_FREQ, "-"))) *
# ##'    outer(sign(d$BETA), sign(d$BETA), "*")
# ##' dimnames(LD)=list(d$SNPID,d$SNPID)
# ##' susie_pip(d, LD)
# ##' @author Chris Wallace, Guillermo Reales

susie_pip <- function(d,LD,nref=503,pi_i=1e-4, sd.prior=NULL, max_it=100) {
  ## checks
  nd <- names(d)
  if(!all(c("SNPID","BETA","SE") %in% nd))
    stop("require columns SNP, BETA, SE")
  ## snps should be unique
  if("SNPID" %in% nd && is.factor(d$SNPID))
    stop("dataset ",suffix,": SNPID should be a character vector but is a factor")
  if("SNPID" %in% nd && any(duplicated(d$SNPID)))
    stop("dataset ",suffix,": duplicated SNPIDs found")
  ## SNPIDs should be in LD matrix
  if(!is.matrix(LD) ||
     !all(d$SNPID %in% colnames(LD)) ||
     !identical(rownames(LD),colnames(LD)) ||
     !all(LD >= -1 && LD <= 1))
    stop("LD should be a correlation matrix containing all SNPIDs listed in d$SNPID (match by rownames, colnames)")
  
  z=d$BETA/d$SE
  if(is.null(sd.prior)){
    res=susie_rss(z, LD[d$SNPID,d$SNPID], z_ld_weight = 1/nref,
                  ,null_weight=1 - length(d$SNPID)*pi_i
                  ,estimate_prior_method="EM"
                  ## ,verbose=TRUE
                  ,max_iter=max_it )
  } else{
    res=susie_rss(z, LD[d$SNPID,d$SNPID], z_ld_weight = 1/nref,
                  ,null_weight=1 - length(d$SNPID)*pi_i
                  ,estimate_prior_method="EM"
                  ,estimate_prior_variance=FALSE  ### <- stop internal estimation
                  ,prior_variance = sd.prior^2 ### <- use custom sd.prior
                  ## ,verbose=TRUE
                  ,max_iter=max_it )
  }
  pip <- res$pip[ -length(res$pip) ]
  names(pip) <- d$SNPID
  pip
}

# ##' Compute PGS from GWAS summary statistics using Bayesian sum of single-effect 
# ##' (SuSiE) linear regression using z scores
# ##' 
# ##' '\code{rapidopgs_mult} computes PGS from a from GWAS summary statistics 
# ##' using Bayesian sum of single-effect (SuSiE) linear regression using z scores
# ##' 
# ##' This function will take a GWAS summary statistic dataset as an input,
# ##' will assign LD blocks to it, then use a preset reference panel in Plink format 
# ##' to compute LD matrices for each block. Then SuSiE method will be used to 
# ##' compute posterior probabilities of variants to be causal and generate PGS
# ##' weights by multiplying those posteriors by effect sizes (\eqn{\beta}). 
# ##' Unlike \code{rapidopgs.single}, this approach will assume one or more causal variants
# ##' Also, \code{rapidopgs.multi} does not require allele frequencies nor sample sizes
# ##' 
# ##' The GWAS summary statistics file to compute PGS using our method must contain
# ##' the following minimum columns, with these exact column names:
# ##' \describe{
# ##'   \item{CHR}{Chromosome}
# ##'   \item{BP}{Base position (in GRCh37/hg19).
# ##'   \item{REF}{Reference, or non-effect allele}
# ##'   \item{ALT}{Alternative, or effect allele, the one \eqn{\beta} refers to}
# ##'   \item{BETA}{\eqn{\beta} (or log(OR)), or effect sizes}
# ##'   \item{SE}{standard error of \eqn{\beta}}
# ##'   \item{P}{P-value for the association test}
# ##' }
# ##' Other columns are allowed, and will be ignored.
# ##' 
# ##' Reference panel should be divided by chromosome, in Plink format.
# ##' Both reference panel and summary statistic dataset should be in GRCh37/hg19. 
# ##' 
# ##' @param data a data.table containing GWAS summary statistic dataset
# ##'   with all required information.
# ##' @param reference a string representing the path to the directory containing 
# ##'   the reference panel (eg. "../ref-data/").
# ##' @param ancestry a string indicating the ancestral population (DEFAULT: "EUR")
# ##' @param pi_i a scalar representing the prior probability (DEFAULT:
# ##'   \eqn{1 \times 10^{-4}}).
# ##' @param iterations number of maximum iterations at SuSie step. DEFAULT:100, 
# ##'    which should be plenty to ensure convergence. 
# ##' @param ncores a numeric specifying the number of cores to be used.
# ##'    If not set, it won't use parallelisation.
# ##' @param alpha.block a numeric threshold for minimum P-value in LD blocks.
# ##'    Blocks with minimum P above \code{alpha.block} will be skipped. Default: 1e-4.
# ##' @param alpha.snp a numeric threshold for P-value pruning within LD block.
# ##'    SNPs with P above \{alpha.snp} will be removed. Default: 0.01.
# ##' @return a data.table containing the sumstats dataset with
# ##'   computed PGS weights.
# ##' @import data.table 
# ##' @importFrom bigsnpr snp_match snp_cor snp_readBed snp_attach
# ##' @importFrom GenomicRanges GRanges findOverlaps
# ##' @importFrom IRanges IRanges
# ##' @export
# ##' @author Guillermo Reales, Chris Wallace
# ##' @examples
# ##' \dontrun{
# ##' sumstats <- data.table(CHR=c("4","20","14","2","4","6","6","21","13"), 
# ##'			BP=c(1479959, 13000913, 29107209, 203573414, 57331393, 11003529, 149256398, 
# ##'					25630085, 79166661), 
# ##'			REF=c("C","C","C","T","G","C","C","G","T"), 
# ##'			ALT=c("A","T","T","A","A","A","T","A","C"), 
# ##'			ALT_FREQ=c(0.2611,0.4482,0.0321,0.0538,0.574,0.0174,0.0084,0.0304,0.7528),
# ##'			BETA=c(0.012,0.0079,0.0224,0.0033,0.0153,0.058,0.0742,0.001,-0.0131),
# ##'			SE=c(0.0099,0.0066,0.0203,0.0171,0.0063,0.0255,0.043,0.0188,0.0074),
# ##'			P=c(0.2237,0.2316,0.2682,0.8477,0.01473,0.02298,0.08472,0.9573,0.07535))
# ##' PGS  <- rapidopgs_multi(sumstats, reference = "ref-data/")
# ##'}

rapidopgs_multi <- function(data, reference, ancestry="EUR", pi_i = 1e-04, iterations = 100, ncores=1, alpha.block=1e-4, alpha.snp=0.01, sd.prior=NULL){

  mincol <- c("CHR","BP", "REF","ALT","BETA", "SE","P")
  if(!all(mincol %in% names(data)))
    stop("All minimum columns should be present in the summary statistics dataset. Please check, missing columns: ",
         paste(setdiff(mincol,names(data)), collapse=", "))
  if(!file.exists(paste(reference,"chr1.bed", sep=""))) # Temp fix to detect panel		
    stop("No reference panel detected. Please check.")
  
  message("Running RápidoPGS with multiple causal variant assumption.")
  
  ds <- copy(data) # avoid modifying input data.table
  ds[,SNPID:=paste(CHR,BP, sep = ":")]
  ds  <- na.omit(ds, cols=c("CHR","BP", "BETA", "SE", "P"))
  # Assign ld.blocks, in case they werent computed yet.
  # Note that in this case only hg19 is admitted.
  if(!"ld.block" %in% names(ds)){
    blranges <- RapidoPGS::EUR_ld.blocks	
    message("Assigning LD blocks...")
    snpranges <- GRanges(seqnames=paste("chr",ds$CHR, sep=""), ranges=IRanges(start=ds$BP, end=ds$BP, names=ds$SNPID), strand="*")
    ds[,ld.block:=findOverlaps(snpranges, blranges, select='last')]
    message("Done!")
  }
  ## We ensure that rows without assigned block are removed
  if(any(is.na(ds$ld.block))){
    # warning(length(ds$ld.block[is.na(ds$ld.block)]), " SNPs didn't have any LD block assigned and will be removed")
    ds <- na.omit(ds, cols="ld.block")
  }
  
  # We don't need all columns, so we can filter them by name
  ds <- ds[,c("CHR","BP","SNPID","REF","ALT","BETA","SE","P","n_eff", "ld.block")]
  # Dear snp_match require specific names, so let's abide
  names(ds)[c(1:2,4:7)]  <- c("chr","pos","a0","a1","beta", "beta_se")
  results <- data.table()
  
  for(chrs in 1:22){
    message("Working on chromosome ", chrs,".")
    # First thing is to check if we already have .rds files for our chr (ie. if we have read the bed files before). If not, we'll read it. This will create a .bk file, but it won't be a problem since we won't call this function again.
    if(!file.exists(paste0(reference,"chr",chrs,".rds"))){
      snp_readBed(paste0(reference,"chr",chrs,".bed"))
    }
    # Attach the "bigSNP" object in R session
    # This object contain SNP data from a single chromosome from 1000 genomes phase 1.
    # This object contains a matrix of genotypes (analogous to bed), as well as fam and map slots, PLINK format.
    obj.bigSNP <- snp_attach(paste(reference,"chr",chrs, ".rds", sep=""))
    # In this case, we'll focus only on individuals of european ancestry
    euridx  <- grep(ancestry, obj.bigSNP$fam$family.ID)
    # Get aliases for useful slots
    G   <- obj.bigSNP$genotypes
    # Filter sumstats and panel by the SNPs that we're going to use
    ds.chr <- as.data.frame(ds[ds$chr == chrs,])
    map.chr <- obj.bigSNP$map[-3]
    names(map.chr) <- c("chr", "SNPID", "pos", "a0", "a1")
    map.chr$SNPID <- paste(map.chr$chr,map.chr$pos, sep=":") # We're using CHR:BP to match, rather than rsIDs
    
    map.chr <- map.chr[map.chr$SNPID %in% ds.chr$SNPID,]	
    message("Matching and aligning SNPs in chr",chrs," to the reference")
    # Align sumstats to panel
    snp.chr <- snp_match(ds.chr, map.chr)
    
    pb <- txtProgressBar(min = 0, max = length(unique(snp.chr$ld.block)), style = 3) # Show a nice progress bar
    
    for(block in unique(snp.chr$ld.block)){
      #			message("This is block ",block)
      
      idxbar <- which(unique(snp.chr$ld.block) == block)
      
      snp.block <- snp.chr[snp.chr$ld.block == block,]
      # Skip block if min P is above a certain threshold
      if(min(snp.block$P) > alpha.block){
        #				message ("\nBlock ", block," has no SNP with P-value below ",alpha.block," threshold. Skipping block.")
        setTxtProgressBar(pb,idxbar)
        next
      }
      # Remove SNPs with P above a certain threshold
      snp.block <- snp.block[snp.block$P < alpha.snp,]
      
      # Skip block if has only one SNP after filtering
      if(nrow(snp.block) < 2){
        #				message ("\nWarning: Block ", block," has only one SNP. Skipping...")
        setTxtProgressBar(pb,idxbar)
        next
      }
      
      # Recover SNP indices	
      snp.idx <- which(paste(obj.bigSNP$map$chromosome,obj.bigSNP$map$physical.pos, sep=":") %in% snp.block$SNPID) 
      # Remove duplicates, which sometimes appear
      if(length(snp.idx) != length(snp.block$SNPID)){
        du <- paste(obj.bigSNP$map$chromosome,obj.bigSNP$map$physical.pos, sep=":") 
        dup <- du[du %in% snp.block$SNPID]
        dup <- dup[duplicated(dup)]
        snp.block <- snp.block[!snp.block$SNPID %in% dup,]
        snp.idx <- which(du %in% snp.block$SNPID)
      }

      # Compute LD matrix for that block
      LDmat.block <- snp_cor(G, ind.col = snp.idx, ind.row= euridx, ncores = ncores, size = length(snp.idx))
      dimnames(LDmat.block) <- list(snp.block$SNPID,snp.block$SNPID)
      LDmat.block <- as.matrix(LDmat.block)	
      
      snp.block <- snp.block[,c("chr","pos", "a0","a1","SNPID","beta","beta_se")]
      names(snp.block) <- c("CHR","BP","REF","ALT","SNPID","BETA","SE")
      # Compute Susie
      snp.block$ppi_susie <- susie_pip(snp.block,LDmat.block,nref=length(euridx), pi_i=pi_i, max_it = iterations, sd.prior = sd.prior)
      # Append to results
      results <- rbind(results, snp.block)
      
      # Progress bar
      
      setTxtProgressBar(pb, idxbar)
    }
    close(pb)
  }
  results[,weight:=BETA*ppi_susie]
  return(results)
  
}


###############################################
#### GENERATING THE PGS MODELS   ##############
###############################################


# Set path to the datasets and reference, and their names
dspath <- "../datasets/"
datasets <- c("AST_Demenais_29273806_1-qcfilt.tsv.gz","RA_Okada_24390342_3-qcfilt.tsv.gz", "T1D_Cooper_doi101101120022_1-qcfilt.tsv.gz", "T2D_Scott_28566273_1-qcfilt.tsv.gz", "BRCA_Michailidou_29059683_1-qcfilt.tsv.gz","PRCA_Schumacher_29892016_1-qcfilt.tsv.gz","CAD_Nikpay_26343387_1-qcfilt.tsv.gz","MDD_Wray_29700475_1-qcfilt.tsv.gz")
# RapidoPGS-single needs sample sizes, so we provide them here
N0 <- c(107715,61565, 8828,132532,119078,61106,123504,113154)
N1 <- c(19954,19234,5913,26676,137045,79148,60801,59851)
# Note, these sample sizes don't include MDD, which had per-variant sample size, so we do a little trick here (see below)
ref  <- "../references/RapidoPGS-ref/"
ncores <- bigparallelr::nb_cores()

###################################################################
##### RÁPIDOPGS-MULTI MODEL CREATION (CASE-CONTROL DATASETS) ######
###################################################################


# Here we apply RápidoPGSmult to case-control datasets, with alpha.block = 1e-3 and alpha.snp=0.1
# NOTE: Unlike in previous versions, here we'll compute an extra layer of models for binary traits, setting the sd.prior to 0.2. 

sapply(datasets, function(x){
	filebasename <- strsplit(x, split="-")[[1]][1]
	ds <- fread(paste0(dspath,x))
# Automatically computed prior
#	fullresauto  <- rapidopgs_multi(ds,reference=ref, ncores=ncores,iterations=1000, alpha.block = 1e-3, alpha.snp = 0.1)
#	fwrite(fullresauto, paste("../models/",filebasename, "_qcfilt_RapidoPGSmulti1e301_autoprior.full.model",sep=""), sep="\t")
	fullres02  <- rapidopgs_multi(ds,reference=ref, ncores=ncores,iterations=1000, alpha.block = 1e-3, alpha.snp = 0.1, sd.prior=0.2)
	fwrite(fullres02, paste("../models/",filebasename, "_qcfilt_RapidoPGSmulti1e301_prior02.full.model",sep=""), sep="\t")
	    })


# Here we apply RápidoPGSmult to case-control datasets, with alpha.block = 1e-4 and alpha.snp=0.01
# NOTE: Unlike in previous versions, here we'll compute an extra layer of models for binary traits, setting the sd.prior to 0.2. 

sapply(datasets, function(x){
	filebasename <- strsplit(x, split="-")[[1]][1]
	ds <- fread(paste0(dspath,x))
#	fullresauto  <- rapidopgs_multi(ds,reference=ref, ncores=ncores,iterations=1000, alpha.block = 1e-4, alpha.snp = 0.01)
#	fwrite(fullresauto, paste("../models/",filebasename, "_qcfilt_RapidoPGSmulti1e4001_autoprior.full.model",sep=""), sep="\t")
	fullres02  <- rapidopgs_multi(ds,reference=ref, ncores=ncores,iterations=1000, alpha.block = 1e-4, alpha.snp = 0.01, sd.prior=0.2)
	fwrite(fullres02, paste("../models/",filebasename, "_qcfilt_RapidoPGSmulti1e4001_prior02.full.model",sep=""), sep="\t")
	    })

####################################################################
##### RÁPIDOPGS-SINGLE MODEL CREATION (CASE-CONTROL DATASETS) ######
####################################################################

# Here we apply RápidoPGS-single (aka Wakefield approach) to all case-control datasets, specifying their control (N0) and case (N1) number.
# For case-control datasets, default sd.prior = 0.2 will apply


sapply(seq_along(datasets), function(x){
	log.p <- FALSE # Usual default for RápidoPGS-single
	filebasename <- strsplit(datasets[x], split="-")[[1]][1]
	ds <- fread(paste0(dspath,datasets[x]))
	# Set N0 and N1 to column names for MDD dataset
	if(x == 8){ 
		N_controls <- "N0"
		N_cases <- "N1"
	}else{
		N_controls <- N0[x]
		N_cases <- N1[x]
	}
	if(x == 3){ # For T1D we'll use log.p
		log.p <- TRUE
	}
	wakefield.model  <- computePGS(ds,N0=N_controls,N1=N_cases, log.p=log.p) # We use RapidoPGS current function, which is equivalent to rapidopgs_single
	fwrite(wakefield.model, paste0("../models/",filebasename, "_qcfilt_RapidoPGSsingle_prior02.full.model"), sep="\t")

	    })

###################################################################
##### RÁPIDOPGS MODEL CREATION (QUANTITATIVE TRAIT DATASETS) ######
###################################################################

# The following part is intended to run RápidoPGS for BMI and height with custom priors, obtained from the computation of var_b.

priors <- c(0.08245557,0.09135172) # These are priors for var_b, see "Compute_custom_priors_BMI_Height_clean_20201112.R" to find out how they were calculated
datasets <- c("BMI_Locke_25673413_1-qcfilt.tsv.gz","HEIGHT_Wood_25282103_1-qcfilt.tsv.gz")

sapply(seq_along(datasets), function(x){
	filebasename <- strsplit(datasets[x], split="-")[[1]][1]
	prior <- priors[x]
	ds <- fread(paste0(dspath,datasets[x]))
	
	susie.model  <- rapidopgs_multi(ds,reference=ref, ncores=ncores,iterations=1000, alpha.block = 1e-3, alpha.snp = 0.1, sd.prior=prior)
	fwrite(susie.model, paste("../models/",filebasename, "_qcfilt_RapidoPGSmulti1e301_customprior.full.model",sep=""), sep="\t")
	
	N0 <- "N"
	wakefield.model <- computePGS(ds,N0, sd.prior = prior)
	fwrite(wakefield.model, paste("../models/",filebasename, "_qcfilt_RapidoPGSsingle_customprior.full.model",sep=""), sep="\t")
	    })


