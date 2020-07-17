# Benchmarking RapidoPGS vs. LDpred2

# Here we'll compare the timing for PGS generation using rapidoPGS and LDpred2 on a curated Asthma dataset

# Load required libraries
library(bigsnpr)
library(microbenchmark)
library(data.table)
library(GenomicRanges)

### RapidoPGS functions

# Load some helper functions
### g.complement, as found in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)
g.complement <- function (x) {
  x <- toupper(x)
  switches <- c(A = "t", T = "a", C = "g", G = "c")
  for (i in seq_along(switches)) x <- sub(names(switches)[i], switches[i], x)
  toupper(x)
}

### g.rev, as found in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)
g.rev <- function (x, sep = "/") {
  sapply(strsplit(x, sep), function(g) paste(rev(g), collapse = "/"))
}

### g.class, as found in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)
g.class <- function (x, y) {
  if (!identical(names(x), names(y))) 
    stop("x and y must relate to same SNPs")
  mat <- matrix(FALSE, length(x), 4, dimnames = list(names(x), c("nochange", "rev", "comp", "revcomp")))
  mat[, "nochange"] <- x == y
  mat[, "rev"] <- x == g.rev(y)
  mat[, "comp"] <- x == g.complement(y)
  mat[, "revcomp"] <- x == g.rev(g.complement(y))
  indels <- x %in% c("I/D", "D/I")
  if (any(indels)) 
    mat[indels, c("comp", "revcomp")] <- FALSE
  ret <- character(nrow(mat))
  rs <- rowSums(mat)
  if (length(wh <- which(rs > 1))) 
    ret[wh] <- "ambig"
  if (length(wh <- which(rs == 0))) 
    ret[wh] <- "impossible"
  if (length(wh <- which(rs == 1))) 
    ret[wh] <- colnames(mat)[apply(mat[wh, , drop = FALSE], 1, which)]
  return(ret)
}

# We need to define logsum in cupcake because for some reason it does not define it.
logsum <- function(x) {
    my.max <- max(x) ##take out the maximum value in log form)
    my.res <- my.max + log(sum(exp(x - my.max )))
    return(my.res)
}

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

wakefield_pp_quant <- function(beta, se.beta, sdY, sd.prior=0.15, pi_i=1e-4) { 
  # compute V
  V <- se.beta^2
  # Compute z too
  z <- beta/se.beta
  # Multiply prior by sdY
  prior <- sd.prior*sdY
  ## Shrinkage factor: ratio of the prior variance to the total variance
  r <- prior^2 / (prior^2 + V)
  ## Approximate BF
  lABF = 0.5 * (log(1-r) + (r * z^2))
  ## tABF - to add one we create another element at the end of zero for which pi_i is 1
  tABF <- c(lABF,0)
  vpi_i<-c(rep(pi_i,length(lABF)),1)
  sBF <- logsum(tABF + log(vpi_i))
  exp(lABF+log(pi_i)-sBF)
}

wakefield_pp <- function(p,f, N, s,pi_i=1e-4,sd.prior=0.2) {
    if(length(p) != length(f))
      stop("p and f must be vectors of the same size")
    # compute V
    V <- 1 / (2 * N * f * (1 - f) * s * (1 - s))
    # convert p vals to z
    z <- stats::qnorm(0.5 * p, lower.tail = FALSE)
    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF
    lABF = 0.5 * (log(1-r) + (r * z^2))
    ## tABF - to add one we create another element at the end of zero for which pi_i is 1
    tABF <- c(lABF,0)
    vpi_i<-c(rep(pi_i,length(lABF)),1)
    sBF <- logsum(tABF + log(vpi_i))
    exp(lABF+log(pi_i)-sBF)
}



# NOTE: We now included the reference set as argument, so it can be changed!
pgs.file.preprocess  <- function(dataset, blockfile="ld.blocks.RDS", ref=NULL){

	message("Loading dataset...")
	ds <- fread(dataset)
	# Here's a list of columns that the dataset must have, otherwise it will fail
	mincol <- c("CHR19","BP19", "REF","ALT","BETA", "SE", "P", "ALT_FREQ")
	if(!all(mincol %in% names(ds))) stop("All minimum columns should be present. Please check, missing columns: ", setdiff(mincol,names(ds)))
	
	if(!is.null(ref)){
	# Next step is to filter and align our alleles and their effects to the hapmap3 reference, which I have already formatted for our purposes.
	message("Filtering SNPs...")
	refset <- fread(ref)
	refset[, alleles:=paste(REF,ALT, sep="/")][,pid:=paste(CHR19, BP19, sep=":")]
	ds[,alleles:=paste(REF,ALT, sep="/")][,pid:=paste(CHR19, BP19, sep=":")]
	ds <- ds[pid %in% refset$pid,]
	ds  <- merge(ds, refset[,.(SNPID, pid, alleles)], by ='pid', suffixes=c("", ".reference"))
	ds[, alleles:=toupper(alleles)][, c("REF", "ALT"):=list(toupper(REF), toupper(ALT))][, SNPID:=SNPID.reference]

	message("Aligning alleles...")
	# Check if everything is alright  
	if(!all(g.class(ds$alleles.reference, ds$alleles)== "nochange")){
	    allele_diagnostics <- g.class(ds$alleles.reference, ds$alleles)
	    alleles_to_flip <-  allele_diagnostics == "rev"
	    alleles_to_comp <- allele_diagnostics == "comp"
	    alleles_to_revcomp <- allele_diagnostics == "revcomp"
	    cat("Some SNPs have to be flipped. ", sum(alleles_to_flip), " to flip, ", sum(alleles_to_comp), " to find their complement, and ", sum(alleles_to_revcomp), " to find their reverse complement.\n")
	    ds$alleles[alleles_to_flip] <- unlist(g.rev(ds$alleles[alleles_to_flip]))
	    ds$alleles[alleles_to_comp] <- g.complement(ds$alleles[alleles_to_comp])
	    ds$alleles[alleles_to_revcomp] <- unlist(g.rev(g.complement(ds$alleles[alleles_to_revcomp])))
	    ds$REF <- sapply(strsplit(ds$alleles, "/"), `[`, 1)
	    ds$ALT <- sapply(strsplit(ds$alleles, "/"), `[`, 2)
	    ds$BETA[alleles_to_flip] <- ds$BETA[alleles_to_flip]*-1
	    ds$BETA[alleles_to_revcomp] <- ds$BETA[alleles_to_revcomp]*-1
	   # NOTE: I introduced the following bit from milcytokine_basis on to guarantee that we have no ambiguity nor duplicated SNPs
	#if(!all(g.class(ds$alleles.reference, ds$alleles)== "nochange")){
	#	ds  <- ds[g.class(ds$alleles.reference, ds$alleles)== "nochange",]
	#	}
	
	    rm(alleles_to_flip, alleles_to_comp, alleles_to_revcomp)
	  }
	ds[, c("alleles", "alleles.reference"):=NULL]
	}

	ds <- unique(ds)

	# I made ld.blocks file from fourier_ls-all.bed, so now I can simply load the object
	# blocks <- fread("fourier_ls-all.bed")
	# blranges <- GRanges(seqnames=blocks$chr, ranges=IRanges(blocks$start, end=blocks$stop, names=1:nrow(blocks)), strand="*") 
	# saveRDS(blranges, "ld.blocks.RDS")
	message("Assigning LD blocks...")
	blranges <- readRDS(blockfile)
	ds  <- ds[!is.na(CHR19) & !is.na(BP19),] 
	snpranges <- GRanges(seqnames=paste("chr",ds$CHR19, sep=""), ranges=IRanges(ds$BP19, end=ds$BP19, names=ds$SNPID), strand="*")
        ds[,ld.block:=findOverlaps(snpranges, blranges, select='last')]
	ds  <- na.omit(ds)
	message("Done!")
	return(ds)
}

computePGS <- function(ds, N0,N1=NULL,pi_i= 1e-04, sd.prior=0.2, filt_threshold = NULL, recalc=FALSE, forsAUC=FALSE, altformat=FALSE){

	if(is.null(N1)){
		if(is.numeric(N0) && length(N0) == 1) { # In case N0 is supplied
		message("N1 not supplied. Assuming quantitative trait with ", N0, " individuals. Computing PGS.")
		ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=N0)][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]
	if(!is.null(filt_threshold)){
		ds  <- ds[ds$ppi > filt_threshold,]
		if(recalc){
		ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=N0)][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]		
		}
	}

		} else{ # In case column name is supplied
		if(is.character(N0) && length(N0) == 1){
		Nco <- N0
		message("N1 not supplied.  Assuming quantitative trait with multiple N, supplied by column ", Nco,". Mmax N: ", max(ds[,get(Nco)]), ", min N = ", min(ds[,get(Nco)]), ", and mean N: ",mean(ds[,get(Nco)]), ". Computing PGS.")
		ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=get(Nco))][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]
	if(!is.null(filt_threshold)){
		ds  <- ds[ds$ppi > filt_threshold,]
		if(recalc){
		ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=get(Nco))][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]		
		}

	}
		}
		}
	} else  { # If both are supplied
		if(is.character(N0) && is.character(N1)){
			message("Computing PGS for a case-control dataset. Both N0 and N1 columns provided.")
			Nco  <- N0
			Nca <- N1
			ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N=get(Nco)+get(Nca), s = get(Nca)/(get(Nco)+get(Nca)), pi_i, sd.prior), by = "ld.block"][, weight:=ppi*BETA]	
				if(!is.null(filt_threshold)){
					ds  <- ds[ds$ppi > filt_threshold,]
				if(recalc){
					ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N= get(Nco)+get(Nca) , s = get(Nca)/(get(Nca)+get(Nca)), pi_i = 1e-04, sd.prior=0.2), by = "ld.block"][, weight:=ppi*BETA]
				}
				}
		}else{
			message("Computing PGS for a case-control dataset, with ", N0," controls, and ", N1, " cases.")
			ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N= N0+N1, s = N1/(N0+N1), pi_i, sd.prior), by = "ld.block"][, weight:=ppi*BETA]	
				if(!is.null(filt_threshold)){
					ds  <- ds[ds$ppi > filt_threshold,]
				if(recalc){
					ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N= N0+N1 , s = N1/(N0+N1), pi_i = 1e-04, sd.prior=0.2), by = "ld.block"][, weight:=ppi*BETA]
				}
				}
		  }	
		}
	if(forsAUC){
		ds[,pid:=paste(CHR19,BP19,sep=":")]
		ds <- ds[,c("pid", "ALT", "weight")]
	}
	if(altformat){
		ds  <- ds[,c("SNPID", "ALT", "weight")]
	}

	return(ds)
}	

## LDPred2 function

computePGS.ldpred2.auto  <- function(sumstats, N0,N1=NULL, ancestry="EUR", iterations=500, initial_p=1e-4){
	
	message("Computing PGS using LDpred2 auto.")
	message("Loading summary statistics dataset: ", sumstats,".")
	NCORES <- nb_cores()
	sumstats <- bigreadr::fread2(sumstats)

	if(is.null(N1)){ # Only N0 is supplied
		message("No N1 supplied. Quantitative trait assumed")
		if(is.numeric(N0) && length(N0)==1){ # N0 is supplied as a number
			message("Sample size supplied (",N0,").")
			if(all(c("ALT_FREQ","SE") %in% names(sumstats))){
				sdy  <- sdY.est(sumstats$SE^2,sumstats$ALT_FREQ, N0)
				message("sdY for this dataset is ", sdy,".")
				sumstats$n_eff  <- 4 * N0 / sdy^2
			} else stop("ALT_FREQ or SE missing! Please check!")	
		} else if(is.character(N0) && length(N0) ==1){
			message("Sample size column supplied (",N0,").")
			N0  <- sumstats[,N0]
			if(all(c("ALT_FREQ","SE") %in% names(sumstats))){
				sdy  <- sdY.est(sumstats$SE^2,sumstats$ALT_FREQ, N0)
				message("sdY for this dataset is ", sdy,".")
				sumstats$n_eff  <- 4 * N0 / sdy^2
			} else stop("ALT_FREQ or SE missing! Please check!")	
		} else stop("Please provide a valid sample size for N0, either a column name or a atomic numeric")
	} else if(is.character(N0) && is.character(N1)){ # Both columns are supplied as column names
		message("Sample size column supplied, controls: ",N0,", and cases: ",N1,".")
		N0 <- sumstats[,N0]
		N1 <- sumstats[,N1]
		sumstats$n_eff <- 4 / (1 / N1 + 1 / N0)
	} else if(is.numeric(N0) && is.numeric(N1)){ # Both columns are supplied as numbers
		message("Sample size supplied. Controls:",N0,", Cases:",N1,".")
		sumstats$n_eff <- 4 / (1 / N1 + 1 / N0)
	} else { stop("Please provide valid sample sizes: Either N0 (quantitative traits) or N0 and N1 (case-control studies) as either atomic numeric or column names containing them in the dataset.")
	}

	# We don't need all columns, so we can filter them by name
	sumstats <- sumstats[,c("CHR19","BP19","SNPID","REF","ALT","BETA","SE","P","n_eff")]
	# Dear snp_match require specific names, so let's abide
	names(sumstats)[c(1:2,4:7)]  <- c("chr","pos","a0","a1","beta", "beta_se")

	# We don't want NA in our dataset, because this will give us problems down the road.
	# We'll also filter out sites with less than half the Neff, as LDpred2 guys do in their vignette. This will affect only files with columns for sample sizes, since for the rest n_eff will be equivalent for all sites.
	sumstats <- na.omit(sumstats)
	sumstats <- sumstats[sumstats$n_eff > (0.5 * max(sumstats$n_eff)),]

	results <- data.frame()
	for(chr in 1:22){
		
		message("Working on chromosome ", chr)
		snp_readBed(paste("../LDpred2/ref-data/chr",chr,".bed", sep=""))
		# Attach the "bigSNP" object in R session
		# This object contain SNP data from a single chromosome from 1000 genomes phase 1.
		# This object contains a matrix of genotypes (analogous to bed), as well as fam and map slots, PLINK format.
		obj.bigSNP <- snp_attach(paste("../LDpred2/ref-data/chr",chr, ".rds", sep=""))

		# See how the file looks like
		#str(obj.bigSNP, max.level = 2, strict.width = "cut")

		# In this case, we'll focus only on individuals of european ancestry
		euridx  <- grep(ancestry, obj.bigSNP$fam$family.ID)
		eurset  <- subset(obj.bigSNP, ind.row=euridx)
		obj.bigSNP <- snp_attach(eurset)

		# Get aliases for useful slots
		G   <- obj.bigSNP$genotypes
		CHR <- obj.bigSNP$map$chromosome
		POS <- obj.bigSNP$map$physical.pos
		#y   <- obj.bigSNP$fam$affection - 1
		#NCORES <- nb_cores()

		sumstats.chr <- sumstats[sumstats$chr == chr,]

		# Also, we must remove all positions with se = 0 (otherwise it will complain
		sumstats.chr <- sumstats.chr[sumstats.chr$beta_se > 0,]
		# Matching variants between genotype data and summary statistics
		# To match variants contained in genotype data and summary statistics, the variables "chr" (chromosome number), "pos" (genetic position), "a0" (reference allele) and "a1" (derived allele) should be available in the summary statistics and in the genotype data. These 4 variables are used to match variants between the two data frames.
		# This step is analogous to the alignment we do when reducing
		# We need to calculate the Neff from the summary statistics file to be used later by snp_ldsc2, which performs a LDscore regression. 
		#sumstats$n_eff <- 4 / (1 / sumstats$n_case + 1 / sumstats$n_control)
		#sumstats$n_case <- sumstats$n_control <- NULL
		map <- obj.bigSNP$map[-(2:3)]
		names(map) <- c("chr", "pos", "a0", "a1")

		message("Matching and aligning SNPs to the reference")
		info_snp.chr <- snp_match(sumstats.chr, map)
		
		message("Original sumstat had ", nrow(sumstats.chr)," for chr",chr,". After matching ", nrow(info_snp.chr)," remained, and will be used for further steps.")
		# Computing LDpred2 scores for one chromosome

		#Correlation

		# First, you need to compute correlations between variants. We recommend to use a window size of 3 cM (see ref).

		# snp_asGeneticPos 
		POS2 <- snp_asGeneticPos(CHR, POS, dir = "../LDpred2/ref-data", ncores = NCORES)

		## indices in info_snp
		#ind.chr <- which(info_snp$chr == 1)
		df_beta <- info_snp.chr[, c("beta", "beta_se", "n_eff")]

		## indices in G
		ind.chr2 <- info_snp.chr$`_NUM_ID_`
		message("Computing correlation matrix for chr",chr)
		
		corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = NCORES, infos.pos = POS2[ind.chr2], size = 3 / 1000)

		corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))

# We need to estimate h2 from the data, using ldsc which performs a LDscore regression. 
		message("Using LD score regression to estimate h2.")
		ldsc <- snp_ldsc2(corr0, df_beta)
		h2_est  <- ldsc[[2]]
		message("Estimated h2:",h2_est)
		if(h2_est < 1e-4){
			message("h2_est is too low, so I can't proceed, skipping chr ", chr,".")
			next
		}
		# We can also use the Automatic model to test many different h2 and p parameters and choose the best model
# takes a few minutes
		message("Computing LDpred2 auto, with initial p = ",initial_p, ", initial h2 = ", h2_est,", and ",iterations," iterations (+1000 burn-in)")
		multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est, vec_p_init = initial_p,  num_iter=iterations, ncores = NCORES)
		message("Done! p_est = ", multi_auto[[1]]$p_est,", and h2_est = ", multi_auto[[1]]$h2_est,".")
		info_snp.chr$beta_auto <- multi_auto[[1]]$beta_est 
		results <- rbind(results,info_snp.chr)
		unlink("../LDpred2/ref-data/*bk")
	}
		return(results)
}

# Load dataset



### BENCHMARKING!!!

N0 <- 107715	
N1 <- 19954 
data <- "../datasets/AST_Demenais_29273806_1-sAUCfilt.tsv.gz"

rapidoPGS <- function(dataset,N0,N1){
ds <- pgs.file.preprocess(dataset)
full.pgs  <- computePGS(ds,N0,N1)
return(full.pgs)
}

mbm <- microbenchmark(rapidoPGS(data,N0,N1),
	       computePGS.ldpred2.auto(data, N0,N1, iterations=500),
	       computePGS.ldpred2.auto(data, N0,N1, iterations=1000),
	       computePGS.ldpred2.auto(data, N0,N1, iterations=3000),
	       times=1)
mbm
