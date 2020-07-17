# This  script is intended to develop PGS scores for the remaining traits, using rapidoPGS

# Here we'll develop PGS scores using the datasets used in previous rounds, but prefiltered by SNPs contained in sAUC .bim files. 

# NOTE: This version of rapidoPGS_model_creation has the only goal to improve inference on T1D dataset.
# This dataset had P = 0 for a number of SNPs, and this made RapidoPGS fail for the entire ld.bloc, so we lost important information.
# We introduced a small change in ComputePGS and wakefield_pp, to use log(P) instead of vanilla P.
# We then created new models for T1D only, with the hope of them being better at evaluations.

# Let's load the libraries
library(data.table)
library(GenomicRanges)


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

wakefield_pp <- function(p,f, N, s,pi_i=1e-4,sd.prior=0.2,log.p=FALSE) {
    if(length(p) != length(f))
      stop("p and f must be vectors of the same size")
    # compute V
    V <- 1 / (2 * N * f * (1 - f) * s * (1 - s))
    # convert p vals to z
    z <- stats::qnorm(0.5 * p, lower.tail = FALSE,log.p=log.p)
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

computePGS <- function(ds, N0,N1=NULL,pi_i= 1e-04, sd.prior=0.2, log.p=FALSE, filt_threshold = NULL, recalc=FALSE, forsAUC=FALSE, altformat=FALSE){

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
			if(log.p){
				ds[,P:=pnorm(-abs(BETA/SE),log.p=TRUE)*2]	
			}
			ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N=get(Nco)+get(Nca), s = get(Nca)/(get(Nco)+get(Nca)), pi_i, sd.prior, log.p), by = "ld.block"][, weight:=ppi*BETA]	
				if(!is.null(filt_threshold)){
					ds  <- ds[ds$ppi > filt_threshold,]
				if(recalc){
					ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N= get(Nco)+get(Nca) , s = get(Nca)/(get(Nca)+get(Nca)), pi_i, sd.prior, log.p), by = "ld.block"][, weight:=ppi*BETA]
				}
				}
		}else{
			message("Computing PGS for a case-control dataset, with ", N0," controls, and ", N1, " cases.")
			if(log.p){
				ds[,P:=pnorm(-abs(BETA/SE),log.p=TRUE)*2]	
			}
			ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N= N0+N1, s = N1/(N0+N1), pi_i, sd.prior, log.p), by = "ld.block"][, weight:=ppi*BETA]	
				if(!is.null(filt_threshold)){
					ds  <- ds[ds$ppi > filt_threshold,]
				if(recalc){
					ds[,ppi:=wakefield_pp(p = P, f = ALT_FREQ, N= N0+N1 , s = N1/(N0+N1), pi_i, sd.prior, log.p), by = "ld.block"][, weight:=ppi*BETA]
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


###############################################
#### GENERATING THE PGS MODELS   #############
#############################################


##### Asthma (AST) dataset, Martin filter
#
## This is a re-make of the asthma models. I could use the old ones, but I think it's better just to remake them from scracth, since I had to do the same for LDpred2
#
#trait <- "AST_Demenais_29273806_1_hapmap3qcfilt"
#N0 <- 107715  
#N1 <- 19954
#
## NOTE: We have to use the reference, since the original file was not aligned.
## This file has some extra columns that we don't need, so we'll remove them before proceeding.
#tmp <- fread("../datasets/AST_Demenais_29273806_1-hapmap3qcfilt.tsv.gz", select = c("SNPID", "CHR19","BP19", "REF","ALT","BETA", "SE", "P", "ALT_FREQ"))
#fwrite(tmp,"tmp.tsv.gz", sep="\t")
#
#ds <- pgs.file.preprocess("tmp.tsv.gz", ref="../references/reference_hapmap3_asthma_qcSNPs.txt")
## Also, let's keep the reference.SNPID
#ds[,SNPID:=SNPID.reference][,SNPID.reference:=NULL]
#
#full.pgs  <- computePGS(ds,N0,N1)
#full.pgs.model  <- computePGS(ds,N0,N1, forsAUC=TRUE)
#full.pgs.martin  <- computePGS(ds,N0,N1, altformat=TRUE)
#pgs.20k.model  <- full.pgs.model[order(-rank(abs(weight))),][1:20000,]
#pgs.20k.martin  <- full.pgs.martin[order(-rank(abs(weight))),][1:20000,]
#
#
#pgs.1e4.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, forsAUC=TRUE)
#pgs.1e3.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, forsAUC=TRUE)
#pgs.1e2.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, forsAUC=TRUE)
#pgs.1e1.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, forsAUC=TRUE)
#
#pgs.1e4recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, forsAUC=TRUE)
#pgs.1e3recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, forsAUC=TRUE)
#pgs.1e2recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, forsAUC=TRUE)
#pgs.1e1recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, forsAUC=TRUE)
#
#
#pgs.1e4.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, altformat=TRUE)
#pgs.1e3.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, altformat=TRUE)
#pgs.1e2.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, altformat=TRUE)
#pgs.1e1.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, altformat=TRUE)
#
#pgs.1e4recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, altformat=TRUE)
#pgs.1e3recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, altformat=TRUE)
#pgs.1e2recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, altformat=TRUE)
#pgs.1e1recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, altformat=TRUE)
#
#
#
#fwrite(pgs.20k.model, paste("../models/",trait, "-20k.prs.model",sep=""), sep="\t")
#fwrite(pgs.1e4.model, paste("../models/",trait, "-1e4.prs.model",sep=""), sep="\t")
#fwrite(pgs.1e4recalc.model, paste("../models/",trait, "-1e4recalc.prs.model",sep=""), sep="\t")
#fwrite(pgs.1e3.model, paste("../models/",trait, "-1e3.prs.model",sep=""), sep="\t")
#fwrite(pgs.1e3recalc.model, paste("../models/",trait, "-1e3recalc.prs.model",sep=""), sep="\t")
#fwrite(pgs.1e2.model, paste("../models/",trait, "-1e2.prs.model",sep=""), sep="\t")
#fwrite(pgs.1e2recalc.model, paste("../models/",trait, "-1e2recalc.prs.model",sep=""), sep="\t")
#fwrite(pgs.1e1.model, paste("../models/",trait, "-1e1.prs.model",sep=""), sep="\t")
#fwrite(pgs.1e1recalc.model, paste("../models/",trait, "-1e1recalc.prs.model",sep=""), sep="\t")
#
#
#fwrite(pgs.20k.martin, paste("../models/",trait, "-20k.prs.martin",sep=""), sep="\t")
#fwrite(pgs.1e4.martin, paste("../models/",trait, "-1e4.prs.martin",sep=""), sep="\t", col.names=F)
#fwrite(pgs.1e4recalc.martin, paste("../models/",trait, "-1e4recalc.prs.martin",sep=""), sep="\t", col.names=F)
#fwrite(pgs.1e3.martin, paste("../models/",trait, "-1e3.prs.martin",sep=""), sep="\t", col.names=F)
#fwrite(pgs.1e3recalc.martin, paste("../models/",trait, "-1e3recalc.prs.martin",sep=""), sep="\t", col.names=F)
#fwrite(pgs.1e2.martin, paste("../models/",trait, "-1e2.prs.martin",sep=""), sep="\t", col.names=F)
#fwrite(pgs.1e2recalc.martin, paste("../models/",trait, "-1e2recalc.prs.martin",sep=""), sep="\t", col.names=F)
#fwrite(pgs.1e1.martin, paste("../models/",trait, "-1e1.prs.martin",sep=""), sep="\t", col.names=F)
#fwrite(pgs.1e1recalc.martin, paste("../models/",trait, "-1e1recalc.prs.martin",sep=""), sep="\t", col.names=F)
#
#



### Asthma (AST) dataset, sAUC filter

trait <- "AST_Demenais_29273806_1_sAUCfilt"
N0 <- 107715  
N1 <- 19954

# NOTE: We have to use the reference, since the original file was not aligned.
# This file has some extra columns that we don't need, so we'll remove them before proceeding.
tmp <- fread("../datasets/AST_Demenais_29273806_1-sAUCfilt.tsv.gz", select = c("SNPID", "CHR19","BP19", "REF","ALT","BETA", "SE", "P", "ALT_FREQ"))
fwrite(tmp,"tmp.tsv.gz", sep="\t")

ds <- pgs.file.preprocess("tmp.tsv.gz", ref="../references/reference_sAUC.txt")
# Also, let's keep the reference.SNPID
ds[,SNPID:=SNPID.reference][,SNPID.reference:=NULL]
unlink("tmp.tsv.gz")

full.pgs  <- computePGS(ds,N0,N1)
full.pgs.model  <- computePGS(ds,N0,N1, forsAUC=TRUE)
full.pgs.martin  <- computePGS(ds,N0,N1, altformat=TRUE)
pgs.20k.model  <- full.pgs.model[order(-rank(abs(weight))),][1:20000,]
pgs.20k.martin  <- full.pgs.martin[order(-rank(abs(weight))),][1:20000,]


pgs.1e4.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, forsAUC=TRUE)
pgs.1e3.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, forsAUC=TRUE)
pgs.1e2.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, forsAUC=TRUE)
pgs.1e1.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, forsAUC=TRUE)

pgs.1e4recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, forsAUC=TRUE)
pgs.1e3recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, forsAUC=TRUE)
pgs.1e2recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, forsAUC=TRUE)
pgs.1e1recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, forsAUC=TRUE)


pgs.1e4.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, altformat=TRUE)
pgs.1e3.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, altformat=TRUE)
pgs.1e2.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, altformat=TRUE)
pgs.1e1.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, altformat=TRUE)

pgs.1e4recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, altformat=TRUE)
pgs.1e3recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, altformat=TRUE)
pgs.1e2recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, altformat=TRUE)
pgs.1e1recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, altformat=TRUE)



fwrite(pgs.20k.model, paste("../models/",trait, "-20k.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4.model, paste("../models/",trait, "-1e4.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4recalc.model, paste("../models/",trait, "-1e4recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3.model, paste("../models/",trait, "-1e3.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3recalc.model, paste("../models/",trait, "-1e3recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2.model, paste("../models/",trait, "-1e2.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2recalc.model, paste("../models/",trait, "-1e2recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1.model, paste("../models/",trait, "-1e1.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1recalc.model, paste("../models/",trait, "-1e1recalc.prs.model",sep=""), sep="\t")


fwrite(pgs.20k.martin, paste("../models/",trait, "-20k.prs.martin",sep=""), sep="\t")
fwrite(pgs.1e4.martin, paste("../models/",trait, "-1e4.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e4recalc.martin, paste("../models/",trait, "-1e4recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3.martin, paste("../models/",trait, "-1e3.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3recalc.martin, paste("../models/",trait, "-1e3recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2.martin, paste("../models/",trait, "-1e2.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2recalc.martin, paste("../models/",trait, "-1e2recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1.martin, paste("../models/",trait, "-1e1.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1recalc.martin, paste("../models/",trait, "-1e1recalc.prs.martin",sep=""), sep="\t", col.names=F)



#### Type 1 diabetes (T1D) dataset, sAUC filter


trait <- "T1D_Cooper_1_sAUCfilt"
N0 <- 8828    
N1 <- 5913 

# NOTE: We have to use the reference, since the original file was not aligned.
# This file has some extra columns that we don't need, so we'll remove them before proceeding.
tmp <- fread("../datasets/T1D_Cooper_1-sAUCfilt.tsv.gz", select = c("SNPID", "CHR19","BP19", "REF","ALT","BETA", "SE", "P", "ALT_FREQ"))
fwrite(tmp,"tmp.tsv.gz", sep="\t")

ds <- pgs.file.preprocess("tmp.tsv.gz", ref="../references/reference_sAUC.txt")
# Also, let's keep the reference.SNPID
ds[,SNPID:=SNPID.reference][,SNPID.reference:=NULL]
unlink("tmp.tsv.gz")
# Since it takes a while to preprocess this file, let us save it
fwrite(ds, "../datasets/T1D_Cooper_1-sAUCfilt_fullyprocessed.tsv.gz",sep="\t")

#NOTE: We have fully processed T1D file that we used in previous, so we don't need to re-align it or compute LD.blocks on it
# and save time at the very time-consuming alingment step
ds <- fread("../datasets/T1D_Cooper_1-sAUCfilt_fullyprocessed.tsv.gz")


full.pgs  <- computePGS(ds,N0,N1, log.p=TRUE)
full.pgs.model  <- computePGS(ds,N0,N1, forsAUC=TRUE, log.p=TRUE)
full.pgs.martin  <- computePGS(ds,N0,N1, altformat=TRUE, log.p=TRUE)
pgs.20k.model  <- full.pgs.model[order(-rank(abs(weight))),][1:20000,]
pgs.20k.martin  <- full.pgs.martin[order(-rank(abs(weight))),][1:20000,]
pgs.1e4.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, forsAUC=TRUE, log.p=TRUE)
pgs.1e3.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, forsAUC=TRUE, log.p=TRUE)
pgs.1e2.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, forsAUC=TRUE, log.p=TRUE)
pgs.1e1.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, forsAUC=TRUE, log.p=TRUE)

pgs.1e4recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, forsAUC=TRUE, log.p=TRUE)
pgs.1e3recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, forsAUC=TRUE, log.p=TRUE)
pgs.1e2recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, forsAUC=TRUE, log.p=TRUE)
pgs.1e1recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, forsAUC=TRUE, log.p=TRUE)


pgs.1e4.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, altformat=TRUE, log.p=TRUE)
pgs.1e3.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, altformat=TRUE, log.p=TRUE)
pgs.1e2.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, altformat=TRUE, log.p=TRUE)
pgs.1e1.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, altformat=TRUE, log.p=TRUE)

pgs.1e4recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, altformat=TRUE, log.p=TRUE)
pgs.1e3recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, altformat=TRUE, log.p=TRUE)
pgs.1e2recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, altformat=TRUE, log.p=TRUE)
pgs.1e1recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, altformat=TRUE, log.p=TRUE)



fwrite(pgs.20k.model, paste("../models/",trait, "-20k.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4.model, paste("../models/",trait, "-1e4.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4recalc.model, paste("../models/",trait, "-1e4recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3.model, paste("../models/",trait, "-1e3.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3recalc.model, paste("../models/",trait, "-1e3recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2.model, paste("../models/",trait, "-1e2.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2recalc.model, paste("../models/",trait, "-1e2recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1.model, paste("../models/",trait, "-1e1.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1recalc.model, paste("../models/",trait, "-1e1recalc.prs.model",sep=""), sep="\t")


fwrite(pgs.20k.martin, paste("../models/",trait, "-20k.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e4.martin, paste("../models/",trait, "-1e4.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e4recalc.martin, paste("../models/",trait, "-1e4recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3.martin, paste("../models/",trait, "-1e3.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3recalc.martin, paste("../models/",trait, "-1e3recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2.martin, paste("../models/",trait, "-1e2.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2recalc.martin, paste("../models/",trait, "-1e2recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1.martin, paste("../models/",trait, "-1e1.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1recalc.martin, paste("../models/",trait, "-1e1recalc.prs.martin",sep=""), sep="\t", col.names=F)


#### Type 1 diabetes (T1D) dataset, sAUC filter, noHLA

# This is a re-make of the asthma models. I could use the old ones, but I think it's better just to remake them from scracth, since I had to do the same for LDpred2

trait <- "T1DnoHLA_Cooper_1_sAUCfilt"
N0 <- 8828    
N1 <- 5913 

# NOTE: We have to use the reference, since the original file was not aligned.
# This file has some extra columns that we don't need, so we'll remove them before proceeding.
tmp <- fread("../datasets/T1DnoHLA_Cooper_1-sAUCfilt.tsv.gz", select = c("SNPID", "CHR19","BP19", "REF","ALT","BETA", "SE", "P", "ALT_FREQ"))
fwrite(tmp,"tmp.tsv.gz", sep="\t")

ds <- pgs.file.preprocess("tmp.tsv.gz", ref="../references/reference_sAUC.txt")
# Also, let's keep the reference.SNPID
ds[,SNPID:=SNPID.reference][,SNPID.reference:=NULL]
unlink("tmp.tsv.gz")

full.pgs  <- computePGS(ds,N0,N1)
full.pgs.model  <- computePGS(ds,N0,N1, forsAUC=TRUE)
full.pgs.martin  <- computePGS(ds,N0,N1, altformat=TRUE)
pgs.20k.model  <- full.pgs.model[order(-rank(abs(weight))),][1:20000,]
pgs.20k.martin  <- full.pgs.martin[order(-rank(abs(weight))),][1:20000,]


pgs.1e4.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, forsAUC=TRUE)
pgs.1e3.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, forsAUC=TRUE)
pgs.1e2.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, forsAUC=TRUE)
pgs.1e1.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, forsAUC=TRUE)

pgs.1e4recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, forsAUC=TRUE)
pgs.1e3recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, forsAUC=TRUE)
pgs.1e2recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, forsAUC=TRUE)
pgs.1e1recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, forsAUC=TRUE)


pgs.1e4.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, altformat=TRUE)
pgs.1e3.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, altformat=TRUE)
pgs.1e2.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, altformat=TRUE)
pgs.1e1.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, altformat=TRUE)

pgs.1e4recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, altformat=TRUE)
pgs.1e3recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, altformat=TRUE)
pgs.1e2recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, altformat=TRUE)
pgs.1e1recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, altformat=TRUE)



fwrite(pgs.20k.model, paste("../models/",trait, "-20k.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4.model, paste("../models/",trait, "-1e4.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4recalc.model, paste("../models/",trait, "-1e4recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3.model, paste("../models/",trait, "-1e3.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3recalc.model, paste("../models/",trait, "-1e3recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2.model, paste("../models/",trait, "-1e2.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2recalc.model, paste("../models/",trait, "-1e2recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1.model, paste("../models/",trait, "-1e1.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1recalc.model, paste("../models/",trait, "-1e1recalc.prs.model",sep=""), sep="\t")


fwrite(pgs.20k.martin, paste("../models/",trait, "-20k.prs.martin",sep=""), sep="\t")
fwrite(pgs.1e4.martin, paste("../models/",trait, "-1e4.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e4recalc.martin, paste("../models/",trait, "-1e4recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3.martin, paste("../models/",trait, "-1e3.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3recalc.martin, paste("../models/",trait, "-1e3recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2.martin, paste("../models/",trait, "-1e2.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2recalc.martin, paste("../models/",trait, "-1e2recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1.martin, paste("../models/",trait, "-1e1.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1recalc.martin, paste("../models/",trait, "-1e1recalc.prs.martin",sep=""), sep="\t", col.names=F)




#### Type 2 diabetes (T2D) dataset, sAUC filter

# This is a re-make of the asthma models. I could use the old ones, but I think it's better just to remake them from scracth, since I had to do the same for LDpred2

trait <- "T2D_Scott_1_sAUCfilt"
N0 <- 132532	
N1  <- 26676 

# NOTE: We have to use the reference, since the original file was not aligned.
# This file has some extra columns that we don't need, so we'll remove them before proceeding.
tmp <- fread("../datasets/T2D_Scott_28566273_1-sAUCfilt.tsv.gz", select = c("SNPID", "CHR19","BP19", "REF","ALT","BETA", "SE", "P", "ALT_FREQ"))
fwrite(tmp,"tmp.tsv.gz", sep="\t")

ds <- pgs.file.preprocess("tmp.tsv.gz", ref="../references/reference_sAUC.txt")
# Also, let's keep the reference.SNPID
ds[,SNPID:=SNPID.reference][,SNPID.reference:=NULL]
unlink("tmp.tsv.gz")

full.pgs  <- computePGS(ds,N0,N1)
full.pgs.model  <- computePGS(ds,N0,N1, forsAUC=TRUE)
full.pgs.martin  <- computePGS(ds,N0,N1, altformat=TRUE)
pgs.20k.model  <- full.pgs.model[order(-rank(abs(weight))),][1:20000,]
pgs.20k.martin  <- full.pgs.martin[order(-rank(abs(weight))),][1:20000,]


pgs.1e4.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, forsAUC=TRUE)
pgs.1e3.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, forsAUC=TRUE)
pgs.1e2.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, forsAUC=TRUE)
pgs.1e1.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, forsAUC=TRUE)

pgs.1e4recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, forsAUC=TRUE)
pgs.1e3recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, forsAUC=TRUE)
pgs.1e2recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, forsAUC=TRUE)
pgs.1e1recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, forsAUC=TRUE)


pgs.1e4.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, altformat=TRUE)
pgs.1e3.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, altformat=TRUE)
pgs.1e2.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, altformat=TRUE)
pgs.1e1.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, altformat=TRUE)

pgs.1e4recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, altformat=TRUE)
pgs.1e3recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, altformat=TRUE)
pgs.1e2recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, altformat=TRUE)
pgs.1e1recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, altformat=TRUE)



fwrite(pgs.20k.model, paste("../models/",trait, "-20k.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4.model, paste("../models/",trait, "-1e4.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4recalc.model, paste("../models/",trait, "-1e4recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3.model, paste("../models/",trait, "-1e3.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3recalc.model, paste("../models/",trait, "-1e3recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2.model, paste("../models/",trait, "-1e2.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2recalc.model, paste("../models/",trait, "-1e2recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1.model, paste("../models/",trait, "-1e1.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1recalc.model, paste("../models/",trait, "-1e1recalc.prs.model",sep=""), sep="\t")


fwrite(pgs.20k.martin, paste("../models/",trait, "-20k.prs.martin",sep=""), sep="\t")
fwrite(pgs.1e4.martin, paste("../models/",trait, "-1e4.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e4recalc.martin, paste("../models/",trait, "-1e4recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3.martin, paste("../models/",trait, "-1e3.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3recalc.martin, paste("../models/",trait, "-1e3recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2.martin, paste("../models/",trait, "-1e2.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2recalc.martin, paste("../models/",trait, "-1e2recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1.martin, paste("../models/",trait, "-1e1.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1recalc.martin, paste("../models/",trait, "-1e1recalc.prs.martin",sep=""), sep="\t", col.names=F)






#### Rheumathoid Arthritis (RA) Dataset 3

# I will generate a new rapidoPGS for this trait, as I had chosen the wrong dataset before
# NOTE: Case and control numbers in LDpred2 paper don't match those in the dataset website and/or paper. Using here the numbers found in the paper and at http://plaza.umin.ac.jp/~yokada/datasource/software.htm

trait <- "RA_Okada_24390342_3_sAUCfilt" 
N0 <- 61565	
N1 <- 19234	
 

ds <- pgs.file.preprocess("../datasets/RA_Okada_24390342_3-sAUCfilt.tsv.gz")

full.pgs  <- computePGS(ds,N0,N1)
full.pgs.model  <- computePGS(ds,N0,N1, forsAUC=TRUE)
full.pgs.martin  <- computePGS(ds,N0,N1, altformat=TRUE)
pgs.20k.model  <- full.pgs.model[order(-rank(abs(weight))),][1:20000,]
pgs.20k.martin  <- full.pgs.martin[order(-rank(abs(weight))),][1:20000,]


pgs.1e4.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, forsAUC=TRUE)
pgs.1e3.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, forsAUC=TRUE)
pgs.1e2.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, forsAUC=TRUE)
pgs.1e1.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, forsAUC=TRUE)

pgs.1e4recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, forsAUC=TRUE)
pgs.1e3recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, forsAUC=TRUE)
pgs.1e2recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, forsAUC=TRUE)
pgs.1e1recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, forsAUC=TRUE)


pgs.1e4.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, altformat=TRUE)
pgs.1e3.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, altformat=TRUE)
pgs.1e2.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, altformat=TRUE)
pgs.1e1.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, altformat=TRUE)

pgs.1e4recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, altformat=TRUE)
pgs.1e3recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, altformat=TRUE)
pgs.1e2recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, altformat=TRUE)
pgs.1e1recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, altformat=TRUE)



fwrite(pgs.20k.model, paste("../models/",trait, "-20k.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4.model, paste("../models/",trait, "-1e4.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4recalc.model, paste("../models/",trait, "-1e4recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3.model, paste("../models/",trait, "-1e3.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3recalc.model, paste("../models/",trait, "-1e3recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2.model, paste("../models/",trait, "-1e2.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2recalc.model, paste("../models/",trait, "-1e2recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1.model, paste("../models/",trait, "-1e1.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1recalc.model, paste("../models/",trait, "-1e1recalc.prs.model",sep=""), sep="\t")


fwrite(pgs.20k.martin, paste("../models/",trait, "-20k.prs.martin",sep=""), sep="\t")
fwrite(pgs.1e4.martin, paste("../models/",trait, "-1e4.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e4recalc.martin, paste("../models/",trait, "-1e4recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3.martin, paste("../models/",trait, "-1e3.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3recalc.martin, paste("../models/",trait, "-1e3recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2.martin, paste("../models/",trait, "-1e2.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2recalc.martin, paste("../models/",trait, "-1e2recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1.martin, paste("../models/",trait, "-1e1.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1recalc.martin, paste("../models/",trait, "-1e1recalc.prs.martin",sep=""), sep="\t", col.names=F)


#### Breast cancer (BRCA) Dataset


trait <- "BRCA_Michailidou_29059683_1_sAUCfilt" 
N0 <- 119078	
N1 <- 137045 

ds <- pgs.file.preprocess("../datasets/BRCA_Michailidou_29059683_1-sAUCfilt.tsv.gz")

full.pgs  <- computePGS(ds,N0,N1)
full.pgs.model  <- computePGS(ds,N0,N1, forsAUC=TRUE)
full.pgs.martin  <- computePGS(ds,N0,N1, altformat=TRUE)
pgs.20k.model  <- full.pgs.model[order(-rank(abs(weight))),][1:20000,]
pgs.20k.martin  <- full.pgs.martin[order(-rank(abs(weight))),][1:20000,]


pgs.1e4.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, forsAUC=TRUE)
pgs.1e3.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, forsAUC=TRUE)
pgs.1e2.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, forsAUC=TRUE)
pgs.1e1.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, forsAUC=TRUE)

pgs.1e4recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, forsAUC=TRUE)
pgs.1e3recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, forsAUC=TRUE)
pgs.1e2recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, forsAUC=TRUE)
pgs.1e1recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, forsAUC=TRUE)


pgs.1e4.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, altformat=TRUE)
pgs.1e3.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, altformat=TRUE)
pgs.1e2.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, altformat=TRUE)
pgs.1e1.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, altformat=TRUE)

pgs.1e4recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, altformat=TRUE)
pgs.1e3recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, altformat=TRUE)
pgs.1e2recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, altformat=TRUE)
pgs.1e1recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, altformat=TRUE)



fwrite(pgs.20k.model, paste("../models/",trait, "-20k.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4.model, paste("../models/",trait, "-1e4.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4recalc.model, paste("../models/",trait, "-1e4recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3.model, paste("../models/",trait, "-1e3.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3recalc.model, paste("../models/",trait, "-1e3recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2.model, paste("../models/",trait, "-1e2.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2recalc.model, paste("../models/",trait, "-1e2recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1.model, paste("../models/",trait, "-1e1.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1recalc.model, paste("../models/",trait, "-1e1recalc.prs.model",sep=""), sep="\t")


fwrite(pgs.20k.martin, paste("../models/",trait, "-20k.prs.martin",sep=""), sep="\t")
fwrite(pgs.1e4.martin, paste("../models/",trait, "-1e4.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e4recalc.martin, paste("../models/",trait, "-1e4recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3.martin, paste("../models/",trait, "-1e3.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3recalc.martin, paste("../models/",trait, "-1e3recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2.martin, paste("../models/",trait, "-1e2.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2recalc.martin, paste("../models/",trait, "-1e2recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1.martin, paste("../models/",trait, "-1e1.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1recalc.martin, paste("../models/",trait, "-1e1recalc.prs.martin",sep=""), sep="\t", col.names=F)


#### Prostate cancer (PRCA) Dataset


trait <- "PRCA_Schumacher_29892016_1_sAUCfilt" 
N0 <- 61106	
N1 <- 79148

ds <- pgs.file.preprocess("../datasets/PRCA_Schumacher_29892016_1-sAUCfilt.tsv.gz")

full.pgs  <- computePGS(ds,N0,N1)
full.pgs.model  <- computePGS(ds,N0,N1, forsAUC=TRUE)
full.pgs.martin  <- computePGS(ds,N0,N1, altformat=TRUE)
pgs.20k.model  <- full.pgs.model[order(-rank(abs(weight))),][1:20000,]
pgs.20k.martin  <- full.pgs.martin[order(-rank(abs(weight))),][1:20000,]


pgs.1e4.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, forsAUC=TRUE)
pgs.1e3.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, forsAUC=TRUE)
pgs.1e2.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, forsAUC=TRUE)
pgs.1e1.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, forsAUC=TRUE)

pgs.1e4recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, forsAUC=TRUE)
pgs.1e3recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, forsAUC=TRUE)
pgs.1e2recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, forsAUC=TRUE)
pgs.1e1recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, forsAUC=TRUE)


pgs.1e4.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, altformat=TRUE)
pgs.1e3.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, altformat=TRUE)
pgs.1e2.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, altformat=TRUE)
pgs.1e1.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, altformat=TRUE)

pgs.1e4recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, altformat=TRUE)
pgs.1e3recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, altformat=TRUE)
pgs.1e2recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, altformat=TRUE)
pgs.1e1recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, altformat=TRUE)



fwrite(pgs.20k.model, paste("../models/",trait, "-20k.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4.model, paste("../models/",trait, "-1e4.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4recalc.model, paste("../models/",trait, "-1e4recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3.model, paste("../models/",trait, "-1e3.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3recalc.model, paste("../models/",trait, "-1e3recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2.model, paste("../models/",trait, "-1e2.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2recalc.model, paste("../models/",trait, "-1e2recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1.model, paste("../models/",trait, "-1e1.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1recalc.model, paste("../models/",trait, "-1e1recalc.prs.model",sep=""), sep="\t")


fwrite(pgs.20k.martin, paste("../models/",trait, "-20k.prs.martin",sep=""), sep="\t")
fwrite(pgs.1e4.martin, paste("../models/",trait, "-1e4.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e4recalc.martin, paste("../models/",trait, "-1e4recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3.martin, paste("../models/",trait, "-1e3.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3recalc.martin, paste("../models/",trait, "-1e3recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2.martin, paste("../models/",trait, "-1e2.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2recalc.martin, paste("../models/",trait, "-1e2recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1.martin, paste("../models/",trait, "-1e1.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1recalc.martin, paste("../models/",trait, "-1e1recalc.prs.martin",sep=""), sep="\t", col.names=F)


#### Coronary artery disease (CAD)Dataset


trait <- "CAD_Nikpay_26343387_1_sAUCfilt"
N0 <- 123504	
N1 <- 60801
ds <- pgs.file.preprocess("../datasets/CAD_Nikpay_26343387_1-sAUCfilt.tsv.gz")

full.pgs  <- computePGS(ds,N0,N1)
full.pgs.model  <- computePGS(ds,N0,N1, forsAUC=TRUE)
full.pgs.martin  <- computePGS(ds,N0,N1, altformat=TRUE)
pgs.20k.model  <- full.pgs.model[order(-rank(abs(weight))),][1:20000,]
pgs.20k.martin  <- full.pgs.martin[order(-rank(abs(weight))),][1:20000,]


pgs.1e4.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, forsAUC=TRUE)
pgs.1e3.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, forsAUC=TRUE)
pgs.1e2.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, forsAUC=TRUE)
pgs.1e1.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, forsAUC=TRUE)

pgs.1e4recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, forsAUC=TRUE)
pgs.1e3recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, forsAUC=TRUE)
pgs.1e2recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, forsAUC=TRUE)
pgs.1e1recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, forsAUC=TRUE)


pgs.1e4.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, altformat=TRUE)
pgs.1e3.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, altformat=TRUE)
pgs.1e2.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, altformat=TRUE)
pgs.1e1.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, altformat=TRUE)

pgs.1e4recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, altformat=TRUE)
pgs.1e3recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, altformat=TRUE)
pgs.1e2recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, altformat=TRUE)
pgs.1e1recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, altformat=TRUE)



fwrite(pgs.20k.model, paste("../models/",trait, "-20k.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4.model, paste("../models/",trait, "-1e4.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4recalc.model, paste("../models/",trait, "-1e4recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3.model, paste("../models/",trait, "-1e3.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3recalc.model, paste("../models/",trait, "-1e3recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2.model, paste("../models/",trait, "-1e2.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2recalc.model, paste("../models/",trait, "-1e2recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1.model, paste("../models/",trait, "-1e1.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1recalc.model, paste("../models/",trait, "-1e1recalc.prs.model",sep=""), sep="\t")


fwrite(pgs.20k.martin, paste("../models/",trait, "-20k.prs.martin",sep=""), sep="\t")
fwrite(pgs.1e4.martin, paste("../models/",trait, "-1e4.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e4recalc.martin, paste("../models/",trait, "-1e4recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3.martin, paste("../models/",trait, "-1e3.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3recalc.martin, paste("../models/",trait, "-1e3recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2.martin, paste("../models/",trait, "-1e2.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2recalc.martin, paste("../models/",trait, "-1e2recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1.martin, paste("../models/",trait, "-1e1.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1recalc.martin, paste("../models/",trait, "-1e1recalc.prs.martin",sep=""), sep="\t", col.names=F)


#### Major depression disorder (MDD) Dataset


trait <- "MDD_Wray_29700475_1_sAUCfilt"
N0 <- "N0"	
N1 <- "N1"
ds <- pgs.file.preprocess("../datasets/MDD_Wray_29700475_1-sAUCfilt.tsv.gz")

full.pgs  <- computePGS(ds,N0,N1)
full.pgs.model  <- computePGS(ds,N0,N1, forsAUC=TRUE)
full.pgs.martin  <- computePGS(ds,N0,N1, altformat=TRUE)
pgs.20k.model  <- full.pgs.model[order(-rank(abs(weight))),][1:20000,]
pgs.20k.martin  <- full.pgs.martin[order(-rank(abs(weight))),][1:20000,]


pgs.1e4.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, forsAUC=TRUE)
pgs.1e3.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, forsAUC=TRUE)
pgs.1e2.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, forsAUC=TRUE)
pgs.1e1.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, forsAUC=TRUE)

pgs.1e4recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, forsAUC=TRUE)
pgs.1e3recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, forsAUC=TRUE)
pgs.1e2recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, forsAUC=TRUE)
pgs.1e1recalc.model <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, forsAUC=TRUE)


pgs.1e4.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, altformat=TRUE)
pgs.1e3.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, altformat=TRUE)
pgs.1e2.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, altformat=TRUE)
pgs.1e1.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, altformat=TRUE)

pgs.1e4recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-04, recalc=TRUE, altformat=TRUE)
pgs.1e3recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-03, recalc=TRUE, altformat=TRUE)
pgs.1e2recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-02, recalc=TRUE, altformat=TRUE)
pgs.1e1recalc.martin <- computePGS(ds,N0,N1, filt_threshold = 1e-01, recalc=TRUE, altformat=TRUE)



fwrite(pgs.20k.model, paste("../models/",trait, "-20k.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4.model, paste("../models/",trait, "-1e4.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4recalc.model, paste("../models/",trait, "-1e4recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3.model, paste("../models/",trait, "-1e3.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3recalc.model, paste("../models/",trait, "-1e3recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2.model, paste("../models/",trait, "-1e2.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2recalc.model, paste("../models/",trait, "-1e2recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1.model, paste("../models/",trait, "-1e1.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1recalc.model, paste("../models/",trait, "-1e1recalc.prs.model",sep=""), sep="\t")


fwrite(pgs.20k.martin, paste("../models/",trait, "-20k.prs.martin",sep=""), sep="\t")
fwrite(pgs.1e4.martin, paste("../models/",trait, "-1e4.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e4recalc.martin, paste("../models/",trait, "-1e4recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3.martin, paste("../models/",trait, "-1e3.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3recalc.martin, paste("../models/",trait, "-1e3recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2.martin, paste("../models/",trait, "-1e2.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2recalc.martin, paste("../models/",trait, "-1e2recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1.martin, paste("../models/",trait, "-1e1.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1recalc.martin, paste("../models/",trait, "-1e1recalc.prs.martin",sep=""), sep="\t", col.names=F)


#### BMI dataset, custom qc filter


trait <- "BMI_Locke_25673413_1_qcfilt"


ds <- pgs.file.preprocess("../datasets/BMI_Locke_25673413_1-qcfilt.tsv.gz", ref="../references/reference_hapmap3_BMI_qcSNPs.txt")
# Also, let's keep the reference.SNPID
ds[,SNPID:=SNPID.reference][,SNPID.reference:=NULL]

full.pgs  <- computePGS(ds,N0="N",sd.prior=0.15)
full.pgs.model  <- computePGS(ds,N0="N",sd.prior=0.15, forsAUC=TRUE)
full.pgs.martin  <- computePGS(ds,N0="N",sd.prior=0.15, altformat=TRUE)
pgs.20k.model  <- full.pgs.model[order(-rank(abs(weight))),][1:20000,]
pgs.20k.martin  <- full.pgs.martin[order(-rank(abs(weight))),][1:20000,]


pgs.1e4.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-04, forsAUC=TRUE)
pgs.1e3.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-03, forsAUC=TRUE)
pgs.1e2.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-02, forsAUC=TRUE)
pgs.1e1.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-01, forsAUC=TRUE)

pgs.1e4recalc.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-04, recalc=TRUE, forsAUC=TRUE)
pgs.1e3recalc.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-03, recalc=TRUE, forsAUC=TRUE)
pgs.1e2recalc.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-02, recalc=TRUE, forsAUC=TRUE)
pgs.1e1recalc.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-01, recalc=TRUE, forsAUC=TRUE)


pgs.1e4.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-04, altformat=TRUE)
pgs.1e3.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-03, altformat=TRUE)
pgs.1e2.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-02, altformat=TRUE)
pgs.1e1.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-01, altformat=TRUE)

pgs.1e4recalc.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-04, recalc=TRUE, altformat=TRUE)
pgs.1e3recalc.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-03, recalc=TRUE, altformat=TRUE)
pgs.1e2recalc.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-02, recalc=TRUE, altformat=TRUE)
pgs.1e1recalc.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-01, recalc=TRUE, altformat=TRUE)



fwrite(pgs.20k.model, paste("../models/",trait, "-20k.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4.model, paste("../models/",trait, "-1e4.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4recalc.model, paste("../models/",trait, "-1e4recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3.model, paste("../models/",trait, "-1e3.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3recalc.model, paste("../models/",trait, "-1e3recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2.model, paste("../models/",trait, "-1e2.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2recalc.model, paste("../models/",trait, "-1e2recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1.model, paste("../models/",trait, "-1e1.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1recalc.model, paste("../models/",trait, "-1e1recalc.prs.model",sep=""), sep="\t")


fwrite(pgs.20k.martin, paste("../models/",trait, "-20k.prs.martin",sep=""), sep="\t")
fwrite(pgs.1e4.martin, paste("../models/",trait, "-1e4.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e4recalc.martin, paste("../models/",trait, "-1e4recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3.martin, paste("../models/",trait, "-1e3.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3recalc.martin, paste("../models/",trait, "-1e3recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2.martin, paste("../models/",trait, "-1e2.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2recalc.martin, paste("../models/",trait, "-1e2recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1.martin, paste("../models/",trait, "-1e1.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1recalc.martin, paste("../models/",trait, "-1e1recalc.prs.martin",sep=""), sep="\t", col.names=F)


### HEIGHT dataset, custom qc filter


trait <- "HEIGHT_Wood_25282103_1_qcfilt"

ds <- pgs.file.preprocess("../datasets/HEIGHT_Wood_25282103_1-qcfilt.tsv.gz", ref="../references/reference_hapmap3_Height_qcSNPs.txt")
# Also, let's keep the reference.SNPID
ds[,SNPID:=SNPID.reference][,SNPID.reference:=NULL]

full.pgs  <- computePGS(ds,N0="N",sd.prior=0.15)
full.pgs.model  <- computePGS(ds,N0="N",sd.prior=0.15, forsAUC=TRUE)
full.pgs.martin  <- computePGS(ds,N0="N",sd.prior=0.15, altformat=TRUE)
pgs.20k.model  <- full.pgs.model[order(-rank(abs(weight))),][1:20000,]
pgs.20k.martin  <- full.pgs.martin[order(-rank(abs(weight))),][1:20000,]


pgs.1e4.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-04, forsAUC=TRUE)
pgs.1e3.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-03, forsAUC=TRUE)
pgs.1e2.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-02, forsAUC=TRUE)
pgs.1e1.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-01, forsAUC=TRUE)

pgs.1e4recalc.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-04, recalc=TRUE, forsAUC=TRUE)
pgs.1e3recalc.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-03, recalc=TRUE, forsAUC=TRUE)
pgs.1e2recalc.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-02, recalc=TRUE, forsAUC=TRUE)
pgs.1e1recalc.model <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-01, recalc=TRUE, forsAUC=TRUE)


pgs.1e4.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-04, altformat=TRUE)
pgs.1e3.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-03, altformat=TRUE)
pgs.1e2.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-02, altformat=TRUE)
pgs.1e1.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-01, altformat=TRUE)

pgs.1e4recalc.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-04, recalc=TRUE, altformat=TRUE)
pgs.1e3recalc.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-03, recalc=TRUE, altformat=TRUE)
pgs.1e2recalc.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-02, recalc=TRUE, altformat=TRUE)
pgs.1e1recalc.martin <- computePGS(ds,N0="N",sd.prior=0.15, filt_threshold = 1e-01, recalc=TRUE, altformat=TRUE)



fwrite(pgs.20k.model, paste("../models/",trait, "-20k.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4.model, paste("../models/",trait, "-1e4.prs.model",sep=""), sep="\t")
fwrite(pgs.1e4recalc.model, paste("../models/",trait, "-1e4recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3.model, paste("../models/",trait, "-1e3.prs.model",sep=""), sep="\t")
fwrite(pgs.1e3recalc.model, paste("../models/",trait, "-1e3recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2.model, paste("../models/",trait, "-1e2.prs.model",sep=""), sep="\t")
fwrite(pgs.1e2recalc.model, paste("../models/",trait, "-1e2recalc.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1.model, paste("../models/",trait, "-1e1.prs.model",sep=""), sep="\t")
fwrite(pgs.1e1recalc.model, paste("../models/",trait, "-1e1recalc.prs.model",sep=""), sep="\t")


fwrite(pgs.20k.martin, paste("../models/",trait, "-20k.prs.martin",sep=""), sep="\t")
fwrite(pgs.1e4.martin, paste("../models/",trait, "-1e4.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e4recalc.martin, paste("../models/",trait, "-1e4recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3.martin, paste("../models/",trait, "-1e3.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e3recalc.martin, paste("../models/",trait, "-1e3recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2.martin, paste("../models/",trait, "-1e2.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e2recalc.martin, paste("../models/",trait, "-1e2recalc.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1.martin, paste("../models/",trait, "-1e1.prs.martin",sep=""), sep="\t", col.names=F)
fwrite(pgs.1e1recalc.martin, paste("../models/",trait, "-1e1recalc.prs.martin",sep=""), sep="\t", col.names=F)

