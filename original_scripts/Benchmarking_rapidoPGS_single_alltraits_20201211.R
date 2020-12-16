# Benchmarking RapidoPGS for all traits
# 2020/12/11

# This script will benchmark RápidoPGS-single (new function with no ALT_FREQ requirement) for all traits, using default (0.2) for case-control datasets, and custom priors for quantitative datasets.

# Load required libraries
#library(RapidoPGS)
library(bigsnpr)
library(microbenchmark)
library(susieR)
library(GenomicRanges)
library(data.table)

#############################################
########   Relevant functions   #############
#############################################

 ###' Helper function to sum logs without loss of precision
 ###'
 ###' Sums logs without loss of precision
 ###' This function is verbatim of its namesake in cupcake package (github.com/ollyburren/cupcake/)
 ###'
 ###' @param x a vector of logs to sum
 ###' @return a scalar
 ###' @author Chris Wallace
logsum <- function(x) {
    my.max <- max(x) ##take out the maximum value in log form)
    my.res <- my.max + log(sum(exp(x - my.max )))
    return(my.res)
}


# ##' Estimate trait standard deviation given vectors of variance of coefficients,  MAF and sample size
# ##'
# ##' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
# ##' var(X) = 2*maf*(1-maf)
# ##' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
# ##' This function is verbatim from its namesake in coloc package (github.com/chr1swallace/coloc/), by Chris Wallace
# ##' 
# ##' @title Estimate trait variance, internal function
# ##' @param vbeta vector of variance of coefficients
# ##' @param maf vector of MAF (same length as vbeta)
# ##' @param n sample size
# ##' @return estimated standard deviation of Y
# ##' @author Chris Wallace
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


# ##' Compute posterior probabilities using Wakefield's approximate Bayes Factors for quantitative traits
# ##'
# ##' \code{wakefield_pp_quant} computes posterior probabilities for a given SNP to be causal for a given SNP under the assumption of a single causal variant.
# ##'
# ##' This function was adapted from \code{wakefield_pp} in cupcake package (github.com/ollyburren/cupcake/)
# ##'
# ##' @param beta a vector of effect sizes (\eqn{\beta}) from a quantitative trait GWAS
# ##' @param se.beta vector of standard errors of effect sizes (\eqn{\beta})
# ##' @param sdY a scalar of the standard deviation given vectors of variance of coefficients,  MAF and sample size. Can be calculated using \code{sdY.est}
# ##' @param sd.prior a scalar representing our prior expectation of \eqn{\beta} (DEFAULT 0.15).
# ##' @param pi_i a scalar representing the prior probability (DEFAULT \eqn{1 \times 10^{-4}})
# ##' The method assumes a normal prior on the population log relative risk centred at 0 and the DEFAULT
# ##' value sets the variance of this distribution to 0.04, equivalent to a 95\%  belief that the true relative risk
# ##' is in the range of 0.66-1.5 at any causal variant.
# ##' @return a vector of posterior probabilities.
# ##' @author Guillermo Reales, Chris Wallace
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


# ##' compute posterior probabilities using Wakefield's approximate Bayes Factors
# ##' \code{wakefield_pp} computes posterior probabilities for a given SNP to be causal for a given SNP under the assumption of a single causal variant.
# ##'
# ##' This function was adapted from its namesake in cupcake package (github.com/ollyburren/cupcake/) to no longer require allele frequencies.
# ##'
# ##' @param z a vector of univariate Z scores from a GWAS
# ##' @param f a vector of minor allele frequencies taken from some reference population.
# ##' @param N a scalar or vector for total sample size of GWAS
# ##' @param s a scalar representing the proportion of cases (n.cases/N)
# ##' @param pi_i a scalar representing the prior probability (DEFAULT \eqn{1 \times 10^{-4}})
# ##' @param sd.prior a scalar representing our prior expectation of \eqn{\beta} (DEFAULT 0.2).
# ##' The method assumes a normal prior on the population log relative risk centred at 0 and the DEFAULT
# ##' value sets the variance of this distribution to 0.04, equivalent to a 95\%  belief that the true relative risk
# ##' is in the range of 0.66-1.5 at any causal variant.
# ##' @param log.p if FALSE (DEFAULT), p is a p value. If TRUE, p is a log(p) value.  Use this if your dataset holds p values too small to be accurately stored without using logs
# ##' @return a vector of posterior probabilities.
# ##' @author Olly Burren, Chris Wallace, Guillermo Reales
wakefield_pp <- function(beta, se, pi_i=1e-4,sd.prior=0.2) {
  if(length(beta) != length(se))
    stop("beta and se must be vectors of the same size")
    z=beta/se
    # compute V
    V  <- se^2
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


# ##' Compute PGS from GWAS summary statistics using posteriors from Wakefield's approximate Bayes Factors
# ##' 
# ##' '\code{rapidopgs_single} computes PGS from a from GWAS summary statistics using posteriors from Wakefield's approximate Bayes Factors
# ##' 
# ##' Main RápidoPGS function. This function will take a GWAS summary statistic dataset as an input,
# ##' will assign align it to a reference panel file (if provided), then it will assign 
# ##' SNPs to LD blocks and compute Wakefield's ppi by LD block, then will use it 
# ##' to generate PGS weights by multiplying those posteriors by effect sizes (\eqn{\beta}). 
# ##' Optionally, it will filter SNPs by a custom filter on ppi and then recalculate weights, to improve accuracy.
# ##' 
# ##' Alternatively, if filt_threshold is larger than one, RápidoPGS will select the top
# ##' \code{filt_threshold} SNPs by absolute weights (note, not ppi but weights).
# ##' 
# ##' The GWAS summary statistics file to compute PGS using our method must contain the following minimum columns, with these exact column names:
# ##' \describe{
# ##'   \item{CHR}{Chromosome}
# ##'   \item{BP}{Base position (in GRCh37/hg19 or GRCh38/hg38). If using hg38, use build = "hg38" in parameters}
# ##'   \item{SNPID}{rsids, or SNP identifiers. If not available, they can be anything (eg. CHR_BP)}
# ##'   \item{REF}{Reference, or non-effect allele}
# ##'   \item{ALT}{Alternative, or effect allele, the one \eqn{\beta} refers to}
# ##'   \item{ALT_FREQ}{Minor/ALT allele frequency in the tested population, or in a close population from a reference panel}
# ##'   \item{BETA}{\eqn{\beta} (or log(OR)), or effect sizes}
# ##'   \item{SE}{standard error of \eqn{\beta}}
# ##' }
# ##'
# ##' If a reference is provided. It should have 5 columns: CHR, BP,
# ##' SNPID, REF, and ALT. Also, it should be in the same build as 
# ##' the summary statistics. In both files, column order does not matter.
# ##' @param data a data.table containing GWAS summary statistic dataset
# ##'   with all required information.
# ##' @param N0 a scalar representing the number of controls in the
# ##'   study (or the number of subjects in quantitative trait GWAS),
# ##'   or a string indicating the column name containing it.
# ##' @param N1 a scalar representing the number of cases in the
# ##'   case-control study, or a string indicating the column name containing it. 
# ##'   If NULL (DEFAULT), quantitative trait will be assumed.
# ##' @param build a string containing the genome build of the dataset,
# ##'   either "hg19" (for hg19/GRCh37) or "hg38" (hg38/GRCh38). DEFAULT
# ##'   "hg19".
# ##' @param pi_i a scalar representing the prior probability (DEFAULT:
# ##'   \eqn{1 \times 10^{-4}}).
# ##' @param sd.prior the prior specifies that BETA at causal SNPs
# ##'   follows a centred normal distribution with standard deviation
# ##'   sd.prior. Sensible and widely used DEFAULTs are 0.2 for case
# ##'   control traits, and 0.15 * var(trait) for quantitative (selected
# ##'   if N1 is NULL).
# ##' @param filt_threshold a scalar indicating the ppi threshold (if
# ##'   \code{filt_threshold} < 1) or the number of top SNPs by absolute
# ##'   weights (if \code{filt_threshold} >= 1) to filter the dataset
# ##'   after PGS computation. If NULL (DEFAULT), no thresholding will
# ##'   be applied.
# ##' @param recalc a logical indicating if weights should be
# ##'   recalculated after thresholding. Only relevant if \code{filt_threshold}
# ##'   is defined.
# ##' @param reference a string indicating the path of the reference file 
# ##'   SNPs should be filtered and aligned to, see Details.
# ##' @param forsAUC a logical indicating if output should be in sAUC
# ##'   evaluation format as we used it for the paper.
# ##' @param altformat a logical indicating if output should be in a
# ##'   format containing pid (chr:pos), ALT, and weights only. DEFAULT
# ##'   FALSE
# ##' @return a data.table containing the formatted sumstats dataset with
# ##'   computed PGS weights.
# ##' @import data.table 
# ##' @importFrom bigsnpr snp_match
# ##' @importFrom GenomicRanges GRanges findOverlaps
# ##' @importFrom IRanges IRanges
# ##' @export
# ##' @author Guillermo Reales, Chris Wallace
# ##' @examples
# ##' sumstats <- data.table(SNPID=c("rs139096444","rs3843766","rs61977545", "rs544733737",
# ##'			"rs2177641", "rs183491817", "rs72995775","rs78598863", "rs1411315"), 
# ##'			CHR=c("4","20","14","2","4","6","6","21","13"), 
# ##'			BP=c(1479959, 13000913, 29107209, 203573414, 57331393, 11003529, 149256398, 
# ##'					25630085, 79166661), 
# ##'			REF=c("C","C","C","T","G","C","C","G","T"), 
# ##'			ALT=c("A","T","T","A","A","A","T","A","C"), 
# ##'			ALT_FREQ=c(0.2611,0.4482,0.0321,0.0538,0.574,0.0174,0.0084,0.0304,0.7528),
# ##'			BETA=c(0.012,0.0079,0.0224,0.0033,0.0153,0.058,0.0742,0.001,-0.0131),
# ##'			SE=c(0.0099,0.0066,0.0203,0.0171,0.0063,0.0255,0.043,0.0188,0.0074))
# ##'
# ##' PGS  <- rapidopgs_single(sumstats,  N0= 119078 ,N1=137045, build = "hg38")
# ##'
rapidopgs_single <- function(data,
                       N=NULL,
		       trait=c("cc","quant"),
                       build = "hg19",
                       pi_i= 1e-04,
                       sd.prior=if(trait == "quant") {0.15} else {0.2},
                       filt_threshold = NULL,
                       recalc=TRUE,
                       reference=NULL
                       ){

  if(!"data.table" %in% class(data))
	  data <- as.data.table(data)
  if(!trait %in% c("cc", "quant")) stop("Please, specify your study type, choose case-control ('cc') or quantitative ('quant').")
  if(length(trait) != 1) stop("Please select only one study type")
  if(trait == "quant" && is.null(N)) stop("N (sample size) is required for quantitative traits, please provide them, either as an integer or as column name containing it.")
  if(trait == "quant"){
  	mincol <- c("CHR","BP", "REF","ALT","BETA", "SE", "ALT_FREQ")
	  if(!all(mincol %in% names(data)))
    	stop("All minimum columns should be present in the summary statistics dataset. Please check, missing columns: ",
         paste(setdiff(mincol,names(data)), collapse=", "))  
  	ds <- copy(data) # avoid modifying input data.table
	ds  <- na.omit(ds, cols=c("CHR","BP", "BETA", "SE", "ALT_FREQ")) # Remove NA in relevant columns
  } else{
  	## Here's a list of columns that the dataset must have, otherwise it will fail
 	mincol <- c("CHR","BP", "REF","ALT","BETA", "SE")
  	if(!all(mincol %in% names(data)))
  	  stop("All minimum columns should be present in the summary statistics dataset. Please check, missing columns: ",
  	       paste(setdiff(mincol,names(data)), collapse=", "))
  	ds <- copy(data) # avoid modifying input data.table
  	ds  <- na.omit(ds, cols=c("CHR","BP", "BETA", "SE")) # Remove NA in relevant columns
  }	
  if(!"SNPID" %in% names(ds)){
	ds[,SNPID:=paste(CHR,BP, sep=":")] # SNPID is not strictly required to be provided. If it's not, we create it using CHR:BP
	       }
  ## First step, align to reference, if reference is provided
  if(!is.null(reference)){
          ## Next step is to filter and align our alleles and their effects to the hapmap3 reference, which I have already formatted for our purposes.
	message("Filtering SNPs...")
	refset <- fread(reference)
	mincolref <- c("CHR","BP", "SNPID", "REF","ALT")
	if(!all(mincolref %in% names(refset))) stop("All minimum columns should be present in the reference file. Please check, missing columns: ", paste(setdiff(mincol,names(refset)), collapse=", "))
	
	setnames(refset, old=c("CHR","BP", "SNPID", "REF","ALT"), new=c("chr","pos","id","a0","a1"))
	setnames(ds, old=c("CHR","BP", "BETA", "REF","ALT"), new=c("chr","pos","beta","a0","a1"))

	message("Matching and aligning SNPs to the reference")
	info_snp <- snp_match(ds,refset)
		
	message("Original sumstat had ", nrow(sumstats.chr)," for chr",chr,". After matching ", nrow(info_snp.chr)," remained, and will be used for further steps.")
	ds <- data.table(info_snp)
	ds[,SNPID:=id][, c("_NUM_ID_.ss","_NUM_ID_", "id"):=NULL] # Replaces SNPID in the sumstat by SNPID in the reference, and removes snp_match cols.
	setnames(ds, old=c("chr","pos","beta","a0","a1"), new=c("CHR","BP", "BETA", "REF","ALT"))
   }
	
   # Assign ld.blocks, in case they werent computed yet
   if(!"ld.block" %in% names(ds)){
   	if(build == "hg19"){ 
   		blranges <- RapidoPGS::EUR_ld.blocks
   	}else if(build == "hg38"){ 
   		blranges <- RapidoPGS::EUR_ld.blocks38
   	}else{ 
   		stop("RapidoPGS only accepts hg19 or hg38 at the moment, please check.")
   	}		
   message("Assigning LD blocks...")
   snpranges <- GRanges(seqnames=paste("chr",ds$CHR, sep=""), ranges=IRanges(start=ds$BP, end=ds$BP, names=ds$SNPID), strand="*")
   ds[,ld.block:=findOverlaps(snpranges, blranges, select='last')]
   message("Done!")
   }

  if(trait == "quant"){
    if(is.numeric(N) && length(N) == 1) { # In case N is supplied as a number
      message("Computing a RápidoPGS-single model for a quantitative trait with", N, " individuals...")
      ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=N)][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]
      if(!is.null(filt_threshold)){
        ds  <- ds[ds$ppi > filt_threshold,]
        if(recalc){
          ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=N)][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]		
        }
      }
    } else{ # In case column name is supplied
      if(is.character(N) && length(N) == 1){
        Nco <- N
        message("Computing a RápidoPGS-single model for a quantitative trait multiple with multiple N, supplied by column ", Nco,". Mmax N: ", max(ds[,get(Nco)]), ", min N = ", min(ds[,get(Nco)]), ", and mean N: ",mean(ds[,get(Nco)]), "...")
        ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=get(Nco))][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]
	if(!is.null(filt_threshold)){
          ds  <- ds[ds$ppi > filt_threshold,]
          if(recalc){
            ds[,sdY:=sdY.est(vbeta=SE^2, maf=ALT_FREQ, n=get(Nco))][,ppi:=wakefield_pp_quant(BETA,SE,sdY,sd.prior), by="ld.block"][,weight:=ppi*BETA]		
          }
	}
      }
    }
  }
 if(trait == "cc"){
      message("Computing a RápidoPGS-single model for a case-control dataset...")
      ds[,ppi:=wakefield_pp(beta=BETA, se=SE,pi_i, sd.prior), by = "ld.block"][, weight:=ppi*BETA]
      if(!is.null(filt_threshold)){
        if(filt_threshold < 1){
          ds  <- ds[ds$ppi > filt_threshold,]
        } else {
          ds <- ds[order(-rank(abs(weight))),][1:filt_threshold,] 
        }
        if(recalc){
          ds[,ppi:=wakefield_pp(beta=BETA, se=SE,pi_i, sd.prior), by = "ld.block"][, weight:=ppi*BETA]
        }
      }
  }
  return(ds)
}	




### BENCHMARKING!!!
set.seed(1)
dspath <- "../datasets/"
datasets <- c("AST_Demenais_29273806_1-qcfilt.tsv.gz","RA_Okada_24390342_3-qcfilt.tsv.gz", "T1D_Cooper_doi101101120022_1-qcfilt.tsv.gz", "T2D_Scott_28566273_1-qcfilt.tsv.gz", "BRCA_Michailidou_29059683_1-qcfilt.tsv.gz","PRCA_Schumacher_29892016_1-qcfilt.tsv.gz","CAD_Nikpay_26343387_1-qcfilt.tsv.gz","MDD_Wray_29700475_1-qcfilt.tsv.gz","BMI_Locke_25673413_1-qcfilt.tsv.gz","HEIGHT_Wood_25282103_1-qcfilt.tsv.gz")
data  <- paste0(dspath,datasets)


rapidoPGS_single <- function(data, N=NULL, trait="cc", sd.prior=0.2){
ds <- fread(data)
if(is.character(ds$CHR)){
ds[,CHR:=as.numeric(CHR)]
}
if(is.null(N)){
full.pgs  <- rapidopgs_single(ds,trait=trait, sd.prior=sd.prior)
} else{
full.pgs  <- rapidopgs_single(ds,N=N,trait=trait, sd.prior=sd.prior)
	       }
return(full.pgs)
}


results <- data.table()

# First the CC traits

for(i in 1:8){

mbm <- microbenchmark(rapidoPGS_single(data[i],trait="cc"), times=1, unit="s")
mbm <- as.data.table(summary(mbm))
mbm[,expr:=gsub("../datasets/","", data[i])]
results <- rbind(results,mbm)

}

# Then for the quantitative traits, which are a bit different
N <- "N"
priors <- c(0.08245557,0.09135172) # These are sd.priors, see "Compute_custom_priors_BMI_Height_clean_20201112.R" to find out how they were calculated
data <- data[9:10]

for(i in 1:2){

mbm <- microbenchmark(rapidoPGS_single(data[i],N=N, trait="quant", sd.prior=priors[i]), times=1, unit="s")
mbm <- as.data.table(summary(mbm))
mbm[,expr:=gsub("../datasets/","", data[i])]
results <- rbind(results,mbm)
}
fwrite(results, "../evaluations/Results_benchmark_RapidoPGS_single_all_20201211.tsv",sep="\t")

