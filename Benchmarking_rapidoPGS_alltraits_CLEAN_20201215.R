# Benchmarking RapidoPGS with set priors for all traits (Figure S1).
# 2020/12/15

# Here we'll compare the timing for PGS generation using different rapidoPGS methods on all traits
# This script will use a set prior for all traits in rapidopgs_multi, including 0.2 for binary traits

# Load required libraries
library(RapidoPGS)
library(bigsnpr)
library(microbenchmark)
library(susieR)
library(GenomicRanges)


#############################################
########   Relevant functions   #############
#############################################

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

# NOTE: This is the version used in benchmarking. Not the last version, as we later added a bug fix (to filter SE = 0 SNPs) that shouldn't affect run times.
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



######################
###   BENCHMARK    ###
######################
set.seed(1)


dspath <- "../datasets/"
datasets <- c("AST_Demenais_29273806_1-qcfilt.tsv.gz","RA_Okada_24390342_3-qcfilt.tsv.gz", "T1D_Cooper_doi101101120022_1-qcfilt.tsv.gz", "T2D_Scott_28566273_1-qcfilt.tsv.gz", "BRCA_Michailidou_29059683_1-qcfilt.tsv.gz","PRCA_Schumacher_29892016_1-qcfilt.tsv.gz","CAD_Nikpay_26343387_1-qcfilt.tsv.gz","MDD_Wray_29700475_1-qcfilt.tsv.gz","BMI_Locke_25673413_1-qcfilt.tsv.gz","HEIGHT_Wood_25282103_1-qcfilt.tsv.gz")
data  <- paste0(dspath,datasets)
ref  <- "../references/RapidoPGS-ref/"
ncores <- bigparallelr::nb_cores()
tests <- c("RapidoPGS_single", "RapidoPGS_multi,1e-3,0.1,prior0.2","RapidoPGS_multi,1e-4,0.01,prior0.2" )


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

rapidoPGS_mult <- function(data, alpha.block, alpha.snp, sd.prior=NULL){
ds <- fread(data)
if(is.character(ds$CHR)){
ds[,CHR:=as.numeric(CHR)]
}
full.pgs  <- rapidopgs_multi(ds,reference=ref, ncores=ncores, iterations=1000, alpha.block= alpha.block, alpha.snp=alpha.snp, sd.prior=sd.prior)
return(full.pgs)
} 

results <- data.table()

#####################################################
#####          RÁPIDOPGS-SINGLE               #######
#####################################################


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

#####################################################
#####          RÁPIDOPGS-MULTI (SET PRIOR)    #######
#####################################################

results <- data.table()
datasets <- c("AST_Demenais_29273806_1-qcfilt.tsv.gz","RA_Okada_24390342_3-qcfilt.tsv.gz", "T1D_Cooper_doi101101120022_1-qcfilt.tsv.gz", "T2D_Scott_28566273_1-qcfilt.tsv.gz", "BRCA_Michailidou_29059683_1-qcfilt.tsv.gz","PRCA_Schumacher_29892016_1-qcfilt.tsv.gz","CAD_Nikpay_26343387_1-qcfilt.tsv.gz","MDD_Wray_29700475_1-qcfilt.tsv.gz","BMI_Locke_25673413_1-qcfilt.tsv.gz","HEIGHT_Wood_25282103_1-qcfilt.tsv.gz")
data  <- paste0(dspath,datasets)

for(i in 1:8){

mbm <- microbenchmark(rapidoPGS_mult(data[i], alpha.block= 1e-3, alpha.snp=0.1, sd.prior=0.2),
	       rapidoPGS_mult(data[i], alpha.block= 1e-4, alpha.snp=0.01, sd.prior=0.2),
	       times=1, unit="s")
mbm <- as.data.table(summary(mbm))
mbm[,expr:=gsub("../datasets/","", data[i])][,method:=tests][,run:="set_prior"]
results <- rbind(results,mbm)
fwrite(results, "../evaluations/Results_benchmark_RapidoPGS_all_setprior_20201115.tsv",sep="\t")

}

# Then for the quantitative traits, which are a bit different
N0 <- "N"
rm(N1)
priors <- c(0.08245557,0.09135172) # These are sd.priors, see "Compute_custom_priors_BMI_Height_CLEAN_20201215.R" to find out how they were calculated
data <- data[9:10]
tests <- c("RapidoPGS_multi,1e-3,0.1,customprior","RapidoPGS_multi,1e-4,0.01,customprior" )

for(i in 1:2){

mbm <- microbenchmark(rapidoPGS_mult(data[i], alpha.block= 1e-3, alpha.snp=0.1, sd.prior=priors[i]),
	       rapidoPGS_mult(data[i], alpha.block= 1e-4, alpha.snp=0.01, sd.prior=priors[i]),
	       times=1, unit="s")
mbm <- as.data.table(summary(mbm))
mbm[,expr:=gsub("../datasets/","", data[i])][,method:=tests][,run:="set_prior"]
results <- rbind(results,mbm)
fwrite(results, "../evaluations/Results_benchmark_RapidoPGS_all_setprior_20201115.tsv",sep="\t")

}

#####################################################
#####          RÁPIDOPGS-MULTI (AUTO PRIOR)   #######
#####################################################

datasets <- c("AST_Demenais_29273806_1-qcfilt.tsv.gz","RA_Okada_24390342_3-qcfilt.tsv.gz", "T1D_Cooper_doi101101120022_1-qcfilt.tsv.gz", "T2D_Scott_28566273_1-qcfilt.tsv.gz", "BRCA_Michailidou_29059683_1-qcfilt.tsv.gz","PRCA_Schumacher_29892016_1-qcfilt.tsv.gz","CAD_Nikpay_26343387_1-qcfilt.tsv.gz","MDD_Wray_29700475_1-qcfilt.tsv.gz","BMI_Locke_25673413_1-qcfilt.tsv.gz","HEIGHT_Wood_25282103_1-qcfilt.tsv.gz")
data  <- paste0(dspath,datasets)

 for(i in 1:8){
 
 mbm <- microbenchmark(rapidoPGS_mult(data[i], alpha.block= 1e-3, alpha.snp=0.1),
 	       rapidoPGS_mult(data[i], alpha.block= 1e-4, alpha.snp=0.01),
 	       times=1, unit="s")
 mbm <- as.data.table(summary(mbm))
 mbm[,expr:=gsub("../datasets/","", data[i])][,method:=tests][,run:="autoprior"]
 results <- rbind(results,mbm)
 fwrite(results, "../evaluations/Results_benchmark_RapidoPGS_all_autoprior_20201115.tsv",sep="\t")
 
 }

# Then for the quantitative traits, which are a bit different
N0 <- "N"
rm(N1)
priors <- c(0.08245557,0.09135172) # These are priors for var_b, see "Compute_custom_priors_BMI_Height_CLEAN_20201215.R" to find out how they were calculated
data <- data[9:10]
tests <- c("RapidoPGS_multi,1e-3,0.1,customprior","RapidoPGS_multi,1e-4,0.01,customprior")

for(i in 1:2){

mbm <- microbenchmark(rapidoPGS_mult(data[i], alpha.block= 1e-3, alpha.snp=0.1, sd.prior=priors[i]),
	       rapidoPGS_mult(data[i], alpha.block= 1e-4, alpha.snp=0.01, sd.prior=priors[i]),
	       times=1, unit="s")
mbm <- as.data.table(summary(mbm))
mbm[,expr:=gsub("../datasets/","", data[i])][,method:=tests]
results <- rbind(results,mbm)
fwrite(results, "../evaluations/Results_benchmark_RapidoPGS_all_autoprior_20201215.tsv",sep="\t")

}




