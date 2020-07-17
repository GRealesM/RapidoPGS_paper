# Computing polygenic scores using LDpred2 auto for all traits

# This script follows the vignette (https://privefl.github.io/bigsnpr/articles/LDpred2.html)

# Load libraries and data
library(bigsnpr)
library(data.table)
library(ggplot2)

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


computePGS.ldpred2.auto  <- function(sumstats, N0,N1=NULL, ancestry="EUR", iterations=3000, initial_p=1e-4){
	
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

####### COMPUTING!

# BMI (Martin QCfilt)


results  <- computePGS.ldpred2.auto("../datasets/BMI_Locke_25673413_1-qcfilt.tsv.gz",N0="N")
results <- as.data.table(results)
results <- results[,c("chr","pos","a0","a1", "SNPID", "beta","beta_se","P","beta_auto")]
names(results) <- c("CHR19","BP19","REF","ALT","SNPID","BETA","SE","P","weight")
results.martin <- results[,c("SNPID","ALT","weight")]
results.20k.martin  <- results.martin[order(-rank(abs(weight))),][1:20000,]

fwrite(results.martin, "../models/BMI_Locke_25673413_1_LDpred2_auto_1e4_4000_qcfilt.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.20k.martin, "../models/BMI_Locke_25673413_1_LDpred2_auto_1e4_4000_qcfilt-20k.martin",col.names=FALSE,sep="\t", quote=FALSE)


# HEIGHT

results  <- computePGS.ldpred2.auto("../datasets/HEIGHT_Wood_25282103_1-qcfilt.tsv.gz",N0="N")
results <- as.data.table(results)
results <- results[,c("chr","pos","a0","a1", "SNPID", "beta","beta_se","P","beta_auto")]
names(results) <- c("CHR19","BP19","REF","ALT","SNPID","BETA","SE","P","weight")
results.martin <- results[,c("SNPID","ALT","weight")]
results.20k.martin  <- results.martin[order(-rank(abs(weight))),][1:20000,]

fwrite(results.martin, "../models/HEIGHT_Wood_25282103_1_LDpred2_auto_1e4_4000_qcfilt.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.20k.martin, "../models/HEIGHT_Wood_25282103_1_LDpred2_auto_1e4_4000_qcfilt-20k.martin",col.names=FALSE,sep="\t", quote=FALSE)

# ASTHMA (sAUCfilt)

N0 <- 107715	
N1 <- 19954 

results  <- computePGS.ldpred2.auto("../datasets/AST_Demenais_29273806_1-sAUCfilt.tsv.gz", N0,N1)
results <- as.data.table(results)
results <- results[,c("chr","pos","a0","a1", "SNPID", "beta","beta_se","P","beta_auto")]
names(results) <- c("CHR19","BP19","REF","ALT","SNPID","BETA","SE","P","weight")
results.martin <- results[,c("SNPID","ALT","weight")]
results.20k.martin  <- results.martin[order(-rank(abs(weight))),][1:20000,]

results.sauc <- results[,pid:=paste(CHR19,BP19, sep=":")][,c("pid","ALT","weight")]
results.20k.sauc  <- results.sauc[order(-rank(abs(weight))),][1:20000,]

fwrite(results.martin, "../models/AST_Demenais_29273806_1_LDpred2_auto_1e4_4000_sAUCfilt.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.20k.martin, "../models/AST_Demenais_29273806_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.sauc, "../models/AST_Demenais_29273806_1_LDpred2_auto_1e4_4000_sAUCfilt.model",sep="\t", quote=FALSE)
fwrite(results.20k.sauc, "../models/AST_Demenais_29273806_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.model",sep="\t", quote=FALSE)

## ASTHMA (Martinfilt)

#N0 <- 107715	
#N1 <- 19954 
#
#results  <- computePGS.ldpred2.auto("../datasets/AST_Demenais_29273806_1-hapmap3qcfilt.tsv.gz",N0,N1)
#results <- as.data.table(results)
#results <- results[,c("chr","pos","a0","a1", "SNPID", "beta","beta_se","P","beta_auto")]
#names(results) <- c("CHR19","BP19","REF","ALT","SNPID","BETA","SE","P","weight")
#results.martin <- results[,c("SNPID","ALT","weight")]
#results.20k.martin  <- results.martin[order(-rank(abs(weight))),][1:20000,]
#
#results.sauc <- results[,pid:=paste(CHR19,BP19, sep=":")][,c("pid","ALT","weight")]
#results.20k.sauc  <- results.sauc[order(-rank(abs(weight))),][1:20000,]
#
#fwrite(results.martin, "../models/AST_Demenais_29273806_1_LDpred2_auto_1e4_4000_hapmap3qcfilt.martin",col.names=FALSE,sep="\t", quote=FALSE)
#fwrite(results.20k.martin, "../models/AST_Demenais_29273806_1_LDpred2_auto_1e4_4000_hapmap3qcfilt-20k.martin",col.names=FALSE,sep="\t", quote=FALSE)
#fwrite(results.sauc, "../models/AST_Demenais_29273806_1_LDpred2_auto_1e4_4000_hapmap3qcfilt.model",sep="\t", quote=FALSE)
#fwrite(results.20k.sauc, "../models/AST_Demenais_29273806_1_LDpred2_auto_1e4_4000_hapmap3qcfilt-20k.model",sep="\t", quote=FALSE)

# RA (sAUCfilt)

N0 <- 61565	
N1 <- 19234	

results  <- computePGS.ldpred2.auto("../datasets/RA_Okada_24390342_3-sAUCfilt.tsv.gz",N0,N1)
results <- as.data.table(results)
results <- results[,c("chr","pos","a0","a1", "SNPID", "beta","beta_se","P","beta_auto")]
names(results) <- c("CHR19","BP19","REF","ALT","SNPID","BETA","SE","P","weight")
results.martin <- results[,c("SNPID","ALT","weight")]
results.20k.martin  <- results.martin[order(-rank(abs(weight))),][1:20000,]

results.sauc <- results[,pid:=paste(CHR19,BP19, sep=":")][,c("pid","ALT","weight")]
results.20k.sauc  <- results.sauc[order(-rank(abs(weight))),][1:20000,]

fwrite(results.martin, "../models/RA_Okada_24390342_3_LDpred2_auto_1e4_4000_sAUCfilt.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.20k.martin, "../models/RA_Okada_24390342_3_LDpred2_auto_1e4_4000_sAUCfilt-20k.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.sauc, "../models/RA_Okada_24390342_3_LDpred2_auto_1e4_4000_sAUCfilt.model",sep="\t", quote=FALSE)
fwrite(results.20k.sauc, "../models/RA_Okada_24390342_3_LDpred2_auto_1e4_4000_sAUCfilt-20k.model",sep="\t", quote=FALSE)

## T1D (sAUCfilt)

N0 <- 8828
N1 <- 5913

results  <- computePGS.ldpred2.auto("../datasets/T1D_Cooper_1-sAUCfilt.tsv.gz",N0,N1)
results <- as.data.table(results)
results <- results[,c("chr","pos","a0","a1", "SNPID", "beta","beta_se","P","beta_auto")]
names(results) <- c("CHR19","BP19","REF","ALT","SNPID","BETA","SE","P","weight")
results.martin <- results[,c("SNPID","ALT","weight")]
results.20k.martin  <- results.martin[order(-rank(abs(weight))),][1:20000,]

results.sauc <- results[,pid:=paste(CHR19,BP19, sep=":")][,c("pid","ALT","weight")]
results.20k.sauc  <- results.sauc[order(-rank(abs(weight))),][1:20000,]

fwrite(results.martin, "../models/T1D_Cooper_1_LDpred2_auto_1e4_4000_sAUCfilt.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.20k.martin, "../models/T1D_Cooper_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.sauc, "../models/T1D_Cooper_1_LDpred2_auto_1e4_4000_sAUCfilt.model",sep="\t", quote=FALSE)
fwrite(results.20k.sauc, "../models/T1D_Cooper_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.model",sep="\t", quote=FALSE)

# T1D-noHLA (sAUCfilt)

set.seed(1)

N0 <- 8828
N1 <- 5913

results  <- computePGS.ldpred2.auto("../datasets/T1DnoHLA_Cooper_1-sAUCfilt.tsv.gz",N0,N1)
results <- as.data.table(results)
results <- results[,c("chr","pos","a0","a1", "SNPID", "beta","beta_se","P","beta_auto")]
names(results) <- c("CHR19","BP19","REF","ALT","SNPID","BETA","SE","P","weight")
results.martin <- results[,c("SNPID","ALT","weight")]
results.20k.martin  <- results.martin[order(-rank(abs(weight))),][1:20000,]

results.sauc <- results[,pid:=paste(CHR19,BP19, sep=":")][,c("pid","ALT","weight")]
results.20k.sauc  <- results.sauc[order(-rank(abs(weight))),][1:20000,]

fwrite(results.martin, "../models/T1DnoHLA_Cooper_1_LDpred2_auto_1e4_4000_sAUCfilt.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.20k.martin, "../models/T1DnoHLA_Cooper_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.sauc, "../models/T1DnoHLA_Cooper_1_LDpred2_auto_1e4_4000_sAUCfilt.model",sep="\t", quote=FALSE)
fwrite(results.20k.sauc, "../models/T1DnoHLA_Cooper_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.model",sep="\t", quote=FALSE)


## T2D (sAUCfilt)

N0 <- 132532
N1 <- 26676

results  <- computePGS.ldpred2.auto("../datasets/T2D_Scott_28566273_1-sAUCfilt.tsv.gz",N0,N1)
results <- as.data.table(results)
results <- results[,c("chr","pos","a0","a1", "SNPID", "beta","beta_se","P","beta_auto")]
names(results) <- c("CHR19","BP19","REF","ALT","SNPID","BETA","SE","P","weight")
results.martin <- results[,c("SNPID","ALT","weight")]
results.20k.martin  <- results.martin[order(-rank(abs(weight))),][1:20000,]

results.sauc <- results[,pid:=paste(CHR19,BP19, sep=":")][,c("pid","ALT","weight")]
results.20k.sauc  <- results.sauc[order(-rank(abs(weight))),][1:20000,]

fwrite(results.martin, "../models/T2D_Scott_28566273_1_LDpred2_auto_1e4_4000_sAUCfilt.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.20k.martin, "../models/T2D_Scott_28566273_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.sauc, "../models/T2D_Scott_28566273_1_LDpred2_auto_1e4_4000_sAUCfilt.model",sep="\t", quote=FALSE)
fwrite(results.20k.sauc, "../models/T2D_Scott_28566273_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.model",sep="\t", quote=FALSE)


## BRCA (sAUCfilt)

N0 <- 119078	
N1 <- 137045 

results  <- computePGS.ldpred2.auto("../datasets/BRCA_Michailidou_29059683_1-sAUCfilt.tsv.gz",N0,N1)
results <- as.data.table(results)
results <- results[,c("chr","pos","a0","a1", "SNPID", "beta","beta_se","P","beta_auto")]
names(results) <- c("CHR19","BP19","REF","ALT","SNPID","BETA","SE","P","weight")
results.martin <- results[,c("SNPID","ALT","weight")]
results.20k.martin  <- results.martin[order(-rank(abs(weight))),][1:20000,]

results.sauc <- results[,pid:=paste(CHR19,BP19, sep=":")][,c("pid","ALT","weight")]
results.20k.sauc  <- results.sauc[order(-rank(abs(weight))),][1:20000,]

fwrite(results.martin, "../models/BRCA_Michailidou_29059683_1_LDpred2_auto_1e4_4000_sAUCfilt.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.20k.martin, "../models/BRCA_Michailidou_29059683_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.sauc, "../models/BRCA_Michailidou_29059683_1_LDpred2_auto_1e4_4000_sAUCfilt.model",sep="\t", quote=FALSE)
fwrite(results.20k.sauc, "../models/BRCA_Michailidou_29059683_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.model",sep="\t", quote=FALSE)


## PRCA (sAUCfilt)

N0 <- 61106	
N1 <- 79148

results  <- computePGS.ldpred2.auto("../datasets/PRCA_Schumacher_29892016_1-sAUCfilt.tsv.gz",N0,N1)
results <- as.data.table(results)
results <- results[,c("chr","pos","a0","a1", "SNPID", "beta","beta_se","P","beta_auto")]
names(results) <- c("CHR19","BP19","REF","ALT","SNPID","BETA","SE","P","weight")
results.martin <- results[,c("SNPID","ALT","weight")]
results.20k.martin  <- results.martin[order(-rank(abs(weight))),][1:20000,]

results.sauc <- results[,pid:=paste(CHR19,BP19, sep=":")][,c("pid","ALT","weight")]
results.20k.sauc  <- results.sauc[order(-rank(abs(weight))),][1:20000,]

fwrite(results.martin, "../models/PRCA_Schumacher_29892016_1_LDpred2_auto_1e4_4000_sAUCfilt.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.20k.martin, "../models/PRCA_Schumacher_29892016_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.sauc, "../models/PRCA_Schumacher_29892016_1_LDpred2_auto_1e4_4000_sAUCfilt.model",sep="\t", quote=FALSE)
fwrite(results.20k.sauc, "../models/PRCA_Schumacher_29892016_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.model",sep="\t", quote=FALSE)



## CAD (sAUCfilt)

N0 <- 123504	
N1 <- 60801

results  <- computePGS.ldpred2.auto("../datasets/CAD_Nikpay_26343387_1-sAUCfilt.tsv.gz",N0,N1)
results <- as.data.table(results)
results <- results[,c("chr","pos","a0","a1", "SNPID", "beta","beta_se","P","beta_auto")]
names(results) <- c("CHR19","BP19","REF","ALT","SNPID","BETA","SE","P","weight")
results.martin <- results[,c("SNPID","ALT","weight")]
results.20k.martin  <- results.martin[order(-rank(abs(weight))),][1:20000,]

results.sauc <- results[,pid:=paste(CHR19,BP19, sep=":")][,c("pid","ALT","weight")]
results.20k.sauc  <- results.sauc[order(-rank(abs(weight))),][1:20000,]

fwrite(results.martin, "../models/CAD_Nikpay_26343387_1_LDpred2_auto_1e4_4000_sAUCfilt.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.20k.martin, "../models/CAD_Nikpay_26343387_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.sauc, "../models/CAD_Nikpay_26343387_1_LDpred2_auto_1e4_4000_sAUCfilt.model",sep="\t", quote=FALSE)
fwrite(results.20k.sauc, "../models/CAD_Nikpay_26343387_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.model",sep="\t", quote=FALSE)


# MDD (sAUCfilt)


results  <- computePGS.ldpred2.auto("../datasets/MDD_Wray_29700475_1-sAUCfilt.tsv.gz",N0="N0",N1="N1")
results <- as.data.table(results)
results <- results[,c("chr","pos","a0","a1", "SNPID", "beta","beta_se","P","beta_auto")]
names(results) <- c("CHR19","BP19","REF","ALT","SNPID","BETA","SE","P","weight")
results.martin <- results[,c("SNPID","ALT","weight")]
results.20k.martin  <- results.martin[order(-rank(abs(weight))),][1:20000,]

results.sauc <- results[,pid:=paste(CHR19,BP19, sep=":")][,c("pid","ALT","weight")]
results.20k.sauc  <- results.sauc[order(-rank(abs(weight))),][1:20000,]

fwrite(results.martin, "../models/MDD_Wray_29700475_1_LDpred2_auto_1e4_4000_sAUCfilt.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.20k.martin, "../models/MDD_Wray_29700475_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.martin",col.names=FALSE,sep="\t", quote=FALSE)
fwrite(results.sauc, "../models/MDD_Wray_29700475_1_LDpred2_auto_1e4_4000_sAUCfilt.model",sep="\t", quote=FALSE)
fwrite(results.20k.sauc, "../models/MDD_Wray_29700475_1_LDpred2_auto_1e4_4000_sAUCfilt-20k.model",sep="\t", quote=FALSE)

