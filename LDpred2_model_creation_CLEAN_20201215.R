# Computing polygenic scores using LDpred2 auto for all traits

# 2020/12/15
# This script follows the vignette (https://privefl.github.io/bigsnpr/articles/LDpred2.html) to generate PGS using LDpred2.
# We use 10 datasets (8 CC + 2 quantitative) which were previously QC'ed (see Preparing_datasets...R file).
# We computed n_eff at QC step, so files should have everything they need.

# NOTE: This is a clean version of the original code used to run LDpred2-auto.
# We splitted this code per dataset prior to submission to the HPC, due to runtime constraints.
# We also set here the paremeter match.min.prop=0.05 to snp_match function. This prevented errors from arising in T1D dataset,
# but has no other side effects on the other datasets.

# Load libraries and data
library(bigsnpr)
library(data.table)
## To make them reproducible!
set.seed(1)

## Split the runs by supplying the numbers of the datasets (1 to 10)
## Example run: Rscript --vanilla LDpred2_model_creation_clean_20201215.R 1 
args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)


computePGS.ldpred2.auto  <- function(sumstats, ancestry="EUR", iterations=3000, initial_p= seq_log(1e-4, 0.9, 15)){
	
	message("Computing PGS using LDpred2 auto.")
	NCORES <- nb_cores()

	# Check for minimum columns
    	mincol <- c("CHR","BP","SNPID","REF","ALT","BETA","SE","P","n_eff")
    	if (!all(mincol %in% names(sumstats))) 
    	    stop("All minimum columns should be present in the summary statistics da
tase	t. Please check, missing columns: ", 
    	        paste(setdiff(mincol, names(sumstats)), collapse = ", "))
	# We don't need all columns, so we can filter them by name
	sumstats <- sumstats[,c("CHR","BP","SNPID","REF","ALT","BETA","SE","P","n_eff")]
	# Dear snp_match require specific names, so let's abide
	names(sumstats)[c(1:2,4:7)]  <- c("chr","pos","a0","a1","beta", "beta_se")

	sumstats <- na.omit(sumstats)

	results <- data.frame()
	for(chrs in 1:22){
		
		message("Working on chromosome ", chrs)
		
		# First thing is to check if we already have .rds files for our chr (ie. if we have read the bed files before). If not, we'll read it. This will create a .bk file, but it won't be a problem since we won't call this function again.
		if(!file.exists(paste0("../references/LDpred2-ref/chr",chrs,".rds"))){
		snp_readBed(paste("../references/LDpred2-ref/chr",chrs,".bed", sep=""))
		}
		# Attach the "bigSNP" object in R session
		# This object contain SNP data from a single chromosome from 1000 genomes phase 1.
		# This object contains a matrix of genotypes (analogous to bed), as well as fam and map slots, PLINK format.
		obj.bigSNP <- snp_attach(paste("../references/LDpred2-ref/chr",chrs, ".rds", sep=""))

		# See how the file looks like
		#str(obj.bigSNP, max.level = 2, strict.width = "cut")

		# In this case, we'll focus only on individuals of european ancestry
		euridx  <- grep(ancestry, obj.bigSNP$fam$family.ID)
		
		#eurset  <- subset(obj.bigSNP, ind.row=euridx)
		#obj.bigSNP <- snp_attach(eurset)

		# Get aliases for useful slots
		G   <- obj.bigSNP$genotypes
		CHR <- obj.bigSNP$map$chromosome
		POS <- obj.bigSNP$map$physical.pos
		#y   <- obj.bigSNP$fam$affection - 1
		#NCORES <- nb_cores()

		sumstats.chr <- sumstats[sumstats$chr == chrs,]

		# Also, we must remove all positions with se = 0 (otherwise it will complain)
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
		info_snp.chr <- snp_match(sumstats.chr, map, match.min.prop=0.05)
		
		message("Original sumstat had ", nrow(sumstats.chr)," for chr",chrs,". After matching ", nrow(info_snp.chr)," remained, and will be used for further steps.")
		# Computing LDpred2 scores for one chromosome

		#Correlation

		# First, you need to compute correlations between variants. We recommend to use a window size of 3 cM (see ref).

		# snp_asGeneticPos 
		POS2 <- snp_asGeneticPos(CHR, POS, dir = "../references/LDpred2-ref", ncores = NCORES)

		## indices in info_snp
		#ind.chr <- which(info_snp$chr == 1)
		df_beta <- info_snp.chr[, c("beta", "beta_se", "n_eff")]

		## indices in G
		ind.chr2 <- info_snp.chr$`_NUM_ID_`
		message("Computing correlation matrix for chr",chrs)
	
		## Create some temp file to store the matrix
		temp_file  <-  tempfile()

		corr0 <- snp_cor(G, ind.col = ind.chr2, ind.row= euridx, ncores = NCORES, infos.pos = POS2[ind.chr2], size = 3 / 1000)
		corr <- as_SFBM(corr0, backingfile=temp_file)

# We need to estimate h2 from the data, using ldsc which performs a LDscore regression. 
		message("Using LD score regression to estimate h2.")
		ldsc <- snp_ldsc2(corr0, df_beta)
		h2_est  <- ldsc[[2]]
		message("Estimated h2:",h2_est)
		if(h2_est < 1e-4){
			message("h2_est is too low, so I can't proceed. Skipping chr ", chrs,".")
			# Remove temp files!
			unlink(paste0(temp_file,"*"))
			next
		}
		# We can also use the Automatic model to test many different h2 and p parameters and choose the best model
# takes a few minutes
		message("Computing LDpred2 auto, with initial p = ",paste(initial_p, collapse=", "), ", initial h2 = ", h2_est,", and ",iterations," iterations (+1000 burn-in)")
		multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est, vec_p_init = initial_p,  num_iter=iterations, ncores = NCORES)
		beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
		# We average the resulting betas
		info_snp.chr$beta_auto <- rowMeans(beta_auto)
		results <- rbind(results,info_snp.chr)
		rm(info_snp.chr, corr, corr0, POS2, obj.bigSNP)
		# Remove temp files!
		unlink(paste0(temp_file,"*"))
	}

	return(results)
}



###################################
#### GENERATING LDPRED2 MODELS ####
###################################

dspath <- "../datasets/"
datasets <- c("AST_Demenais_29273806_1-qcfilt.tsv.gz","RA_Okada_24390342_3-qcfilt.tsv.gz", "T1D_Cooper_doi101101120022_1-qcfilt.tsv.gz", "T2D_Scott_28566273_1-qcfilt.tsv.gz", "BRCA_Michailidou_29059683_1-qcfilt.tsv.gz","PRCA_Schumacher_29892016_1-qcfilt.tsv.gz","CAD_Nikpay_26343387_1-qcfilt.tsv.gz","MDD_Wray_29700475_1-qcfilt.tsv.gz","BMI_Locke_25673413_1-qcfilt.tsv.gz","HEIGHT_Wood_25282103_1-qcfilt.tsv.gz")

x <- args
	filebasename <- strsplit(datasets[x], split="-")[[1]][1]
	ds <- bigreadr::fread2(paste0(dspath,datasets[x]))
	
	ldpred2res <- computePGS.ldpred2.auto(ds)
	# Transform and rename as desired
	ldpred2res <- as.data.table(ldpred2res)
	ldpred2res <- ldpred2res[,c("chr","pos","a0","a1", "SNPID", "beta","beta_se","P","beta_auto")]
	names(ldpred2res) <- c("CHR","BP","REF","ALT","SNPID","BETA","SE","P","weight")
	fwrite(ldpred2res, paste("../models/",filebasename, "_qcfilt_LDpred2.full.model",sep=""), sep="\t")

