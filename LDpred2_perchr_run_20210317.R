# Computing polygenic scores using LDpred2-auto per-chr for all traits

# 2020/02/26
# We use 10 datasets (8 CC + 2 quantitative) which were previously QC'ed (see Preparing_datasets...R file), and aligned to HapMap3 manifest.
# We computed n_eff at QC step, so files should have everything they need.

# Load libraries and data
library(bigsnpr)
library(data.table)
## To make them reproducible!
set.seed(1)


ldpred2.auto.perchr <- function(sumstats, iterations=3000, initial_p=seq_log(1e-4,0.9,15), cores){
 
	message("Running LDpred2-auto per chromosome...")
	# Load map, which will help us trace indices in sumstats to those in matrices
	map  <- as.data.table(readRDS("../references/UKBB_mat/map.rds"))
	map[,SNPID:=paste(chr,pos, sep=":")]
	
	# Prepare sumstats
	setnames(sumstats, c("BETA", "SE"), c("beta", "beta_se")) 
#	if(!"SNPID" %in% names(sumstats))
		sumstats[,SNPID:=paste(CHR,BP,sep=":")]
	# Check if dataset is aligned to map, and align if not
	fullid.map <- paste(map$chr, map$pos, map$a0, map$a1, sep =":")
	fullid.sumstats <- paste(sumstats$CHR, sumstats$BP, sumstats$REF, sumstats$ALT, sep=":")
	if(!all(fullid.sumstats %in% fullid.map))
		stop("Sumstats not aligned to map. Please check.")
	
	results <- data.frame()
	
	# Main loop
	for(chrs in 1:22){
	    
	    message("Working on chromosome ", chrs)
	    sumstats.chr <- sumstats[CHR == chrs,]
	    map.chr <- map[chr == chrs,]
	    idx.chr  <- match(sumstats.chr$SNPID, map.chr$SNPID)
            df_beta <- as.data.frame(sumstats.chr[, c("beta", "beta_se", "n_eff")])
	    # Load matrix
	    corr <- readRDS(paste0("../references/UKBB_mat/LD_chr",chrs,".rds"))[idx.chr, idx.chr]
            ldsc <- snp_ldsc2(corr, df_beta) 
            h2_est <- ldsc[["h2"]]
	
	       if(h2_est < 1e-4){
			message("h2_est is too low, so I can't proceed. Skipping chr ", chrs,".")
			next
		}
            tmp <- tempfile(tmpdir = "tmp-data")
            corr <- bigsparser::as_SFBM(as(corr, "dgCMatrix"), tmp)
            on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

            message("Computing LDpred2 auto, with initial p = ",paste(initial_p, collapse=", "), ", initial h2 = ", h2_est,", and ",iterations," iterations (+1000 burn-in)")
     
             multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                                 vec_p_init = initial_p, num_iter=iterations,
                                 ncores = cores)
             beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
             sumstats.chr$weight <- rowMeans(beta_auto)
	     results <- rbind(results,sumstats.chr)

	} # End main loop
        
	setnames(results, c("beta","beta_se"), c("BETA","SE"))
	return(results)
}


#######################

# Here we're using an array job to run the models, so we'll create an array in bash containing
# the indexes to select the traits.
# And we'll use those to generate the models
args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args) # Array index will be a 1-10 integer.

traitlist <- c("Asthma","RA","T1D","T2D","BRCA","PRCA","CAD","MDD","BMI","Height", "RAint")
args <- traitlist[args] # This way we transform the index number into the real trait name

dataset <- paste0(args, "-hm3.tsv.gz")
cores <- 15 # Use all cores available (15 in our case)
setDTthreads(15)

ds <- fread(paste0("../datasets/",dataset))

ldpred2.model <- ldpred2.auto.perchr(ds, iterations=3000, initial_p=seq_log(1e-4,0.9,15), cores=cores)
fwrite(ldpred2.model, paste0("../models/",args, "_hm3_LDpred2auto_perchr.full.model"), sep="\t")

