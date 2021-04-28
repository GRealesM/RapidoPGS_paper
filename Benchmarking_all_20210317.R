# Benchmarking different methods for all traits
# 2020/03/17

# Here we'll compare the timing for PGS generation using different rapidoPGS with a single or multiple cores
# Let's load the libraries
library(RapidoPGS)
#library(data.table)
#library(GenomicRanges)
## remotes::install_github("stephenslab/susieR")
#library(susieR)
library(bigsnpr)
## remotes::install_github("chr1swallace/coloc", ref="susie")
#library(coloc)


##########################
## DEFINE LDPRED2 ########
##########################


ldpred2.auto.perchr <- function(sumstats, iterations=3000, initial_p=seq_log(1e-4,0.9,15), cores){
 
	message("Running LDpred2-auto per chromosome...")
	# Load map, which will help us trace indices in sumstats to those in matrices
	map  <- as.data.table(readRDS("../references/UKBB_mat/map.rds"))
	map[,SNPID:=paste(chr,pos, sep=":")]
	
	# Prepare sumstats
        sumstats <- copy(sumstats)
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

##########################
#### BENCHMARKING     ####
##########################

# Here we're using an array job to run the models, so we'll create an array in bash containing
# the indexes to select the traits.
# And we'll use those to generate the models
args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args) # Array index will be a 1-10 integer.

traitlist <- c("Asthma","RA","T1D","T2D","BRCA","PRCA","CAD","MDD","BMI","Height", "RAint")
tt <- traitlist[args] # This way we transform the index number into the real trait name

dataset <- paste0(tt, "-hm3.tsv.gz")
dataset.1kg <- paste0(tt, "-1kg.tsv.gz")
ncores <- 15 # Use all cores available (15 in our case)
setDTthreads(0) # Allow all cores for data.table


ds <- fread(paste0("../datasets/",dataset))
ds.1kg <- fread(paste0("../datasets/",dataset.1kg))
N_LD=362320 # N individuals to generate UKBB LD matrices, provided by Privé et al.
LDmatrices <- "../references/UKBB_mat"
reference <- "../references/RapidoPGS-ref/"

set.seed(1)

# Check if the datasets are case-controls or quantitative, and adapt the parameters to it
if(tt == "BMI" || tt == "Height"){
 trait.type <- "quant"
 N <- "N" 
 # We set the sd priors for quantitative traits, for both RápidoPGS single and multi.
 if(tt == "BMI"){
   priors  <- 0.1272537 
 } else if(tt == "Height"){
   priors  <- 0.1415274
 }
} else {
  trait.type <- "cc"
  N  <-  NULL
  priors <- 0.2 # For RápidoPGS single and multiple (set prior)
}


tests <- c("SBayesR_hm3", "LDpred2_auto_perchr_hm3", "RapidoPGS_single_hm3", "RapidoPGS_multi1e-3,0.1,setprior_hm3","RapidoPGS_multi,1e-4,0.01,setprior_hm3", "RapidoPGS_multi1e-3,0.1,autoprior_hm3","RapidoPGS_multi,1e-4,0.01,autoprior_hm3","RapidoPGS_single_1kg","RapidoPGS_multi1e-3,0.1,setprior_1kg","RapidoPGS_multi,1e-4,0.01,setprior_1kg", "RapidoPGS_multi1e-3,0.1,autoprior_1kg","RapidoPGS_multi,1e-4,0.01,autoprior_1kg")

time.seconds <- rep(NA, 12)

#######################################################
## Measure tests    ###################################
#######################################################

message("Measuring SBayesR")

timing <- system.time({ system(paste0("./gctb_2.02_Linux/gctb --sbayes R --mldm ../references/SBayesR/ukbEURu_hm3_sparse_mldm_list.txt --pi 0.95,0.02,0.02,0.01 --gamma 0.0,0.01,0.1,1 --gwas-summary ../datasets/",tt,"_SBayesR.ma --chain-length 10000  --burn-in 2000 --out-freq 10 --impute-n --p-value 0.5 --rsq 0.99 --out ../models/",tt,"_SBayesR_benchmark"))})

time.seconds[1] <- timing[3]


message("Measuring LDpred2-auto")
	
timing <- system.time({ldpred2.model <- ldpred2.auto.perchr(ds, iterations=3000, initial_p=seq_log(1e-4,0.9,15), cores=ncores)})

time.seconds[2] <- timing[3]


message("Measuring RápidoPGS-single (HapMap3)")

timing <- system.time({wakefield.model  <- rapidopgs_single(ds, N=N, trait=trait.type, sd.prior=priors)})

time.seconds[3] <- timing[3]


message("Measuring RápidoPGS-multi 1e03,0.1 set prior (HapMap3)")

timing <- system.time({susie1e301 <- rapidopgs_multi(ds, trait=trait.type, reference=NULL, LDmatrices=LDmatrices, N_LD = N_LD, N=N, ncores=ncores, alpha.block = 1e-3, alpha.snp = 0.1, sd.prior=priors) })

time.seconds[4] <- timing[3]


message("Measuring RápidoPGS-multi 1e-4,0.01 set prior (HapMap3)")

timing <- system.time({ susie1e4001  <- rapidopgs_multi(ds, trait=trait.type, reference=NULL,LDmatrices=LDmatrices, N_LD = N_LD, N=N, ncores=ncores, alpha.block = 1e-4, alpha.snp = 0.01, sd.prior=priors)})

time.seconds[5] <- timing[3]


message("Measuring RápidoPGS-multi 1e03,0.1 auto prior (HapMap3)")

timing <- system.time({susie1e301 <- rapidopgs_multi(ds, trait=trait.type, reference=NULL, LDmatrices=LDmatrices, N_LD = N_LD, N=N, ncores=ncores, alpha.block = 1e-3, alpha.snp = 0.1) })

time.seconds[6] <- timing[3]


message("Measuring RápidoPGS-multi 1e-4,0.01 auto prior (HapMap3)")

timing <- system.time({ susie1e4001  <- rapidopgs_multi(ds, trait=trait.type, reference=NULL,LDmatrices=LDmatrices, N_LD = N_LD, N=N, ncores=ncores, alpha.block = 1e-4, alpha.snp = 0.01)})

time.seconds[7] <- timing[3]


message("Measuring RápidoPGS-single (1000G)")

timing <- system.time({wakefield.model  <- rapidopgs_single(ds, N=N, trait=trait.type, sd.prior=priors)})

time.seconds[8] <- timing[3]


message("Measuring RápidoPGS-multi 1e03,0.1 set prior (1000G)")

timing <- system.time({susie1e301 <- rapidopgs_multi(ds, trait=trait.type, reference=reference, LDmatrices=NULL, N_LD = NULL, N=N, ncores=ncores, alpha.block = 1e-3, alpha.snp = 0.1, sd.prior=priors) })

time.seconds[9] <- timing[3]


message("Measuring RápidoPGS-multi 1e-4,0.01 set prior (1000G)")

timing <- system.time({ susie1e4001  <- rapidopgs_multi(ds, trait=trait.type, reference=reference,LDmatrices=NULL, N_LD = NULL, N=N, ncores=ncores, alpha.block = 1e-4, alpha.snp = 0.01, sd.prior=priors)})

time.seconds[10] <- timing[3]


message("Measuring RápidoPGS-multi 1e03,0.1 auto prior (1000G)")

timing <- system.time({susie1e301 <- rapidopgs_multi(ds, trait=trait.type, reference=reference, LDmatrices=NULL, N_LD = NULL, N=N, ncores=ncores, alpha.block = 1e-3, alpha.snp = 0.1) })

time.seconds[11] <- timing[3]


message("Measuring RápidoPGS-multi 1e-4,0.01 auto prior (1000G)")

timing <- system.time({ susie1e4001  <- rapidopgs_multi(ds, trait=trait.type, reference=reference,LDmatrices=NULL, N_LD = NULL, N=N, ncores=ncores, alpha.block = 1e-4, alpha.snp = 0.01)})

time.seconds[12] <- timing[3]

results <- data.table(Trait=tt, Test=tests, Cores=15, Time.seconds=time.seconds)

fwrite(results, paste0("../results/Benchmark_all_",tt,".tsv"),sep="\t") 

