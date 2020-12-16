# This script will create RápidoPGS-single models using a modification of the original code to drop for 8 case-control and 2 quantitative trait datasets using single and multiple causal variant approach in rapidoPGS.
# 2020/12/10
# Note that the functions below are experimental, which will be added to RapidoPGS R package

# Let's load the libraries
# install.packages("RapidoPGS")
#library(RapidoPGS) # We comment out the package, since we'll need to define new functions here
library(data.table)
library(GenomicRanges)
# remotes::install_github("stephenslab/susieR")
#library(susieR) # We won't need this either
library(bigsnpr)

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
	ds <- ds[SE != 0,] # Remove SNPs with exactly zero SE (sign of problems in the data)
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
  ds <- ds[!is.na(ld.block),] # Remove SNPs without assigned block

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




###############################################
#### GENERATING THE PGS MODELS   ##############
###############################################


# Set path to the datasets and reference, and their names
dspath <- "../datasets/"
datasets <- c("AST_Demenais_29273806_1-qcfilt.tsv.gz","RA_Okada_24390342_3-qcfilt.tsv.gz", "T1D_Cooper_doi101101120022_1-qcfilt.tsv.gz", "T2D_Scott_28566273_1-qcfilt.tsv.gz", "BRCA_Michailidou_29059683_1-qcfilt.tsv.gz","PRCA_Schumacher_29892016_1-qcfilt.tsv.gz","CAD_Nikpay_26343387_1-qcfilt.tsv.gz","MDD_Wray_29700475_1-qcfilt.tsv.gz")
# RapidoPGS-single needs sample sizes, so we provide them here
#N0 <- c(107715,61565, 8828,132532,119078,61106,123504,113154)
#N1 <- c(19954,19234,5913,26676,137045,79148,60801,59851)
# Note, these sample sizes don't include MDD, which had per-variant sample size, so we do a little trick here (see below)
#ref  <- "../references/RapidoPGS-ref/"
#ncores <- bigparallelr::nb_cores()

####################################################################
##### RÁPIDOPGS-SINGLE MODEL CREATION (CASE-CONTROL DATASETS) ######
####################################################################

# Here we apply RápidoPGS-single (aka Wakefield approach) to all case-control datasets, specifying their control (N0) and case (N1) number.
# For case-control datasets, default sd.prior = 0.2 will apply


#sapply(seq_along(datasets), function(x){
sapply(c(2,6), function(x){
	filebasename <- strsplit(datasets[x], split="-")[[1]][1]
	ds <- fread(paste0(dspath,datasets[x]))
	wakefield.model  <- rapidopgs_single(ds,trait="cc") # We use RapidoPGS current function, which is equivalent to rapidopgs_single
	fwrite(wakefield.model, paste0("../models/",filebasename, "_qcfilt_RapidoPGSsingle2_prior02.full.model"), sep="\t")

	
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
	
#	susie.model  <- rapidopgs_multi(ds,reference=ref, ncores=ncores,iterations=1000, alpha.block = 1e-3, alpha.snp = 0.1, sd.prior=prior)
#	fwrite(susie.model, paste("../models/",filebasename, "_qcfilt_RapidoPGSmulti1e301_customprior.full.model",sep=""), sep="\t")
	
	N <- "N"
	wakefield.model <- rapidopgs_single(ds,N, trait = "quant", sd.prior = prior)
	fwrite(wakefield.model, paste("../models/",filebasename, "_qcfilt_RapidoPGSsingle2_customprior.full.model",sep=""), sep="\t")
	    })




###################################################################
##### DEBUGGING SITE   					     ######
###################################################################

# Something wrong in PRCA and RA computations lead to NA AUC in evaluations.
# I identified several blocks with no ppi nor weight in both datasets
# Load all functions first
library(magrittr)
mPRCA <- fread("../models/PRCA_Schumacher_29892016_1_qcfilt_RapidoPGSsingle2_prior02.full.model")
mPRCA[is.na(weight),ld.block] %>% unique()
#  [1]   10  181  346  575  581  882 1013 1105 1179 1264 1281 1300 1480 1503
# [15] 1683
# Multiple blocks are missing data here...
mRA <- fread("../models/RA_Okada_24390342_3_qcfilt_RapidoPGSsingle2_prior02.full.model")
mRA[is.na(weight),ld.block] %>% unique()
#    [1] 1019 1072 1012 1073 1074 1020 1075 1076 1077 1078 1079 1011 1080
#   [14] 1081 1082 1083 1084 1085 1086 1087 1021 1088 1089 1090 1091 1092
#   [27] 1093 1022 1094 1095 1023 1013 1024 1025 1026 1014 1027 1028 1029
#   [40] 1030 1031 1032 1033 1034 1015 1035 1036 1037 1038 1039 1040 1041
#   [53] 1016 1042 1043 1044 1045 1046 1047 1048 1017 1049 1050 1051 1052
#   [66] 1053 1054 1055 1056 1018 1057 1058 1059 1060 1061 1062 1063 1064
#   [79] 1065 1066 1067 1068 1069 1070 1071 1152 1103 1096 1153 1154 1155
#   [92] 1156 1157 1158 1159 1104 1160 1161 1162 1163 1164 1165 1166 1167
#  [105] 1168 1169 1097 1170 1171 1172 1105 1173 1174 1175 1176 1177 1178
#  [118] 1106 1179 1107 1108 1109 1110 1111 1112 1113 1114 1115 1116 1117
#  [131] 1118 1119 1098 1120 1121 1122 1123 1124 1099 1125 1126 1127 1128
#  [144] 1129 1130 1100 1131 1132 1133 1134 1135 1136 1101 1137 1138 1139
#  [157] 1140 1141 1102 1142 1143 1144 1145 1146 1147 1148 1149 1150 1151
#  [170] 1238 1188 1180 1239 1240 1241 1189 1242 1243 1244 1181 1245 1246
#  [183] 1247 1248 1249 1250 1251 1190 1252 1253 1254 1255 1256 1191 1257
#  [196] 1258 1259 1260 1261 1192 1193 1194 1195 1196 1197 1182 1198 1199
#  [209] 1200 1201 1202 1183 1203 1204 1205 1206 1207 1184 1208 1209 1210
#  [222] 1211 1212 1185 1213 1214 1215 1216 1217 1218 1186 1219 1220 1221
#  [235] 1222 1223 1224 1225 1226 1227 1228 1229 1230 1187 1231 1232 1233
#  [248] 1234 1235 1236 1237 1311 1312 1313 1314 1315 1316 1317 1318 1319
#  [261] 1320 1321 1322 1323 1262 1263 1264 1265 1266 1267 1268 1269 1270
#  [274] 1271 1272 1273 1274 1275 1276 1277 1278 1279 1280 1281 1282 1283
#  [287] 1284 1285 1286 1287 1288 1289 1290 1291 1292 1293 1294 1295 1296
#  [300] 1297 1298 1299 1300 1301 1302 1303 1304 1305 1306 1307 1308 1309
#  [313] 1310 1375 1376 1377 1378 1379 1324 1325 1326 1327 1328 1329 1330
#  [326] 1331 1332 1333 1334 1335 1336 1337 1338 1339 1340 1341 1342 1343
#  [339] 1344 1345 1346 1347 1348 1349 1350 1351 1352 1353 1354 1355 1356
#  [352] 1357 1358 1359 1360 1361 1362 1363 1364 1365 1366 1367 1368 1369
#  [365] 1370 1371 1372 1373 1374 1427 1428 1429 1380 1381 1382 1383 1384
#  [378] 1385 1386 1387 1388 1389 1390 1391 1392 1393 1394 1395 1396 1397
#  [391] 1398 1399 1400 1401 1402 1403 1404 1405 1406 1407 1408 1409 1410
#  [404] 1411 1412 1413 1414 1415 1416 1417 1418 1419 1420 1421 1422 1423
#  [417] 1424 1425 1426 1439 1430 1440 1441 1431 1442 1443 1444 1445 1446
#  [430] 1447 1448 1449 1450 1451 1452 1432 1453 1454 1433 1455 1456 1434
#  [443] 1457 1458 1459 1460 1435 1461 1462 1463 1464 1465 1466 1436 1467
#  [456] 1468 1469 1437 1470 1471 1472 1473 1474 1475 1476 1438 1477 1478
#  [469] 1479 1480 1481 1482 1483 1484 1492 1493 1485 1494 1495 1496 1497
#  [482] 1498 1499 1486 1500 1501 1502 1503 1504 1505 1506 1487 1507 1508
#  [495] 1509 1510 1511 1488 1512 1513 1514 1515 1516 1517 1489 1518 1519
#  [508] 1520 1521 1522 1523 1524 1525 1490 1526 1527 1528 1529 1530 1491
#  [521] 1537 1532 1538 1539 1531 1540 1533 1541 1542 1543 1544 1545 1546
#  [534] 1547 1548 1549 1550 1551 1552 1553 1534 1554 1555 1556 1557 1558
#  [547] 1559 1560 1561 1562 1563 1535 1564 1565 1566 1567 1568 1569 1570
#  [560] 1571 1572 1536 1573 1574 1575 1576 1577 1578 1587 1580 1588 1589
#  [573] 1590 1591 1592 1593 1594 1595 1581 1596 1597 1579 1598 1599 1582
#  [586] 1600 1601 1602 1603 1604 1605 1606 1607 1583 1608 1609 1610 1611
#  [599] 1612 1613 1614 1615 1616 1617 1584 1585 1586   61    7    1   62
#  [612]   63   64   65   66    8   67   68   69   70   71   72    9   73
#  [625]   10   74   11   75   76   77   78   79   80   81   82   83   84
#  [638]   85   12   86   87   88   89   91   92   93   94   13   95   96
#  [651]    2   97   98   99  100  101  102  103  104  105  106   14  107
#  [664]  108  109  110  111  112   15  113  114  115  116  117  118  119
#  [677]  120   16  121  122  123  124  125  126   17  127  128  129  130
#  [690]  131  132  133   18   19   20   21   22    3   23   24   25   26
#  [703]   27    4   28   29   30   31   32   33   34   35   36   37    5
#  [716]   38   39   40   41   42   43   44   45   46    6   47   48   49
#  [729]   50   51   52   53   54   55   56   57   58   59   60 1625 1619
#  [742] 1618 1626 1627 1628 1629 1630 1631 1632 1633 1634 1620 1635 1636
#  [755] 1637 1638 1639 1640 1621 1641 1642 1643 1644 1645 1646 1647 1648
#  [768] 1649 1650 1622 1651 1652 1653 1654 1655 1623 1624 1656 1657 1658
#  [781] 1659 1660 1661 1662 1663 1664 1665 1666 1667 1668 1669 1670 1671
#  [794] 1672 1673 1674 1675 1676 1677 1678 1679 1680 1681 1682 1683 1684
#  [807] 1685 1686 1687 1688 1689 1690 1691 1692 1693 1694 1695 1696 1697
#  [820] 1698 1699 1700 1701 1702 1703  191  139  134  192  193  140  194
#  [833]  195  196  197  198  199  200  201  202  203  204  141  205  206
#  [846]  207  208  209  210  211  212  213  214  215  216  217  218  142
#  [859]  219  220  221  222  223  224  225  226  227  228  229  230  231
#  [872]  143  232  233  234  235  236  237  238  239  240  135  241  242
#  [885]  243  244  245  144  246  247  248  249  145  250  251  252  253
#  [898]  254  255  256  257  146  258  259  260  261  262  263  264  265
#  [911]  266  267  268  269  147  270  271  272  273  274  275  276  277
#  [924]  148  149  150  151  152  153  154  155  156  157  158  136  159
#  [937]  160  161  162  163  164  165  166  167  168  169  170  137  171
#  [950]  172  173  174  175  176  177  178  138  179  180  181  182  183
#  [963]  184  185  186  187  188  190  340  285  278  341  342  343  344
#  [976]  345  286  346  347  348  349  350  351  352  353  354  355  356
#  [989]  357  358  287  359  360  361  362  363  364  365  279  366  367
# [1002]  368  288  369  370  371  372  373  374  375  376  377  289  378
# [1015]  379  380  381  382  383  384  385  386  387  290  388  389  390
# [1028]  391  392  393  394  291  395  396  397  398  399  292  293  294
# [1041]  295  296  297  298  299  280  300  301  302  303  304  305  306
# [1054]  307  308  281  309  310  311  312  313  314  315  316  317  318
# [1067]  282  319  320  321  322  323  324  325  283  326  327  328  329
# [1080]  330  331  332  333  334  284  335  336  337  338  339  465  410
# [1093]  401  466  400  467  468  411  469  470  471  472  473  474  475
# [1106]  476  477  478  412  479  480  481  482  483  484  485  486  487
# [1119]  488  489  490  491  492  413  493  494  495  402  496  414  497
# [1132]  498  499  500  501  415  502  503  504  505  506  507  508  509
# [1145]  510  511  416  512  513  514  515  516  517  518  519  417  520
# [1158]  521  418  419  420  421  422  423  403  424  425  426  427  428
# [1171]  429  430  404  431  432  433  434  435  405  436  437  438  406
# [1184]  439  440  441  442  443  444  445  407  446  447  448  449  450
# [1197]  408  451  452  453  409  454  455  456  457  458  459  460  461
# [1210]  462  463  464  530  580  523  522  581  582  583  584  585  586
# [1223]  587  531  588  589  590  591  592  593  594  595  596  597  598
# [1236]  599  600  601  532  602  603  604  605  606  607  608  533  609
# [1249]  610  611  612  613  614  615  616  617  618  619  620  621  622
# [1262]  623  534  624  625  626  627  628  629  630  535  631  536  524
# [1275]  537  538  539  540  541  542  543  544  525  545  546  547  548
# [1288]  549  550  526  551  552  553  554  555  556  557  527  558  559
# [1301]  560  561  562  563  564  528  565  566  567  568  569  570  571
# [1314]  529  572  573  574  575  576  577  578  579  640  698  632  699
# [1327]  700  701  641  702  703  704  705  706  707  708  709  642  710
# [1340]  711  712  713  714  715  716  717  718  643  719  720  721  722
# [1353]  723  724  725  633  726  727  644  728  729  730  731  732  733
# [1366]  734  735  736  737  738  739  740  741  742  743  645  646  647
# [1379]  648  649  634  650  651  652  653  654  655  656  657  658  659
# [1392]  660  635  661  662  663  664  665  666  636  667  668  669  670
# [1405]  671  672  673  674  637  675  676  677  678  638  679  680  681
# [1418]  682  683  684  639  685  686  687  688  689  690  691  692  693
# [1431]  694  695  696  697  804  755  744  745  805  806  807  808  809
# [1444]  810  811  812  756  813  814  815  816  817  818  819  757  820
# [1457]  821  822  823  824  825  826  746  827  828  758  829  830  831
# [1470]  832  833  834  835  836  837  759  838  839  840  841  842  760
# [1483]  761  762  747  763  764  765  766  748  767  768  769  770  771
# [1496]  772  773  774  775  776  749  777  778  779  780  781  750  782
# [1509]  783  784  751  785  786  787  752  788  789  790  791  792  753
# [1522]  793  794  795  796  797  798  799  754  800  801  802  803  910
# [1535]  855  843  911  912  844  856  913  914  915  916  917  918  857
# [1548]  919  920  845  921  922  923  924  925  926  927  928  929  930
# [1561]  858  931  932  933  859  934  935  936  860  861  862  863  864
# [1574]  846  865  866  867  868  847  869  870  871  872  873  874  848
# [1587]  875  876  849  877  878  879  880  850  881  882  851  883  884
# [1600]  885  886  887  888  852  889  890  891  892  893  894  895  896
# [1613]  897  898  899  853  900  901  902  903  904  905  906  854  907
# [1626]  908  909  985  946  937  986  987  988  989  938  947  990  991
# [1639]  992  993  994  995  996  997  998  948  999 1000 1001 1002 1003
# [1652] 1004 1005 1006 1007 1008 1009 1010  949  950  951  939  952  953
# [1665]  954  955  956  957  958  959  960  940  961  962  963  964  965
# [1678]  941  942  943  966  944  967  968  969  970  971  972  973  974
# [1691]  975  976  945  977  978  979  980  981  982  983  984
# Looks like in most LD.blocks there are missing SNPs, probably due to the big amount (~1M) SNPs with BETA and SE == 0.

# Investigating PRCA
filebasename <- strsplit(datasets[6], split="-")[[1]][1]
ds <- fread(paste0(dspath,datasets[6]))
   # Assign ld.blocks, in case they werent computed yet
   		blranges <- RapidoPGS::EUR_ld.blocks
   message("Assigning LD blocks...")
   snpranges <- GRanges(seqnames=paste("chr",ds$CHR, sep=""), ranges=IRanges(start=ds$BP, end=ds$BP, names=ds$SNPID), strand="*")
   ds[,ld.block:=findOverlaps(snpranges, blranges, select='last')]
  sd.prior=0.2
  pi_i = 1e-4
   ds[,ppi:=wakefield_pp(beta=BETA, se=SE,pi_i, sd.prior), by = "ld.block"][, weight:=ppi*BETA]
ds[is.na(weight),]
#        CHR       BP ALT REF    BETA     SE       P ALT_FREQ    n_eff
#     1:   1 12779560   T   C  0.0209 0.0112 0.06177   0.1778 137933.1
#     2:   1 12779563   T   C -0.0568 0.0304 0.06160   0.0200 137933.1
#     3:   1 12779618   T   C  0.0213 0.0110 0.05391   0.1806 137933.1
#     4:   1 12780021   A   T  0.0043 0.0090 0.63680   0.2952 137933.1
#     5:   1 12781142   T   C  0.0100 0.0092 0.27480   0.2740 137933.1
#    ---                                                              
# 67734:  22 22354559   T   C -0.0117 0.0199 0.55660   0.0594 137933.1
# 67735:  22 22355368   T   C -0.0139 0.0118 0.23700   0.1492 137933.1
# 67736:  22 22355496   A   G  0.0131 0.0116 0.25850   0.8464 137933.1
# 67737:  22 22355640   C   G -0.0137 0.0117 0.23930   0.1523 137933.1
# 67738:  22 22356241   G  GA -0.0148 0.0119 0.21680   0.1492 137933.1
#              SNPID ld.block ppi weight
#     1:  1:12779560       10 NaN    NaN
#     2:  1:12779563       10 NaN    NaN
#     3:  1:12779618       10 NaN    NaN
#     4:  1:12780021       10 NaN    NaN
#     5:  1:12781142       10 NaN    NaN
#    ---                                
# 67734: 22:22354559     1683 NaN    NaN
# 67735: 22:22355368     1683 NaN    NaN
# 67736: 22:22355496     1683 NaN    NaN
# 67737: 22:22355640     1683 NaN    NaN
# 67738: 22:22356241     1683 NaN    NaN

# As expected, some blocks simply fail, let's investigate one
summary(ds[ld.block == 10,])

ds10 <- ds[ld.block == 10,]
# Opening wakefield_pp function open
beta <- ds10$BETA
se <- ds10$SE

    z=beta/se
    # compute V
    V  <- se^2
    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF
    lABF = 0.5 * (log(1-r) + (r * z^2))
    ## tABF - to add one we create another element at the end of zero for which pi_i is 1
    tABF <- c(lABF,0) # Contains NA
    vpi_i<-c(rep(pi_i,length(lABF)),1) 
    sBF <- logsum(tABF + log(vpi_i)) # NaN
    exp(lABF+log(pi_i)-sBF) # all NaN

# What if we remove all SE == 0?
ds <- ds[SE != 0,]
ds[,ppi:=wakefield_pp(beta=BETA, se=SE,pi_i, sd.prior), by = "ld.block"][, weight:=ppi*BETA]
ds[is.na(weight),]
# Empty data.table (0 rows and 13 cols): CHR,BP,ALT,REF,BETA,SE...
# Removing SE == 0 solves the issue (apparently)

# Trying to fix RA, too
filebasename <- strsplit(datasets[2], split="-")[[1]][1]
ds <- fread(paste0(dspath,datasets[2]))
ds <- ds[SE != 0,]
   snpranges <- GRanges(seqnames=paste("chr",ds$CHR, sep=""), ranges=IRanges(start=ds$BP, end=ds$BP, names=ds$SNPID), strand="*")
   ds[,ld.block:=findOverlaps(snpranges, blranges, select='last')]
ds[,ppi:=wakefield_pp(beta=BETA, se=SE,pi_i, sd.prior), by = "ld.block"][, weight:=ppi*BETA]
ds[is.na(weight),]
# Empty data.table (0 rows and 13 cols): SNPID,CHR,BP,ALT,REF,P...


filebasename <- strsplit(datasets[6], split="-")[[1]][1]
ds <- fread(paste0(dspath,datasets[6]))
   # Assign ld.blocks, in case they werent computed yet
   		blranges <- RapidoPGS::EUR_ld.blocks
   message("Assigning LD blocks...")
   snpranges <- GRanges(seqnames=paste("chr",ds$CHR, sep=""), ranges=IRanges(start=ds$BP, end=ds$BP, names=ds$SNPID), strand="*")
   ds[,ld.block:=findOverlaps(snpranges, blranges, select='last')]
  sd.prior=0.2
  pi_i = 1e-4
	wakefield.model  <- rapidopgs_single(ds,trait="cc") # We use RapidoPGS current function, which is equivalent to rapidopgs_single





