# Computing custom priors for BMI and Height
# 2020/11/12
# The previous version contained an error in the formula, so we include here a revised version

library(data.table)

##' Compute Standard deviation prior (SD prior) for quantitative traits
##'  using pre-computed heritability.
##' 
##' \code{sd.prior.est} function will take the dataset as an input, a \eqn{h^2} 
##' value obtained from a public repository such as LDhub, 
##' (http://ldsc.broadinstitute.org/ldhub/), sample size and number of variants,
##' and will provide a sd.prior estimate that can be used to improve prediction
##' performance of RápidoPGS functions on quantitative traits.
##' 
##' @param data a data.table containing the GWAS summary statistic input dataset. 
##' Must contain SNPID and SE columns.
##' @param h2 a numeric. Heritability estimate or h^2 (See details).
##' @param pi_i a numeric. Prior that a given variant is causal. DEFAULT = 1e-4.
##' @export
##' @author Guillermo Reales, Elena Vigorito, Chris Wallace
##' @examples 
##' sumstats <- data.table(SNPID=c("4:1479959","20:13000913","14:29107209","2:203573414",
##' "4:57331393","6:11003529","6:149256398","21:25630085","13:79166661"), 
##'			REF=c("C","C","C","T","G","C","C","G","T"), 
##'			ALT=c("A","T","T","A","A","A","T","A","C"), 
##'			ALT_FREQ=c(0.2611,0.4482,0.0321,0.0538,0.574,0.0174,0.0084,0.0304,0.7528),
##'			BETA=c(0.012,0.0079,0.0224,0.0033,0.0153,0.058,0.0742,0.001,-0.0131),
##'			SE=c(0.0099,0.0066,0.0203,0.0171,0.0063,0.0255,0.043,0.0188,0.0074),
##'			P=c(0.2237,0.2316,0.2682,0.8477,0.01473,0.02298,0.08472,0.9573,0.07535)) 
##' sd.prior <- sd.prior.est(sumstats, h2 = 0.2456, N = 45658, pi_i=1e-4)
##' 
##' 
sd.prior.est <- function(data, h2, N, pi_i=1e-4){
  data <- copy(data)
  data <- unique(data, by="SNPID")
  if(is.character(N)){
    sample.size <- N
    var_b <- h2 / (pi_i * sum(1/(data[,get(sample.size)] * data$SE^2)))
  } else{
    var_b <- h2 / (pi_i * sum(1/(N * data$SE^2)))
  }
  return(sqrt(var_b))
  
}

## Custom prior for BMI

bmi <- fread("../datasets/BMI_Locke_25673413_1-qcfilt.tsv.gz")
bmi <- unique(bmi, by="SNPID") # So we ensure there are no duplicated variants at all
bmi <- bmi[ALT_FREQ >= 0.01,]
h2 <- 0.246297646 # heritability estimate (h^2) obtained from LDhub 

sd.prior.bmi <- sd.prior.est(bmi, h2,"N")
# 0.08245557

# Custom prior for Height

height <- fread("../datasets/HEIGHT_Wood_25282103_1-qcfilt.tsv.gz")
height <- unique(height, by="SNPID")
height <- height[ALT_FREQ >= 0.01,]
h2 <- 0.4622973 # heritability estimate (h^2) obtained from LDhub 

sd.prior.height <- sd.prior.est(height,h2,"N")
# 0.09135172  

# Variances
#> c(sd.prior.bmi, sd.prior.height)^2
#[1] 0.006798920 0.008345137



