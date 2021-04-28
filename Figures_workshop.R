# Code for Figures in rapidoPGS paper v3


library(data.table)
library(ggplot2)
#devtools::install_github("chr1swallace/seaborn_palettes", force=TRUE)
library(seaborn)
library(cowplot)
theme_set(theme_cowplot(font_size=10))
library(magrittr)


# Load data
res <- fread("../results/Full_results_table_20210331.tsv")
tim <- fread("../results/Benchmark_results_all_20210331.tsv")

t1d <- data.table(trait="T1D", adjusted=c("no","yes"), auc=0.5, r2=0, N=315306, cases=771, controls=314535, method="SBayesR_hm3") # T1D failed, but to still show the empty space in the plot, we need to do this "trick".
res <- rbind(res,t1d, fill=TRUE)
res <- res[order(trait, method)][,model:=NULL]
res <- res[!grepl("gwide", method)] # Remove T1D LDpred2-auto gwide, the only that worked in time.

res[trait=="RA"]
setnames(res, "trait","Trait")
cols <- (seaborn:::SEABORN_PALETTES$muted)[c(4,2,1,10,3)] # "#D65F5F" "#EE854A" "#4878D0" "#82C6E2" "#6ACC64"

##### Plot 1 - Adjusted, Case-control AUC hm3 (multi 1e-3, autoprior)
p1 <- copy(res)
p1 <- p1[adjusted=="yes"][!grepl("1kg", method)][!grepl("setprior", method)][!grepl("1e4001", method)][,label:=factor(rep(c("LDpred2_auto", "PRScs_auto", "RápidoPGS-multi","RápidoPGS-single", "SBayesR"), 10))][, label:=factor(label, levels(label)[c(1,2,5,3,4)])][, Trait:= factor(Trait, levels(factor(Trait))[c(3,7,1,8,9,10,4,6,2,5)])] 
cols <- cols[c(1,2,5,3,4)] # Re-order colours 

p1a <- p1[!is.na(cases)] ## Case-control only
unique(p1$Trait)

f1 <- ggplot(p1a,aes(fill=label, y=auc-0.5, x=Trait))+
  geom_bar(position="dodge", stat="identity", #alpha=0.8,
           colour = "black", size = 0.8) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_errorbar(aes(ymin=lci-0.5, ymax=hci-0.5), 
	         width=.2,
                 position=position_dodge(.9))+ 
  scale_fill_manual(values = cols)+
  ## coord_cartesian(ylim = c(0.5,0.9))+
  background_grid(major="y") +
  ## scale_y_continuous(minor_breaks = seq(0.5 , 0.9, 0.02), breaks = seq(0.5, 0.9, 0.1))+
  labs(fill = "Method")  +
  theme(legend.position="none")+
  scale_y_continuous("AUC",breaks=seq(0,0.4,by=0.1),
                     labels=seq(0.5,0.9,by=0.1))+ # fool hadley
  ylab("AUC")

##### Plot 2 - Adjusted, All r2 hm3 (multi 1e-3, autoprior)

f2 <- ggplot(p1, aes(fill=label, y=r2, x=Trait))+
  geom_bar(position="dodge", stat="identity", #alpha=0.8,
           colour = "black", size = 0.8) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_errorbar(aes(ymin=r2_lci, ymax=r2_hci), 
	         width=.2,
                 position=position_dodge(.9))+ 
  scale_fill_manual(values = cols)+
  ## coord_cartesian(ylim = c(0.5,0.9))+
  background_grid(major="y") +
  ## scale_y_continuous(minor_breaks = seq(0.5 , 0.9, 0.02), breaks = seq(0.5, 0.9, 0.1))+
  labs(fill = "Method")  +
  theme(legend.position="none")+
  scale_y_continuous(limits=c(0,0.27)) +
  ylab(bquote(r^2))

legend <- get_legend(
  # create some space to the left of the legend
  f1 + theme(legend.position="right", legend.box.margin = margin(0, 0, 0, 12)
  )
)

toprow <- plot_grid(f1, NULL, legend,
		    labels = c("A",""),
		    nrow=1, rel_widths=c(3,0.05, 0.6))
plot_grid(toprow,f2,
          labels=c("","B"),
          nrow=2)
ggsave("../figures/Fig3_20210331.png",width=220, height=100, units="mm")


##### Plot 3 - Timing for all methods (RápidoPGS as Plot 1, ie. 1e-3, 0.1, autoprior for multi)

p2  <- copy(tim)
p2 <- p2[!grepl("1kg", Test)][!grepl("setprior", Test)][!grepl("1e-4", Test)][,Test:=gsub("_hm3","",Test)][,Test:=gsub(",1e-3,0.1,autoprior","",Test)][,Test:=gsub("Rapido", "Rápido", Test)][,Time.minutes:=Time.seconds/60][, Test:=factor(Test, levels(as.factor(Test))[c(1,2,5,3,4)])][, Trait:= factor(Trait, levels(factor(Trait))[c(3,7,1,8,9,10,4,6,2,5)])][, Time.seconds:=ifelse(Time.seconds < 1000, paste0(formatC(Time.seconds, format = "f", big.mark = " ", digits=2), '"'), NA)]
p2

f3 <- ggplot(p2, aes(fill=reorder(Test, Time.minutes),
		     y=Time.minutes, x=forcats::fct_rev(Trait)))+
  scale_fill_manual(values = rev(cols[c(2,1,3,4,5)]),guide = guide_legend(reverse=TRUE, nrow=2, byrow=TRUE))+
  geom_bar(position="dodge", stat="identity", alpha=0.65,
           colour = "black", size = 0.8) +
  geom_text(aes(label= Time.seconds),position = position_dodge(0.9), size=3, color = "black", hjust=-0.1) +
  geom_hline(yintercept = 0, linetype="dashed")+
  labs(fill = "Method")  +
  ylab("Time (minutes)")+
  xlab("Trait")+
  scale_y_continuous(breaks=seq(0,460,20), limits=c(0,480))+
  background_grid(major="x") +
  theme(legend.position = "bottom" ,
	legend.title = element_text(size = 10),
        legend.text = element_text(size = 8, hjust=0))+
  coord_flip()

f3
ggsave("../figures/Fig4_20210406.png",width=180, height=180, units="mm")

##### Plot 4 - RápidoPGS performance comparison (adjusted, hm3) in AUC, of alpha methods and prior status
p3 <- copy(res)
p3 <- p3[adjusted=="yes"][grepl("RapidoPGS.*hm3", method, perl=TRUE)][!grepl("setprior", method)][,method:=factor(method)][, Trait:= factor(Trait, levels(factor(Trait))[c(3,7,1,8,9,10,4,6,2,5)])]
p3a <- p3[!is.na(cases)] ## Case-control only

fmuted <- (seaborn:::SEABORN_PALETTES$muted)
#fpastel <-  (seaborn:::SEABORN_PALETTES$pastel) 
#p3cols <- c(fmuted[1],fpastel[1], fmuted[5],fpastel[5],fmuted[3])
# [1] "#4878D0" "#A1C9F4" "#956CB4" "#D0BBFF" "#6ACC64"
p3cols <- fmuted[c(1,5,3)]
p3cols
# [1] "#4878D0" "#956CB4" "#6ACC64"

p3alabs <- c(
			expression(paste("RápidoPGS-multi,", alpha['block'],"=", 10^-3,", ", alpha['SNP'], " = ", 0.1, ", auto")),
#			expression(paste("RápidoPGS-multi,",alpha['block'],"=", 10^-3,", ", alpha['SNP'], " = ", 0.1, ", sd.prior = 0.2")),
			expression(paste("RápidoPGS-multi,",alpha['block'],"=", 10^-4,", ", alpha['SNP'], "=", 0.01, ", auto")),
#			expression(paste("RápidoPGS-multi,",alpha['block'],"=", 10^-4,", ", alpha['SNP'], "=", 0.01, ", sd.prior = 0.2")),
			"RápidoPGS-single, sd.prior = 0.2")

f4 <- ggplot(p3a,aes(fill=method, y=auc-0.5, x=Trait))+
  scale_fill_manual(values=p3cols, labels = p3alabs, guide = guide_legend(nrow=3, byrow=TRUE))+
  geom_bar(position="dodge", stat="identity", #alpha=0.8,
           colour = "black", size = 0.8) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_errorbar(aes(ymin=lci-0.5, ymax=hci-0.5), 
	         width=.2,
                 position=position_dodge(.9))+ 
  ## coord_cartesian(ylim = c(0.5,0.9))+
  background_grid(major="y") +
  ## scale_y_continuous(minor_breaks = seq(0.5 , 0.9, 0.02), breaks = seq(0.5, 0.9, 0.1))+
  labs(fill = "Method")  +
  theme(legend.position="none")+
#  theme(legend.position="bottom",legend.text = element_text(size = 8, hjust=0) )+
  scale_y_continuous("AUC",breaks=seq(0,0.4,by=0.1),
                     labels=seq(0.5,0.9,by=0.1))+ # fool hadley
  ylab("AUC")
f4

##### Plot 5 - RápidoPGS performance comparison (adjusted, hm3) in r^2, of alpha methods and prior status

f5labs <- c(
			expression(paste("RápidoPGS-multi,", alpha['block'],"=", 10^-3,", ", alpha['SNP'], " = ", 0.1, ", auto")),
#			expression(paste("RápidoPGS-multi,",alpha['block'],"=", 10^-3,", ", alpha['SNP'], " = ", 0.1, ", informed")),
			expression(paste("RápidoPGS-multi,",alpha['block'],"=", 10^-4,", ", alpha['SNP'], "=", 0.01, ", auto")),
#			expression(paste("RápidoPGS-multi,",alpha['block'],"=", 10^-4,", ", alpha['SNP'], "=", 0.01, ", informed")),
			"RápidoPGS-single, informed")

f5 <- ggplot(p3,aes(fill=method, y=r2, x=Trait))+
  scale_fill_manual(values=p3cols, labels = f5labs, guide = guide_legend(nrow=3, byrow=TRUE))+
  geom_bar(position="dodge", stat="identity", #alpha=0.8,
           colour = "black", size = 0.8) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_errorbar(aes(ymin=r2_lci, ymax=r2_hci), 
	         width=.2,
                 position=position_dodge(.9))+ 
  ## coord_cartesian(ylim = c(0.5,0.9))+
  background_grid(major="y") +
  ## scale_y_continuous(minor_breaks = seq(0.5 , 0.9, 0.02), breaks = seq(0.5, 0.9, 0.1))+
  labs(fill = "Method")  +
  theme(legend.position="bottom",legend.text = element_text(hjust=0) )+
  scale_y_continuous(limits=c(0,0.18))+
  ylab(bquote(r^2))
f5

toprow <- plot_grid(f4, NULL,
		    labels = c("A",""),
		    nrow=1, rel_widths=c(3, 0.66))
plot_grid(toprow,f5,
          labels=c("","B"),
          nrow=2)
ggsave("../figures/Fig1_20210401.png",width=150, height=120, units="mm")


##### Plot 6 - RápidoPGS performance comparison (adjusted, 1kg) in AUC, of alpha methods and prior status
p6 <- copy(res)
p6 <- p6[adjusted=="yes"][grepl("RapidoPGS.*1kg", method, perl=TRUE)][!grepl("setprior", method)][,method:=factor(method)][, Trait:= factor(Trait, levels(factor(Trait))[c(3,7,1,8,9,10,4,6,2,5)])]
p6a <- p6[!is.na(cases)] ## Case-control only


f6 <- ggplot(p6a,aes(fill=method, y=auc-0.5, x=Trait))+
  scale_fill_manual(values=p3cols, labels = p3alabs, guide = guide_legend(nrow=3, byrow=TRUE))+
  geom_bar(position="dodge", stat="identity", #alpha=0.8,
           colour = "black", size = 0.8) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_errorbar(aes(ymin=lci-0.5, ymax=hci-0.5), 
	         width=.2,
                 position=position_dodge(.9))+ 
  ## coord_cartesian(ylim = c(0.5,0.9))+
  background_grid(major="y") +
  ## scale_y_continuous(minor_breaks = seq(0.5 , 0.9, 0.02), breaks = seq(0.5, 0.9, 0.1))+
  labs(fill = "Method")  +
  theme(legend.position="none")+
  # theme(legend.position="bottom",legend.text = element_text(size = 8, hjust=0) )+
  scale_y_continuous("AUC",breaks=seq(0,0.4,by=0.1),
                     labels=seq(0.5,0.9,by=0.1))+ # fool hadley
  ylab("AUC")
f6


##### Plot 8 - Timing for all RápidoPGS methods (hm3)

p8  <- copy(tim)
p8 <- p8[grepl("RapidoPGS.*hm3", Test, perl = TRUE)][!grepl("1kg", Test)][!grepl("setprior", Test)][,Test:=gsub("_hm3","",Test)][,Time.minutes:=Time.seconds/60][, Test:=factor(Test, levels(as.factor(Test))[c(3,2,1)])][, Trait:= factor(Trait, levels(factor(Trait))[c(3,7,1,8,9,10,4,6,2,5)])][, Time.seconds:=ifelse(Time.seconds < 1000, paste0(formatC(Time.seconds, format = "f", big.mark = " ", digits=2), '"'), NA)]


f8cols <- rev(p3cols) # Due to the nature of the plot, we need to reverse the colours
f8labs <- rev(f5labs) # Same goes for the labs
f8 <- ggplot(p8, aes(fill=Test, y= Time.minutes, x=forcats::fct_rev(Trait)))+
  scale_fill_manual(values = f8cols,labels =f8labs, guide = guide_legend(reverse=TRUE, nrow=3, byrow=TRUE))+
  geom_bar(position="dodge", stat="identity", alpha=0.7,
           colour = "black", size = 0.8) +
  geom_text(aes(label= Time.seconds),position = position_dodge(0.9), size=3, color = "black", hjust=-0.1) +
  geom_hline(yintercept = 0, linetype="dashed")+
  labs(fill = "Method")  +
  ylab("Time (minutes)")+
  xlab("Trait")+
  scale_y_continuous(breaks=seq(0,180,20), limits=c(0,180))+
  background_grid(major="x") +
  theme(legend.position = "bottom" ,
	legend.title = element_text(size = 10),
        legend.text = element_text(size = 8, hjust=0))+
  coord_flip()

f8
ggsave("../figures/Fig2_20210407.png",width=180, height=180, units="mm")


###### Some extra computations 

## How faster is RápidoPGS-single?

rsbm <- p2[, .(fastest = min(Time.minutes), slowest = max(Time.minutes)), by=Trait][,timesfaster:=slowest/fastest]
rsbm
#      Trait    fastest  slowest timesfaster
#  1: Asthma 0.01700000 393.9793   23175.255
#  2:    BMI 0.01770000 267.4515   15110.254
#  3:   BRCA 0.03663333 466.0127   12721.001
#  4:    CAD 0.05023333 402.0060    8002.774
#  5: Height 0.01698333 186.0492   10954.809
#  6:    MDD 0.03446667 463.7068   13453.776
#  7:   PRCA 0.04790000 440.9133    9204.871
#  8:     RA 0.01136667 273.3043   24044.370
#  9:    T1D 0.01085000 155.1890   14303.134
# 10:    T2D 0.01651667 451.0030   27305.933
rsbm[,.(max(timesfaster), min(timesfaster))]
#          V1       V2
# 1: 27305.93 8002.774

rmbm <- rbind(p2[Test == "RápidoPGS_multi"], p2[,.SD[which.max(Time.minutes)], by=Trait]) 

rmbm <- rmbm[, .(multi = min(Time.minutes), slowest = max(Time.minutes)), by=Trait]
rmbm[, timesfaster:=slowest/multi]
#      Trait     multi  slowest timesfaster
#  1: Asthma  21.40987 393.9793   18.401765
#  2:    BMI  21.01202 267.4515   12.728502
#  3:   BRCA  30.06263 466.0127   15.501392
#  4:    CAD  18.13695 402.0060   22.165028
#  5: Height 164.84118 186.0492    1.128657
#  6:    MDD  21.84338 463.7068   21.228708
#  7:   PRCA  22.89893 440.9133   19.254754
#  8:     RA  17.73243 273.3043   15.412681
#  9:    T1D  16.02418 155.1890    9.684675
# 10:    T2D  20.99872 451.0030   21.477646
rmbm[,.(max(timesfaster), min(timesfaster))]
#          V1       V2
# 1: 22.16503 1.128657
rmbm[!which.min(timesfaster)][, .(mean(timesfaster))]  
#          V1
# 1: 17.31724

# How faster is 1e-04 than 1e-03

rmbm2 <- p8[!grepl("single", Test)]
rmbm2 <- rmbm2[, .(fastest = min(Time.minutes), slowest = max(Time.minutes)), by=Trait][,reduction:=1-(fastest/slowest)]
rmbm2
#      Trait    fastest   slowest reduction
#  1: Asthma   7.791683  21.40987 0.6360704
#  2:    BMI  11.460983  21.01202 0.4545510
#  3:   BRCA  15.910950  30.06263 0.4707400
#  4:    CAD   7.813000  18.13695 0.5692219
#  5: Height 105.733400 164.84118 0.3585741
#  6:    MDD   7.817517  21.84338 0.6421105
#  7:   PRCA  11.712650  22.89893 0.4885067
#  8:     RA   8.700417  17.73243 0.5093501
#  9:    T1D   6.607550  16.02418 0.5876514
# 10:    T2D   8.490617  20.99872 0.5956602
rmbm2[,.(max(reduction), min(reduction))]
#           V1        V2
# 1: 0.6421105 0.3585741



# How does RápidoPGS compare to other methods?

#  Compare median difference in r2 of RápidoPGS methods vs other methods

rp <- res[adjusted == "yes"][!grepl("setprior",	method)][!grepl("1kg", method)]
rp <- rp[,.(Trait,method, r2)]
rp  <- rp[!(method == "SBayesR_hm3" & Trait == "T1D")]
rp2 <- dcast(rp, Trait ~ method)
rp2[,c("RapidoPGS-single_LDpred2","RapidoPGS-single_PRScs","RapidoPGS-single_SBayesR",
      "RapidoPGS-multi1e3_LDpred2","RapidoPGS-multi1e3_PRScs","RapidoPGS-multi1e3_SBayesR",
      "RapidoPGS-multi1e4_LDpred2","RapidoPGS-multi1e4_PRScs","RapidoPGS-multi1e4_SBayesR"):=
list(RapidoPGSsingle_hm3 - LDpred2auto_perchr_hm3, RapidoPGSsingle_hm3 -  PRScs_hm3, RapidoPGSsingle_hm3- SBayesR_hm3,
   RapidoPGSmulti_1e301_autoprior_hm3 - LDpred2auto_perchr_hm3, RapidoPGSmulti_1e301_autoprior_hm3 -  PRScs_hm3, RapidoPGSmulti_1e301_autoprior_hm3 - SBayesR_hm3,
   RapidoPGSmulti_1e4001_autoprior_hm3 - LDpred2auto_perchr_hm3, RapidoPGSmulti_1e4001_autoprior_hm3 -  PRScs_hm3, RapidoPGSmulti_1e4001_autoprior_hm3 - SBayesR_hm3)]
rp2[,type:=ifelse(Trait %in% c("BMI", "Height"), "q","cc")]
rp2 <- rp2[,-c(2:7)] # Exclude first columns
columns= 2:10 # Select the ones to take the median
rp3 <-  rp2[,lapply(.SD, median, na.rm=TRUE), .SDcols=columns, by=type]
rp4 <- rp2[,lapply(.SD, median, na.rm=TRUE), .SDcols=columns]
rp4 <- rbind(rp3,rp4, fill=T)

rp2[,type:=NULL]
fwrite(rp2, "../results/Table_S4_diff.tsv", sep="\t")
fwrite(rp4, "../results/Table_S4_mediandiff.tsv", sep = "\t")



# What about relative, rather than absolute differences?

rr <- res[adjusted == "yes"][!grepl("setprior",	method)][!grepl("1kg", method)]
rr <- rr[,.(Trait,method, r2)]
rr  <- rr[!(method == "SBayesR_hm3" & Trait == "T1D")]
rr2 <- dcast(rr, Trait ~ method)
rr2[,c("RapidoPGS-single_LDpred2","RapidoPGS-single_PRScs","RapidoPGS-single_SBayesR",
      "RapidoPGS-multi1e3_LDpred2","RapidoPGS-multi1e3_PRScs","RapidoPGS-multi1e3_SBayesR",
      "RapidoPGS-multi1e4_LDpred2","RapidoPGS-multi1e4_PRScs","RapidoPGS-multi1e4_SBayesR"):=
list(RapidoPGSsingle_hm3 / LDpred2auto_perchr_hm3, RapidoPGSsingle_hm3 /  PRScs_hm3, RapidoPGSsingle_hm3/ SBayesR_hm3,
   RapidoPGSmulti_1e301_autoprior_hm3 / LDpred2auto_perchr_hm3, RapidoPGSmulti_1e301_autoprior_hm3 /  PRScs_hm3, RapidoPGSmulti_1e301_autoprior_hm3 / SBayesR_hm3,
   RapidoPGSmulti_1e4001_autoprior_hm3 / LDpred2auto_perchr_hm3, RapidoPGSmulti_1e4001_autoprior_hm3 /  PRScs_hm3, RapidoPGSmulti_1e4001_autoprior_hm3 / SBayesR_hm3)]
rr2[,type:=ifelse(Trait %in% c("BMI", "Height"), "q","cc")]
rr2 <- rr2[,-c(2:7)] # Exclude first columns
columns= 2:10 # Select the ones to take the median
rr3 <-  rr2[,lapply(.SD, median, na.rm=TRUE), .SDcols=columns, by=type]
rr4 <- rr2[,lapply(.SD, median, na.rm=TRUE), .SDcols=columns]
rr4 <- rbind(rr3,rr4, fill=T)
rr2[,type:=NULL]
fwrite(rr2, "../results/Table_S5_ratio.tsv", sep="\t")
fwrite(rr4, "../results/Table_S5_medianratio.tsv", sep = "\t")

# How do RápidoPGS HapMap3 and 1000 Genomes approaches compare?

rap <- res[adjusted == "yes"][grepl("RapidoPGS", method)][!grepl("setprior", method)]
hm3 <- rap[grepl("hm3", method)]
kg <- rap[grepl("1kg",method)]
all(gsub("_hm3", "", hm3$method) == gsub("_1kg", "", kg$method))
# [1] TRUE
all(hm3$Trait == kg$Trait)
# [1] TRUE
hm3[,diff_r2:=r2 - kg$r2]
hm3[,.(med=median(diff_r2)), by = method]
#                                 method          med
# 1:  RapidoPGSmulti_1e301_autoprior_hm3  0.002049484
# 2: RapidoPGSmulti_1e4001_autoprior_hm3  0.001685487
# 3:                 RapidoPGSsingle_hm3 -0.001926743


sim  <- as.data.table(system("wc -l ../models/*single*", intern=TRUE))
sim[, V1:=gsub("\\s+","", V1)]
sim <- sim[!grepl("total", V1)]
sim[,c("SNP","model"):=tstrsplit(V1, split="../models/", fixed=TRUE)]
sim[,V1:=NULL]
sim  <- sim[!grepl("RA_", model)]
sim[, SNP:= as.numeric(SNP) - 1] # remove header count
simhm3 <- sim[grepl("hm3", model)]
sim1kg <- sim[grepl("1kg", model)]
median(sim1kg$SNP - simhm3$SNP)
# [1] 3667300
median(sim1kg$SNP / simhm3$SNP)
# [1] 4.861469
simhm3
sim1kg


# How much did gwide T1D differed from per-chr?

fullres <- fread("../results/Full_results_table_20210331.tsv")
names(fullres)
t1d <- fullres[adjusted == "yes" & grepl("LDpred2", method) & trait == "T1D"]
t1d
#    trait adjusted       auc       lci       hci trait_prevalence
# 1:   T1D      yes 0.7661821 0.7460494 0.7843071      0.002445244
# 2:   T1D      yes 0.7537948 0.7345654 0.7720094      0.002445244
#            r2     r2_lci     r2_hci      N cases controls
# 1: 0.05430819 0.04929962 0.05863855 315306   771   314535
# 2: 0.05124990 0.04635865 0.05571962 315306   771   314535
#                                    model                 method
# 1:  T1D_hm3_LDpred2auto_gwide.full.model  LDpred2auto_gwide_hm3
# 2: T1D_hm3_LDpred2auto_perchr.full.model LDpred2auto_perchr_hm3
t1d[,.(-diff(r2))]
#             V1
# 1: 0.003058282

t1d$auc - t1d$lci
t1d$hci - t1d$auc
