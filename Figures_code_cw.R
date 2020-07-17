# Figures for rapidoPGS paper

library(data.table)
library(ggplot2)
#devtools::install_github("chr1swallace/seaborn_palettes")
library(seaborn)
library(cowplot)
theme_set(theme_cowplot(font_size=10))

# Load data
ds <- fread("Result_summary.tsv")
martinds <- fread("Martin_results.tsv", na.strings = "")

# Figure 1 - barplot showing LDPRED2 20k, PGS 20k (no recalc), PGS 4k (recalc).

fig1ds <- ds[ SNP.prefilt == "sAUCfilt" & Threshold == "20k"] #| Threshold == "1,00E-04" | Threshold == "1e-4"
## fig1sd <- fig1ds[-which(Method.PGS == "RapidoPGS" & Threshold != "20k" & Recalc =="N"),]
## fig1ds$Threshold <- gsub(pattern = "1,00E-04", "1e-4", fig1ds$Threshold)
fig1ds$Method.PGS <- gsub(pattern = "LDpred", "LDpred2-auto", fig1ds$Method.PGS)
fig1ds[,label:=paste(Method.PGS,Threshold#,"summaryAUC"
                   , sep="_")]
fig1ds[,label:=gsub(pattern = "1e-4", replacement = "trimmed", label)][,label:=gsub("_20k", "", label)]

# LDpred paper data, obtained from https://github.com/privefl/paper-ldpred2/blob/master/code/plot-res-ldpred.R
ldpred2paper.data <- data.table(Trait=c("Asthma", "BRCA","CAD","MDD","PRCA","RA", "T1D", "T2D"),
																AUC = c(0.582,0.656,0.621,0.589,0.702,0.597,0.766,0.639),
																label = "LDpred2-auto_Privé")
fig1ds <- rbindlist(list(fig1ds,ldpred2paper.data), fill = TRUE)
fig1ds[,label:=factor(label,
                       levels=c("LDpred2-auto_Privé","LDpred2-auto",
                                "RapidoPGS"## ,"RapidoPGS_trimmed"
                                ))]

## add blanks for Privé on the second row
addblank <- function(x, traits, col,val=0) {
  tmp <- x[Trait %in% traits & label=="LDpred2-auto"]
  tmp[,label:="LDpred2-auto_Privé"]
  tmp[[col]] <- val
  rbind(x,tmp)
}
fig1ds <- addblank(fig1ds, c("T1DnoHLA","BMI","Height"), "AUC", 0.5)

cols <- (seaborn:::SEABORN_PALETTES$muted)[c(4,2,1,3)]  #c("#023EFF", "#1AC938", "#E8000B","#FF7C00")

fig1bds <- fig1ds[Trait == "T1DnoHLA",]
fig1ads <- fig1ds[Trait != "T1DnoHLA",]

## basic plot
p <- ggplot() +
  scale_fill_manual(values = cols[1:4])+
  ## coord_cartesian(ylim = c(0.5,0.9))+
  background_grid(major="y") +
  ## scale_y_continuous(minor_breaks = seq(0.5 , 0.9, 0.02), breaks = seq(0.5, 0.9, 0.1))+
  labs(fill = "Method")  +
  theme(legend.position="none")
pab <- p +
  scale_y_continuous("summaryAUC",breaks=seq(0,0.4,by=0.1),
                     labels=seq(0.5,0.9,by=0.1)) # fool hadley


fig1a <- pab +
  geom_bar(aes(fill=label, y=AUC-0.5, x=Trait),
           data=fig1ads,
           position="dodge", stat="identity", #alpha=0.8,
           colour = "black", size = 0.8) +
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("AUC")
fig1a

## extract the legend
legend <- get_legend(
  # create some space to the left of the legend
  fig1a + theme(legend.position="right", legend.box.margin = margin(0, 0, 0, 12))
)

fig1b <- pab +
  geom_bar(aes(fill=label, y=AUC-0.5, x=Trait),
           data=fig1bds,
           position="dodge", stat="identity", #alpha=0.8,
           colour = "black", size = 0.8)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("AUC")
fig1b

martinds[grep("AST_Demenais",PRS),Trait:="Asthma"]
martinds[PRS=="AST_Demenais_29273806_1_hapmap3qcfilt-1e4.prs.martin",Label:="RapidoPGS_1e-4_recalc"]
martinds[PRS=="AST_Demenais_29273806_1_LDpred2_auto_1e4_4000_hapmap3qcfilt.martin",Label:="LDpred2-auto_full"]
fig1cds <- martinds[,.(Trait, r2, AUC=auc,label=Label)]
## fig1cds <- na.omit(fig1cds)
fig1cds <- fig1cds[label == "LDpred2-auto_full" | #Label == "RapidoPGS_1e-4_recalc" |
                   label == "RapidoPGS_1e-4_recalc",]
fig1cds[,label:=gsub("_1e-4_recalc|_full","",label)]
## fig1cds[,label:=sub("_full","",label)]
fig1cds  <- addblank(fig1cds, c("BMI","Height"), "r2")
fig1cds  <- addblank(fig1cds, c("Asthma"), "AUC",val=fig1ds[Trait=="Asthma" & label=="LDpred2-auto_Privé"]$AUC)
fig1cds
fig1cds[,label:=factor(label,
                       levels=c("LDpred2-auto_Privé","LDpred2-auto",
                                "RapidoPGS","RapidoPGS_trimmed"))]


fig1c.cols <- c("#5c82ff","#023eff","#5dea74","#1ac938","#ff4751","#e8000b","#ffab5c","#ff7c00","#a200b8","#390040", "#120014")

fig1c <- p +
  geom_bar(aes(fill=label, y=r2, x=Trait),
           data=fig1cds[Trait!="Asthma"],
           position="dodge", stat="identity",# alpha=0.8,
           colour = "black", size = 0.8)+
  ## coord_cartesian(ylim = c(0,0.072))+
  scale_y_continuous(limits=c(0,0.1)) +
  geom_hline(yintercept = 0.5, linetype="dashed")+
  ylab(bquote(r^2))
fig1c
  

fig1d <- p +
  geom_bar(aes(fill=label, y=AUC-0.5, x=Trait),
           data=fig1cds[Trait=="Asthma"],
           position="dodge", stat="identity",# alpha=0.8,
           colour = "black", size = 0.8)+
  scale_y_continuous("AUC",limits=c(0,0.4),breaks=seq(0,0.4,by=0.1),
                     labels=seq(0.5,0.9,by=0.1)) + # fool hadley
  geom_hline(yintercept = 0, linetype="dashed") 
fig1d
  

bottomrow <- plot_grid(fig1b,NULL,fig1d,NULL,fig1c,NULL,legend,
                       labels=c("B","","C","","D",""),
                       nrow=1,rel_widths=c(1.4,0.2,1.4,0.2,2,0.2,2))
plot_grid(fig1a,bottomrow,
          labels=c("A",""),
          nrow=2)

ggsave("Fig1.png",width=174, height=100, units="mm")

## gridExtra::grid.arrange(fig1a, gridExtra::arrangeGrob(fig1b, fig1c, ncol = 2), nrow=2)

################################################################################

# Figure 2 - plot of x=threshold (20k, 1e-4, 1e-3, 1e-2), y=rapido, linetype=recalc or not, facet=dataset

fig2ds <- ds[Method.PGS == "RapidoPGS" & SNP.prefilt == "sAUCfilt",]
fig2ds <- fig2ds[order(Trait,-Recalc),]
fig2ds[,column:=ifelse(1:.N <= 45, 1, 2)]
fig2ds[grep("e-",Threshold),Threshold:=paste0("w > ",Threshold)]
fig2ds[-grep("e-",Threshold),Threshold:=paste0("SNPs > ",Threshold)]

theme_set(theme_cowplot(font_size=10))
fig2 <- ggplot(fig2ds[order(Trait,nSNPs)],
               aes(y=AUC, x=nSNPs,
                   pch=Threshold,
                   colour = Recalc, group = Recalc))+
  geom_path()+
  geom_point() +
  scale_x_log10() +
  facet_wrap(~Trait, scales = "free_y",ncol=2)+
  scale_colour_seaborn(palette = "dark")+
  scale_shape_manual(values=c(19,21:24))+
  ylab("summaryAUC") +
  theme_classic()+
  theme(legend.position = "none")+
  guides(color = guide_legend(reverse = TRUE),
         shape=guide_legend(nrow=3,byrow=TRUE))
legend <- get_legend(fig2 + theme(legend.position="bottom",
                                  legend.box="vertical",
                                  legend.spacing=unit(0,"mm"),
                                  legend.key.size=unit(8,"pt"),
                                  legend.title=element_text(size=10),
                                  legend.text=element_text(size=8),
                                  legend.box.margin = margin()))
plot_grid(fig2) +
  draw_grob(legend, x=0.5, y=0.05, width=0.5, height=0.1,scale=0.5)

ggsave("Fig2.png",width=174, height=174, units="mm")

## fig2_second <- ggplot(fig2ds[46:81,], aes(y=AUC, x=Threshold, colour = Recalc, group = Recalc))+
## 	geom_line()+
## 	facet_grid(Trait~., scales = "free_y")+
## 	scale_colour_seaborn(palette = "dark")+
## 	theme_classic()+
## 	theme(axis.title.y = element_blank())+
## 	guides(color = guide_legend(reverse = TRUE))


## fig2_first <- ggplot(fig2ds[1:45,], aes(y=AUC, x=Threshold, colour = Recalc, group = Recalc))+
## 	geom_line()+
## 	facet_grid(Trait~., scales = "free_y")+
## 	scale_colour_seaborn(palette = "dark")+
## 	theme_classic()+
## 	theme(legend.position = "none")+
## 	guides(color = guide_legend(reverse = TRUE))

## fig2_second <- ggplot(fig2ds[46:81,], aes(y=AUC, x=Threshold, colour = Recalc, group = Recalc))+
## 	geom_line()+
## 	facet_grid(Trait~., scales = "free_y")+
## 	scale_colour_seaborn(palette = "dark")+
## 	theme_classic()+
## 	theme(axis.title.y = element_blank())+
## 	guides(color = guide_legend(reverse = TRUE))

## gridExtra::grid.arrange(fig2_first, fig2_second, ncol=2)

## Supplemental figure XX - BMI and Height x=threshold (20k, 1e-4, 1e-3, 1e-2), y=rapido, linetype=recalc or not, facet=dataset
SF2 <- martinds[,c("Trait", "r2", "Threshold","Recalc","Label")]
SF2 <- na.omit(SF2)

figS2 <- ggplot(SF2, aes(y=r2, x=Threshold, colour = Recalc, group = Recalc))+
	geom_line()+
	facet_grid(Trait~., scales = "free_y")+
	scale_colour_seaborn(palette = "dark")+
	theme_classic()+
	guides(color = guide_legend(reverse = TRUE))+
	ylab(bquote(r^2))
figS2

################################################################################

## Figure SXX - Benchmark RapidoPGS vs. LDpred2 with different iterations
bm <- data.table(Method=c("RápidoPGS", "LDpred2-auto_500", "LDpred2-auto_1000", "LDpred2-auto_3000"), time.seconds=c(18.03099,997.13944,1120.31012,1503.94528))
bm[,time.minutes:=time.seconds/60]
library(colorspace)

figBM <- ggplot(bm, aes(x=time.seconds, y=reorder(Method, time.seconds), fill=reorder(Method, time.seconds))) + 
  geom_bar(position="dodge", stat="identity", colour = "black", size = 0.8)+
  geom_text(aes(label=round(time.seconds,2)),hjust=-0.1,size=3) +
	#scale_fill_seaborn(palette = "bright6")+
  scale_fill_manual(values =c(cols[3],lighten(fig1c.cols[9:11],0.5)))+
  scale_x_continuous(limits=c(0,2000))+
#	coord_cartesian(ylim = c(0.5,0.9))+
#	geom_hline(yintercept = 0.5, linetype="dashed")+
  theme_cowplot()+
  background_grid(major="x") +
  theme(#panel.grid.major.y = element_line(),
        #panel.grid.minor = element_line(),
    ## axis.title.x = element_blank(),
    axis.text=element_text(size=8),
        legend.position = "none")+
#	scale_y_continuous(minor_breaks = seq(0.5 , 0.9, 0.02), breaks = seq(0.5, 0.9, 0.1))+
	labs(fill = "Method",x="Time (seconds)",y="Method")
  ## coord_flip() 
figBM

ggsave("Fig3.png",height=85,width=85,units="mm")

## Other calculations

# Median SNPs for 1e-4 RapidoPGS
medianSNPs <- median(c(66360, 210638, 459019, 461582, 277034,203129, 260572,135493,231493))

# Median distance between RapidoPGS and LDpred2
k20dis <- ds[Threshold == "20k" & SNP.prefilt == "sAUCfilt",]
k20rap <- k20dis[Method.PGS == "RapidoPGS",]
k20ldpred <- k20dis[Method.PGS == "LDpred",]
median(k20ldpred$AUC - k20rap$AUC)
# Or a data.table alternative
k20dis[,.(Diff = diff(AUC)), .(Trait)][,median(Diff)]
k20dis[,.(Diff = diff(AUC)), .(Trait)][,.(min(Diff), max(Diff))]


# Extra stuff - Stuff I liked but didn't make to the main paper
fig1cds <- martinds[,c("Trait", "r2", "Label")]
fig1cds <- na.omit(fig1cds)

fig1c.cols <- c("#5c82ff","#023eff","#5dea74","#1ac938","#ff4751","#e8000b","#ffab5c","#ff7c00","#a200b8","#390040", "#120014")

fig1c <- ggplot(fig1cds, aes(fill=Label, y=r2, x=Trait)) + 
	geom_bar(position="dodge", stat="identity", alpha=0.8, colour = "black", size = 0.8)+
	#scale_fill_seaborn(palette = "bright6")+
	scale_fill_manual(values = fig1c.cols)+
	coord_cartesian(ylim = c(0,0.072))+
	#	geom_hline(yintercept = 0.5, linetype="dashed")+
	theme_classic()+
	theme(panel.grid.major.y = element_line(),panel.grid.minor = element_line())+
	#	scale_y_continuous(minor_breaks = seq(0.5 , 0.9, 0.02), breaks = seq(0.5, 0.9, 0.1))+
	labs(fill = "Method")
fig1c
