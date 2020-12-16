# From RDS to a beautiful table

library(data.table)

ds <- readRDS("~/rds/rds-cew54-wallace-share/Projects/rapidoPGS/output/auc_r2/all.summary.rds")

quant <- ds[grep("BMI|HEIGHT",names(ds))]
quant <- lapply(quant, function(x) gsub("r\\^2: ","",x))
ds_quant <- data.frame(model=names(quant), r2 = sapply(quant,as.numeric))
ds_quant  <- ds_quant[order(ds_quant$model),]

cc <- ds[grep("BMI|HEIGHT", invert=TRUE, names(ds))]
ds_cc  <- data.frame(model = names(cc), auc = sapply(cc,`[`,1))
ds_cc  <- ds_cc[order(ds_cc$model),]

ds_quant <- as.data.table(ds_quant)
ds_cc <- as.data.table(ds_cc)

ds_res <- rbind(ds_cc, ds_quant, use.names=TRUE, fill=TRUE)
ds_res[,PGS.method:=ifelse(grepl("RapidoPGSmulti",model),"RápidoPGS-multi", ifelse(grepl("RapidoPGSsingle", model), "RápidoPGS-single", "LDpred2-auto"))][, alpha.block:=ifelse(grepl("1e3",model),1e-3, ifelse(grepl("1e4", model), 1e-4, NA))][, alpha.snp:=ifelse(grepl("1e301",model),0.1, ifelse(grepl("1e4001", model), 0.01, NA))][, sd.prior:=ifelse(grepl("autoprior", model), "auto", ifelse(grepl("prior02", model), "0.2", NA))]
ds_res[grepl("BMI", model) & grepl("custom", model), sd.prior:="0.08245557"]
ds_res[grepl("HEIGHT", model) & grepl("custom", model), sd.prior:="0.09135172"]
setnames(ds_res, "auc", "AUC")

qc1 <- fread("~/rds/rds-cew54-wallace-share/Projects/rapidoPGS/output/model_QC/summary_failed.txt")
# This is for table S1
#qc2 <- fread("~/rds/rds-cew54-wallace-share/Projects/rapidoPGS/output/model_QC/summary_N_inds.txt")
ds_res <- merge(ds_res, qc1, by="model")

fwrite(ds_res, "../evaluations/Full_results_table_20201215.tsv", sep = "\t")
