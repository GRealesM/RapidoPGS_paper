# This is a wrapper to semi-automatically retrieve datasets from GWAS catalog


library(data.table)

gwascat.download <- function(ID, hm_only=TRUE){
		gwc.manifest <- fread("https://www.ebi.ac.uk/gwas/api/search/downloads/studies_alternative")
		
		
		study.manifest <- gwc.manifest[ID == PUBMEDID,] 
		
		if(nrow(study.manifest) == 0) stop("Please provide a valid PUBMED ID")
		if(nrow(study.manifest) > 1){
			message("There are multiple datasets associated to this PUBMED ID")
			print(study.manifest[,"STUDY ACCESSION"])
			study_acc  <- readline(prompt=paste("Please select accession (from 1 to ", nrow(study.manifest),"):", sep=""))
			study_acc  <- as.integer(study_acc)
			while(!study_acc %in% 1:nrow(study.manifest)){
				message("Oops! Seems that your number is not in the options. Please try again.")
				print(study.manifest[,"STUDY ACCESSION"])
				study_acc  <- readline(prompt=paste("Please select accession (from 1 to ", nrow(study.manifest),"):", sep=""))
				study_acc  <- as.integer(study_acc)
			}
		}
		study.manifest <- study.manifest[study_acc,]
		url <- paste("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/", gsub(" ","", study.manifest$`FIRST AUTHOR`), "_", ID,"_",study.manifest$`STUDY ACCESSION`, "/harmonised/",ID,"-",study.manifest$`STUDY ACCESSION`,"-",sapply(strsplit(study.manifest$MAPPED_TRAIT_URI, "/"), `[`,5),".h.tsv.gz", sep="")
		if(!RCurl::url.exists(url)) stop("The file you requested is unavailable. This may be due to the fact that public and harmonised summary statistics do not exist. Please check at GWAS catalog website.")
		message("Retrieving dataset for ",study.manifest$`DISEASE/TRAIT`,", by ", study.manifest$`FIRST AUTHOR`,", from ",study.manifest$DATE, ", published at ", study.manifest$JOURNAL,", with accession ", study.manifest$`STUDY ACCESSION`,".")
		
		ds <- fread(url)
		if(hm_only){
		hmcols <- grep("hm_",names(ds), value=TRUE)
		ds  <- ds[,..hmcols]
		}
}

test1 <- gwascat.download(27863252) # Astle, should work
negcontrol <- gwascat.download(01223247236) # The Empress Pub phone number. Shouldn't work!
