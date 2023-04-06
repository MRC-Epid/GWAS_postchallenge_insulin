#!/usr/bin/env Rscript

## script to perform colocalization for two files
## in this case between post-challenge insulin traits and T2D

rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)

## parameterers for the region
chr.s <- as.numeric(args[1])
pos.s <- as.numeric(args[2])

## name of the data set to be used

## just to test
#chr.s <- 8
#pos.s <- 9184231


####################################################################################################################################

## Insulin Fold Change - 

##################################################################################################################################
setwd("/analysis_dir/")
#files to read in as arguments 
pheno <- "T2D" #pheno should be the file path to the phenotype coloc with the glycabolite - in this case just one phenotype so replaced with string and read in file directly
glyc   <- "IFC" # IFC in this case
## terminate if no input
if(is.na(chr.s) | is.na(pos.s) | is.na(glyc) | is.na(pheno)){
  print("terminate R due to missing input")
  q()
}

cat("----------------------------------\n")
cat("performing coloc for", glyc, "and", pheno, "on chromosome", chr.s, "at position", pos.s, "\n")
cat("----------------------------------\n")


## load some packages
require(data.table)
require(coloc)
res.ifc     <- paste0("IFC_MA_input.txt.gz")
res.ifc     <- data.frame(fread(res.ifc, sep="\t"))
## create column for chr and position
#res.ifc$id         <-NULL
res.ifc$p<- as.numeric(res.ifc$Pvalue)
## create MarkerName to ease merging with other input results
res.ifc$chr        <- as.numeric(gsub("X", 23, res.ifc$chr))


## subset to region of interest (reduce to 1MB region)
res.ifc     <- subset(res.ifc, chr == chr.s & (pos >= pos.s-5e5 & pos <= pos.s+5e5))



#gc()
## same for the T2D phenotype -- in this case MVP EUR T2D
res.T2D     <- data.frame(fread(paste0("MVP.T2D.EUR.MAF0.001.combined.dbGaP.txt.gz")))

names(res.T2D)<- c("rsid", "chr", "pos", "Allele1", "Allele2","Freq1", "Effect", "StdErr", "Pvalue" ,"TotalSampleSize","OR")
#"MarkerName"
res.T2D$chr        <- as.numeric(gsub("X", 23, res.T2D$chr))



#check names and change them if needed 
# Pvalue, chr,pos, Allele1, Allele2
res.T2D     <- subset(res.T2D, chr == chr.s & (pos >= pos.s-5e5 & pos <= pos.s+5e5))
res.T2D$MarkerName <- apply(res.T2D[,c("chr","pos", "Allele2", "Allele1")], 1, function(x){
  paste0("chr", as.numeric(x[1]), "_", as.numeric(x[2]), "_", paste(sort(x[3:4]), collapse = "_"))
})
head(res.T2D)
res.T2D$Pvalue <- as.numeric(res.T2D$Pvalue)
#if T2D P > 1e-5 quit 
p<- min(res.T2D$Pvalue)
p
if(p>1e-5){
  print("terminate R due to p value threshold: T2D")
  q()
}

## make alleles to upper
res.ifc$Allele1 <- toupper(res.ifc$Allele1)
res.ifc$Allele2 <- toupper(res.ifc$Allele2)
res.T2D$Allele1 <- toupper(res.T2D$Allele1)
res.T2D$Allele2 <- toupper(res.T2D$Allele2)


## keep only information needed (use p-values to perform colocalisation)
res.T2D <- res.T2D[, c("MarkerName", "chr", "pos", "Allele1", "Allele2", "Effect", "StdErr", "Pvalue", "Freq1","TotalSampleSize","rsid")]
res.T2D$Pvalue <- as.numeric(res.T2D$Pvalue)

#make the overlapping names matching so that the suffixes can be used further down the script. 
#coded so Allele1 is the effect allele
#-----------------------------------------#
##--            combine both           --##
#-----------------------------------------#

## merge both files

res.coloc <- merge(res.ifc, res.T2D, by="MarkerName", suffixes = c(".ifc", ".T2D"))
#names(res.coloc)
res.coloc$Pvalue.T2D <- as.numeric(res.coloc$Pvalue.T2D)
res.coloc$Pvalue.ifc <- as.numeric(res.coloc$Pvalue.ifc)

## drop value with zero standard error
res.coloc <- subset(res.coloc, !is.na(Effect.T2D) & !is.na(StdErr.T2D) & StdErr.T2D > 0) #for T2D
res.coloc <- subset(res.coloc, !is.na(Effect.ifc) & !is.na(StdErr.ifc) & StdErr.ifc > 0) #for IFC
## create log10p values for plotting
res.coloc$log10p.glyc  <- -pchisq((res.coloc$Effect.ifc/res.coloc$StdErr.ifc)^2, df=1, lower.tail=F, log.p=T)/log(10)
res.coloc$log10p.T2D <- -pchisq((res.coloc$Effect.T2D/res.coloc$StdErr.T2D)^2, df=1, lower.tail=F, log.p=T)/log(10)

# ## align alleles
 res.coloc$Effect.ifc.aligned <- ifelse(res.coloc$Allele1.T2D == res.coloc$Allele1.ifc, res.coloc$Effect.ifc, -res.coloc$Effect.ifc)
## create MAF column
 res.coloc$maf.T2D                <- ifelse(res.coloc$Freq1.T2D <= 0.5, res.coloc$Freq1.T2D, 1-res.coloc$Freq1.T2D) #T2D
 res.coloc$maf.ifc                <- ifelse(res.coloc$Freq1.ifc <= 0.5, res.coloc$Freq1.ifc, 1-res.coloc$Freq1.ifc) #ifc
## now perform the actual coloc analysis (assumes effects are coded by minor allele)


##run the coloc analyses
res <- coloc.abf(dataset1=list(beta=res.coloc$Effect.ifc, varbeta=res.coloc$StdErr.ifc^2,pvalues=res.coloc$Pvalue.ifc, N=max(res.coloc$TotalSampleSize.ifc), type="quant", MAF=res.coloc$maf.ifc),
                 ## obtained numbers from different cohorts
                 dataset2=list(beta=res.coloc$Effect.T2D, varbeta=res.coloc$StdErr.T2D^2,pvalues=res.coloc$Pvalue.T2D, N=max(res.coloc$TotalSampleSize.T2D), type="quant", MAF=res.coloc$maf.T2D))
#do some printing to the output of the run
  print(paste0("coloc performed between ",glyc," ",pheno))
print("coloc.abf results:")
print(res$summary)
print("coloc.abf priors:")
print(res$priors)

## write results to file
write.table(t(res$summary), 
            paste0("/analysis_dir/output/", glyc, ".", chr.s, ".", pos.s,".",pheno,"_MVP_IFC.coloc"), row.names=F, sep="\t", col.names=T, quote=F)



#---------------------------#
### Locus compare plots   ###
#---------------------------#

ifc <- res.coloc[,c("rsid","Pvalue.ifc")]
names(ifc) <- c("rsid","pval")

T2D <- res.coloc[,c("rsid","Pvalue.T2D")]
names(T2D) <- c("rsid","pval")
library(locuscomparer)
ifc <- ifc[!is.na(ifc$rsid),]
T2D <- T2D[!is.na(T2D$rsid),]
ifc$pval <- as.numeric(ifc$pval)
T2D$pval <- as.numeric(T2D$pval) 
png(paste0("/analysis_dir/plots/", glyc, ".", chr.s, ".", pos.s,".",pheno,"_MVP.png"),width=700,height=500)
locuscompare(ifc, T2D,title2 = 'T2D EUR', title1 = 'Insulin Fold Change adj.BMI EUR',population = "EUR")
dev.off()

