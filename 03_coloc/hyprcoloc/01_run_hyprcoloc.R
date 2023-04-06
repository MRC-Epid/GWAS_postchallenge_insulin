#!/usr/bin/env Rscript

rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)

setwd("/analysis_dir/")

#library(hyprcoloc)
library(data.table)
## get the input - from the locus file
region <- args[1]
#print(paste0("region passed to R is",region))
chr.s  <- args[2]
#print(paste0("chr passed to R is",chr.s))
pos.s  <- args[3]
#print(paste0("Start position passed to R is",pos.s))
pos.e  <- args[4]
#print(paste0("End position passed to R is",pos.e))

phenotype   <- c("T2DMVP","FCI","FI","WHRadjBMI","Glu120","MS_ISI")

#-----------------------------------------#
##-- 	        import trait data      --##
#-----------------------------------------#

## loop to get all stats needed 
res <- NULL
res  <- lapply(phenotype, function(x){
  ## read the relevant data
 res            <- paste0("zcat /analysis_dir/input/",
                          x,"_MA_input.txt.gz",
                           " | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e,
                           " '$9==chr && $10>=low && $10<=upp {print $0}' > /analysis_dir/tmp.txt")
system(res)
   #print(res)
  res            <- data.table::fread("/analysis_dir/tmp.txt", sep = "\t", header = F, data.table = F) ### here is where it goes wrong? not subsetted?
#cat("number of rows for this phenotype: ",nrow(res),"\n")
names(res)  <-  c( "MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "Pvalue", "TotalSampleSize", "chr","pos")
  ## add names in harmonised way
res$Allele1 <- toupper(res$Allele1)
res$Allele2 <- toupper(res$Allele2)
res$MarkerName <- apply(res[, c("chr", "pos", "Allele2", "Allele1")], 1, function(x){

 return(paste0("chr", as.numeric(x[1]), "_", as.numeric(x[2]), "_", paste(sort(x[3:4]), collapse = "_")))

})

  ## add trait information to reshape the data afterwards
  res$pheno <- x
cat(class(res$MarkerName),"\n")
#cat(head(as.data.frame(res)),"\n")
head(res)
  ## return if lowest p meets threshold - at least 5e-4
   if(min(res$Pvalue) < 5e-4){
     return(res)
print(x)
print(paste0("total number of rows for one phenotype: ", nrow(res)))
   }
else{
  rm(res)
}
  return(res)
})


## reformat to an actual data frame
res     <- do.call(rbind, res)
print(paste0("Phenotypes included, survived p value filtering: "))
table(res$pheno) 
cat("total number of rows for all phenotypes: ", nrow(res))

## do tidyr
library(tidyr)
head(res)
res <- pivot_wider(res, id_cols = c("MarkerName",  "chr", "pos", "Allele1", "Allele2"),
                   values_from = c("Effect", "StdErr", "Pvalue", "TotalSampleSize"),
                   names_from = "pheno", names_sep = ".")
## convert to data frame
res <- data.frame(res)
names(res)
cat("with", nrow(res), "SNPs in total\n")

table(rowSums(!is.na(res))<29)
## omit NAs
res <- res[,colSums(!is.na(res))>0]
res <- na.omit(res)

cat("with", nrow(res), "SNPs in common\n")

#-----------------------------------------#
##-- 	         run HyPrColoc           --##
#-----------------------------------------#

require(hyprcoloc)
print("generate matrix for input")
## create matrices for input
betas           <- as.matrix(res[, grep("^Effect", names(res), value=T)])
ses             <- as.matrix(res[, grep("^StdErr", names(res), value=T)])
## add names
colnames(betas) <- colnames(ses) <- gsub("Effect\\.", "", colnames(betas))

## store SNP-IDs
snp.id          <- res$rsid
print("select prior grid")
## run over a grid of priors
if(ncol(betas) < 400){
  prior.grid      <- expand.grid(p1= 1e-4, p2 = c(0.95, 0.98, 0.99, 0.999), thr = c(0.5,0.6,0.7,0.8,0.9), stringsAsFactors = F)
}else{
  prior.grid      <- expand.grid(p1= 1e-4, p2 = .98, thr = .5, stringsAsFactors = F)
}
  

cat("run hyprcoloc over a grid of priors \n")

## run hypercoloc
res.coloc       <- lapply(1:nrow(prior.grid), function(x){
  print(prior.grid[x,])
  ## run hypcoloc with set of priors
  tmp <- hyprcoloc(betas, ses, trait.names = colnames(betas), snp.id=snp.id,
                   bb.selection = "align",
                   prior.1 = prior.grid$p1[x],
                   prior.2 = prior.grid$p2[x],
                   reg.thresh = prior.grid$thr[x],
                   align.thresh = prior.grid$thr[x],
                   snpscores = T)
  ## get results
  tmp     <- subset(tmp$results, traits != "None")
  if(nrow(tmp) > 0){
    ## add settings
    tmp$p1  <- prior.grid$p1[x]
    tmp$p2  <- prior.grid$p2[x]
    tmp$thr <- prior.grid$thr[x]
    return(tmp)
  }
})
## combine into one data frame
res.coloc <- do.call(rbind, res.coloc)

cat("do sensitivity analysis \n")


if(ncol(betas) < 400){
  res.sensi <- sensitivity.plot(betas, ses, trait.names = colnames(betas), snp.id=res$MarkerName,
                                bb.selection = "align",
                                reg.thresh = c(0.5,0.6,0.7,0.8,0.9),
                                align.thresh = c(0.5,0.6,0.7,0.8,0.9),
                                prior.2 = c(0.95, 0.98, 0.99, 0.999),
                                equal.thresholds = T,similarity.matrix = T)
}

   ## plot sensitivity analysis
 pdf(file = paste0("analysis_dir/results/hyprcoloc.locus.",region,".", chr.s,".",pos.s,".",pos.e,".pdf"),width = 7,height = 5)
  par(font =2)
  res.sensi[[1]]
  dev.off()

  ## store results
  write.table(res.sensi[[2]], paste0("/analysis_dir/results/hyprcoloc.sensimat.locus.",region,".", chr.s,".",pos.s,".",pos.e,".txt"), sep="\t", row.names=F)
  

cat("store general results \n")

## store the information
write.table(res.coloc, paste0("/analysis_dir/results/hyprcoloc.locus.",region,".", chr.s,".",pos.s,".",pos.e,".txt"), sep="\t", row.names=F)
