#!/usr/bin/env Rscript

## script to perform colocalization for two files
## in this case between a IFC and aneQTL

rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)

## set the working directory 
setwd("/analysis_dir/")
## script to perform colocalization for two files
## in this case between a IFC and aneQTL


## load some packages
require(data.table)
require(coloc)

##format eQTL data
locus_name <- args[1]
print(locus_name)
rsid <-  args[2]
print(rsid)
chromosome <- args[3]
print(chromosome)
rsid_pos_38 <- args[5]
print(rsid_pos_38)
rsid_pos_37 <- args[4]
print(rsid_pos_37)
start_b37 <- as.numeric(rsid_pos_37) - 200000
start_b38 <- as.numeric(rsid_pos_38) - 200000
end_b37 <- as.numeric(rsid_pos_37) + 200000
end_b38 <- as.numeric(rsid_pos_38) + 200000

sink(paste0(locus_name,"_EUR.log"), type=c("output"),append=F)


tissues <- c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Pancreas","Muscle_Skeletal","Liver")
#tissues <- "Muscle_Skeletal"
for (tissue in tissues){

##format eQTL data
res.gene     <- paste0("./GTEx_Analysis_v8_EUR_eQTL_all_associations_csv/",tissue,".v8.EUR.allpairs.chr",chromosome,".csv")



res.gene     <- data.frame(fread(res.gene, sep=","))
names(res.gene)

library(tidyr)
library(dplyr)

res.gene <- res.gene %>% separate(variant_id,c("CHR","POS_38","A1","A2","build"))
head(res.gene)
tmp <- res.gene[,c("CHR","POS_38")]
names(tmp) <- c("chromosome","position")
tmp$position <- as.numeric(tmp$position)
tmp$position2 <- tmp$position +1 
tmp$chr_b38 <- tmp$chromosome
tmp$pos_b38 <- tmp$position
head(tmp)



range <- c(start_b38:end_b38)
range <- range[range>0]
res.gene <- res.gene[res.gene$POS_38 %in% range,]
write.table(tmp,"integrated_GWAS_GTEX_sm_LOinput.txt",col.names = F,quote=F,sep="\t",row.names=F)
system(paste0("liftOver integrated_GWAS_GTEX_sm_LOinput.txt hg38ToHg19.over.chain.gz integrated_GWAS_GTEX_sm_LO_b37.bed integrated_GWAS_GTEX_sm_LO_unmatched.bed"))
names_to_merge <- names(tmp)
tmp <- read.table("integrated_GWAS_GTEX_sm_LO_b37.bed",header=F,stringsAsFactors=F)
names(tmp) <- names_to_merge
tmp$CHR <- tmp$chromosome
tmp$POS <- tmp$position

tmp$CHR_POS <- paste(paste0(tmp$chromosome),tmp$position,sep = ":")
 tmp$CHR_POS_b38 <- paste(tmp$chromosome,tmp$pos_b38,sep=":")

res.gene$CHR_POS_b38 <- paste(res.gene$CHR,res.gene$POS,sep = ":")
head(res.gene)

res.gene <- res.gene[,c("A1","A2","slope","slope_se","pval_nominal","phenotype_id","CHR_POS_b38","maf")]
names(res.gene)<- c("ref","alt","beta","se","pvalue","molecular_trait_id","CHR_POS_b38","maf")
#res.gene$N <- 706
res.gene$N <- ifelse(tissue %in% "Muscle_Skeletal",706,ifelse(tissue %in% "Adipose_Subcutaneous",581,ifelse(tissue %in% "Pancreas",305,ifelse(tissue %in% "Liver",208,ifelse(tissue %in% "Kidney",73,ifelse(tissue %in% "Adipose_Visceral_Omentum",469,NA))))))

tmp <- tmp[,c("CHR_POS","CHR_POS_b38")]

tmp2 <- merge(res.gene,tmp,by="CHR_POS_b38")
library(tidyr)
tmp2 <- tmp2 %>% separate(CHR_POS,c("chromosome","position"),sep=":")
names(tmp2)
#some formatting to desired output:
tmp2 <- tmp2[,c("chromosome","position","ref","alt","maf","beta","se","pvalue","molecular_trait_id","N")]
#export
range <- c(start_b38:end_b38)
range <- range[range>0]
gene_list <- tmp2$molecular_trait_id[tmp2$pvalue < 1e-5]
gene_list <- unique(gene_list)
gene_list
for (gene in gene_list){
#extract gene
tmp3 <- tmp2[tmp2$molecular_trait_id %in% gene,]
tmp3 <- unique(tmp3)
print(paste0("running for:  ",gene, " Tissue:",tissue, "Locus: ",locus_name))
write.table(tmp3,paste0("./files/range_",gene,"_",tissue,".v8.EUR_chr",chromosome,".txt"),quote=F,row.names=F,sep="\t")
#chr_pos in build 37




res.qtl     <- data.frame(fread(paste0("./files/range_",gene,"_",tissue,".v8.EUR_chr",chromosome,".txt"), sep="\t"))

res.IFC     <- data.frame(fread(paste0(".inputcoloc_METAANALYSIS_IFC_adjBMI_EUR_DGC.txt"), sep="\t"))

#-----------------------------------------#
##--            combine both           --##
#-----------------------------------------#

## merge both files
res.qtl$CHR_POS <-paste(res.qtl$chromosome,res.qtl$position,sep="_")
res.qtl$CHR_POS <-gsub("chr","",res.qtl$CHR_POS)
res.IFC$CHR_POS <- paste(res.IFC$CHR,res.IFC$POS,sep="_")
res.coloc <- merge(res.qtl, res.IFC, by="CHR_POS", suffixes = c(".QTL", ".IFC"))
colSums(is.na(res.coloc))


## now perform the actual coloc analysis (assumes effects are coded by minor allele)
###dataset 1 = eQTL
##dataset 2=IFC
names(res.coloc)[20]<-"p.IFC"
res <- coloc.abf(dataset1=list(beta=res.coloc$beta.QTL, varbeta=res.coloc$se^2, varbeta=res.coloc$se^2,pvalues=res.coloc$pvalue, N=max(res.coloc$N), type="quant", MAF=res.coloc$maf.QTL),
                 dataset2=list(beta=res.coloc$beta.IFC, varbeta=res.coloc$StdErr^2,pvalues=res.coloc$p.IFC, N=max(res.coloc$TotalSampleSizeAll), type="quant", MAF=res.coloc$maf.IFC))
res$summary
 write.table(res$summary,paste0("./results/coloc_",locus_name,"_",gene,"_",tissue,".colocEUR.txt"),row.names=T,sep="\t",col.names = T,quote = F)
 res.coloc$log10p.IFC <- log10(res.coloc$p.IFC)
res.coloc$log10p.QTL <- log10(res.coloc$pvalue)
}
}

