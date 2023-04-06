#generate the scores and save in file for analysis input
## load the effects for weighting 
setwd("/analysis_dir/01_generate_scores")
library(data.table)
FCI <- fread("../input/FCI_stats.txt",header=T,stringsAsFactors=F,sep="\t")
MS_ISI <- fread("../input/MS_ISI__stats.txt",header=T,stringsAsFactors=F,sep = "\t")

#format weighting files 
names(FCI)
FCI <- FCI[,c(2,3,4,6,10)]
names(FCI)<- c("rsid","ea","oa","effect_adjBMI","effect_noBMI")
summary(FCI$effect_adjBMI)
FCI$ea <- toupper(FCI$ea)
FCI$oa <- toupper(FCI$oa)

names(MS_ISI)
MS_ISI <- MS_ISI[,c(2,3,4,6,10)]
names(MS_ISI)<- c("rsid","ea","oa","effect_adjBMI","effect_noBMI")
summary(MS_ISI$effect_adjBMI)
MS_ISI$ea <- toupper(MS_ISI$ea)
MS_ISI$oa <- toupper(MS_ISI$oa)
#remember to check if the alleles are aligned in the study and the weighting


# UKBB
ukb_FCI <- fread("FCI_UKBB.dosage.gz")
map.mono <- fread("FCI_UKBB.stats.gz")
head(map.mono)
map.mono$alleleB_dosage_code <- as.character(map.mono$alleleB_dosage_code)
map.mono <- merge(map.mono,FCI,by="rsid")
head(map.mono)
map.mono$weight_UA <- ifelse(map.mono$ea == map.mono$alleleB_dosage_code, map.mono$effect_noBMI, -map.mono$effect_noBMI)
map.mono$weight_Adj <- ifelse(map.mono$ea == map.mono$alleleB_dosage_code, map.mono$effect_adjBMI, -map.mono$effect_adjBMI)
snps <- ukb_FCI

id <- ukb_FCI$IID
rsid <- map.mono$rsid
for (i in rsid){
  ukb_FCI[[i]]<- ukb_FCI[[i]]*map.mono$weight_UA[map.mono$rsid %in% i]
}

ukb_FCI$FCI_noBMI <- NA
ukb_FCI$FCI_noBMI <- rowSums(ukb_FCI[,c(2:5)]) 
ukb_FCInoBMI <- ukb_FCI[,c("IID","FCI_noBMI")] 
ukb_FCInoBMI <- merge(ukb_FCInoBMI,snps,by = "IID")

id <- ukb_FCI$IID
rsid <- map.mono$rsid

ukb_FCI <- fread("FCI_MAGICRI_UKBB.dosage.gz")

for (i in rsid){
  ukb_FCI[[i]]<- ukb_FCI[[i]]*map.mono$weight_Adj[map.mono$rsid %in% i]
}

ukb_FCI$FCI_adjBMI <- NA
ukb_FCI$FCI_adjBMI <- rowSums(ukb_FCI[,c(2:5)]) 
ukb_FCIadjBMI <- ukb_FCI[,c("IID","FCI_adjBMI")] 



ukb_MS_ISI <- fread("MS_ISI_MAGICRI_UKBB.dosage.gz")
map.mono <- fread("MS_ISI_MAGICRI_UKBB.stats.gz")
head(map.mono)
map.mono$alleleB_dosage_code <- as.character(map.mono$alleleB_dosage_code)
map.mono <- merge(map.mono,MS_ISI,by="rsid")
head(map.mono)
map.mono$weight_UA <- ifelse(map.mono$ea == map.mono$alleleB_dosage_code, map.mono$effect_noBMI, -map.mono$effect_noBMI)
map.mono$weight_Adj <- ifelse(map.mono$ea == map.mono$alleleB_dosage_code, map.mono$effect_adjBMI, -map.mono$effect_adjBMI)
snps <- ukb_MS_ISI
id <- ukb_MS_ISI$IID
rsid <- map.mono$rsid
for (i in rsid){
  ukb_MS_ISI[[i]]<- ukb_MS_ISI[[i]]*map.mono$weight_UA[map.mono$rsid %in% i]
}

ukb_MS_ISI$MS_ISI_noBMI <- NA
ukb_MS_ISI$MS_ISI_noBMI <- rowSums(ukb_MS_ISI[,c(2:8)]) 
ukb_MS_ISInoBMI <- ukb_MS_ISI[,c("IID","MS_ISI_noBMI")] 
ukb_MS_ISInoBMI <- merge(ukb_MS_ISInoBMI,snps,by="IID")

id <- ukb_MS_ISI$IID
rsid <- map.mono$rsid

ukb_MS_ISI <- fread("MS_ISI_MAGICRI_UKBB.dosage.gz")

for (i in rsid){
  ukb_MS_ISI[[i]]<- ukb_MS_ISI[[i]]*map.mono$weight_Adj[map.mono$rsid %in% i]
}

ukb_MS_ISI$MS_ISI_adjBMI <- NA
ukb_MS_ISI$MS_ISI_adjBMI <- rowSums(ukb_MS_ISI[,c(2:8)]) 
ukb_MS_ISIadjBMI <- ukb_MS_ISI[,c("IID","MS_ISI_adjBMI")] 

UKBB <- merge(ukb_FCIadjBMI,ukb_FCInoBMI,by = "IID")

UKBB <- merge(UKBB,ukb_MS_ISIadjBMI,by = "IID")
UKBB <- merge(UKBB,ukb_MS_ISInoBMI,by = "IID")
head(UKBB)
write.table(UKBB,"IR_scores_UKBB.txt",sep = "\t",quote = F,row.names = F)





