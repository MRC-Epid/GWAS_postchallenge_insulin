#generate the scores and save in file for analysis input
## load the effects for weighting 
setwd("/analysis_dir/01_generate_scores")
library(data.table)

FCI <- fread("../input/FCI_stats.txt",header=T,stringsAsFactors=F,sep="\t")
MS_ISI <- fread("../input/MS_ISI_stats.txt",header=T,stringsAsFactors=F,sep = "\t")

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
##Fenland OMICs


fen_FCI <- fread("FCI_MAGICRI_Fenland-OMICS_hrcuk10k.dosage")
map.mono <- fread("FCI_MAGICRI_Fenland-OMICS_hrcuk10k.dosage_code")
head(map.mono)
map.mono$dosage_coded_allele <- as.character(map.mono$dosage_coded_allele)
map.mono <- merge(map.mono,FCI,by="rsid")
head(map.mono)
map.mono$weight_UA <- ifelse(map.mono$ea==map.mono$dosage_coded_allele, map.mono$effect_noBMI, -map.mono$effect_noBMI)
map.mono$weight_Adj <- ifelse(map.mono$ea==map.mono$dosage_coded_allele, map.mono$effect_adjBMI, -map.mono$effect_adjBMI)

id <- fen_FCI$id
rsid <- map.mono$rsid
for (i in rsid){
  fen_FCI[[i]]<- fen_FCI[[i]]*map.mono$weight_UA[map.mono$rsid %in% i]
}

fen_FCI$FCI_noBMI <- NA
fen_FCI$FCI_noBMI <- rowSums(fen_FCI[,c(2:6)]) 
fen_FCInoBMI <- fen_FCI[,c("id","FCI_noBMI")] 


id <- fen_FCI$id
rsid <- map.mono$rsid

fen_FCI <- fread("FCI_MAGICRI_Fenland-OMICS_hrcuk10k.dosage")

for (i in rsid){
  fen_FCI[[i]]<- fen_FCI[[i]]*map.mono$weight_Adj[map.mono$rsid %in% i]
}

fen_FCI$FCI_adjBMI <- NA
fen_FCI$FCI_adjBMI <- rowSums(fen_FCI[,c(2:6)]) 
fen_FCIadjBMI <- fen_FCI[,c("id","FCI_adjBMI")] 



fen_MS_ISI <- fread("MS_ISI_MAGICRI_Fenland-OMICS_hrcuk10k.dosage")
map.mono <- fread("MS_ISI_MAGICRI_Fenland-OMICS_hrcuk10k.dosage_code")
head(map.mono)
map.mono$dosage_coded_allele <- as.character(map.mono$dosage_coded_allele)
map.mono <- merge(map.mono,MS_ISI,by="rsid")
head(map.mono)
map.mono$weight_UA <- ifelse(map.mono$ea == map.mono$dosage_coded_allele, map.mono$effect_noBMI, -map.mono$effect_noBMI)
map.mono$weight_Adj <- ifelse(map.mono$ea == map.mono$dosage_coded_allele, map.mono$effect_adjBMI, -map.mono$effect_adjBMI)

id <- fen_MS_ISI$id
rsid <- map.mono$rsid
for (i in rsid){
  fen_MS_ISI[[i]]<- fen_MS_ISI[[i]]*map.mono$weight_UA[map.mono$rsid %in% i]
}

fen_MS_ISI$MS_ISI_noBMI <- NA
fen_MS_ISI$MS_ISI_noBMI <- rowSums(fen_MS_ISI[,c(2:9)]) 
fen_MS_ISInoBMI <- fen_MS_ISI[,c("id","MS_ISI_noBMI")] 


id <- fen_MS_ISI$id
rsid <- map.mono$rsid

fen_MS_ISI <- fread("MS_ISI_MAGICRI_Fenland-OMICS_hrcuk10k.dosage")

for (i in rsid){
  fen_MS_ISI[[i]]<- fen_MS_ISI[[i]]*map.mono$weight_Adj[map.mono$rsid %in% i]
}

fen_MS_ISI$MS_ISI_adjBMI <- NA
fen_MS_ISI$MS_ISI_adjBMI <- rowSums(fen_MS_ISI[,c(2:9)]) 
fen_MS_ISIadjBMI <- fen_MS_ISI[,c("id","MS_ISI_adjBMI")] 



Fenland_OMICs <- merge(fen_FCIadjBMI,fen_FCInoBMI,by = "id")
Fenland_OMICs <- merge(Fenland_OMICs,fen_MS_ISIadjBMI,by = "id")
Fenland_OMICs <- merge(Fenland_OMICs,fen_MS_ISInoBMI,by = "id")
head(Fenland_OMICs)
write.table(Fenland_OMICs,"IR_scores_Fenland_OMICS.txt",sep = "\t",quote = F,row.names = F)





