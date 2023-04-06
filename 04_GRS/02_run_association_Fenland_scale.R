setwd("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Alice/IR_scores/02_run_association/")


score <- read.table("../01_generate_scores/IR_scores_Fenland_OMICS.txt",header=T,stringsAsFactors=F)
library("data.table")
library(haven)

#read in pre-prepared phenotype file
pheno <- read_dta("./Fenland_R8a_cleaned.dta")
PCs <- read_dta("./Fenland_PCs.dta")
geno_pheno <- merge(score,pheno,by="id")

#merge with genetic PC and remove ancestry outlier and related
geno_pheno<- merge(geno_pheno,PCs,by="id")
geno_pheno <- geno_pheno[!geno_pheno$relative_excl %in% 1,]
geno_pheno <- geno_pheno[!geno_pheno$ancestry_outlier %in% 1,]
adj.linreg <- paste0("age+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")

snps_list <- c(names(score)[2:7])

sum_geno <- NULL
 for (i in feat){
  
  for (k in snps_list) {
    geno <- paste("scale(",i,")", " ~ ", k," + ",adj.linreg, sep = "")
    print(paste0("running reg: ",geno))
    sum_geno[[k]] <- summary(lm(as.formula(geno),data = geno_pheno))
  }
  
  tmp <- lapply(sum_geno,function(x) x$coefficients[2,])
  
  tmp2 <- names(tmp)
  tmp <- do.call(rbind,tmp)
  tmp <- as.data.frame(tmp)
  names(tmp) <- c("Estimate","SE","tval","pval")
  tmp$score <- tmp2
  file_write <- paste0("./Fenland_results/IRscores_EUR_unrelated_",i,".txt")
file_write
 write.table(tmp,file_write,row.names = F,quote = F,sep = "\t")
  
}


