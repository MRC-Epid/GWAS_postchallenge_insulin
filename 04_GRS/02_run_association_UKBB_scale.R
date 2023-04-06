setwd("/analysis_dirs/02_run_association/")

score <- read.table("../01_generate_scores/IR_scores_UKBB.txt",header=T,stringsAsFactors=F)
library(data.table)

#prepared phentype file clean
pheno <- fread("UKB_data.txt",header = T,stringsAsFactors = F)

geno_pheno <- merge(score,pheno,by="IID")
#remove relatedness
geno_pheno <- geno_pheno[geno_pheno$Ex_for_white_noREL %in% 0,]

#adjust age sex and 1st 10 PCs
adj.linreg <- paste0("age_at_first_check+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")

feat <- c("list of phenotypes")
snps <- c("list of scores")
snps
names(score)
sum_geno <- NULL
 for (i in feat){
  
  for (k in snps) {
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
  file_write <- paste0("./UKBB_results/IRscores_EUR_unrelated_",i,".txt")
file_write
  write.table(tmp,file_write,row.names = F,quote = F,sep = "\t")
  
}

#logistic for T2D
feat <- "T2D_PrevInc_Feb2021"
for (i in feat){
  
  for (k in snps) {
    geno <- paste(i, " ~ ", k," + ",adj.linreg, sep = "")
    print(paste0("running reg: ",geno))
    sum_geno[[k]] <- summary(glm(as.formula(geno),data = geno_pheno,family="binomial"))
  }
  
  tmp <- lapply(sum_geno,function(x) x$coefficients[2,])
  
  tmp2 <- names(tmp)
  tmp <- do.call(rbind,tmp)
  tmp <- as.data.frame(tmp)
  names(tmp) <- c("Estimate","SE","tval","pval")
  tmp$score <- tmp2
  file_write <- paste0("./UKBB_results/IRscores_EUR_unrelated_",i,".txt")
  file_write
  write.table(tmp,file_write,row.names = F,quote = F,sep = "\t")
  
}



