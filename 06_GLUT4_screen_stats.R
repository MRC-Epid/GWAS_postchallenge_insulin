

#t-test to identify sig diff from NT control 
setwd("/analysis_dir/")                                                                
library(data.table)
tmp_combined_res<- fread("combined_data_allreplicates_processesed.txt")
tmp_combined_res <- as.data.frame(tmp_combined_res)


condition <- c("surface_total_glut4","surface_glut4_nuc","total_glut4_nuc")
tmp_combined_res<- tmp_combined_res[tmp_combined_res$Insulin_Dose %in% 0.5,]
i <- condition
for (i in condition){
  #default parameters, pairwise = F & alternative = "two.sided" - so standard two-sided t-test
  #first argument is column containing measure of interest e.g. surface GLUT4/ nuclei
  #second argment - how values are grouped, in this case "name" denotes KD condition
pHA <- as.data.frame(pairwise.t.test(tmp_combined_res[,i],tmp_combined_res$name,p.adjust.method = "none",paired = F,alternative="two.sided")$p.value)
#this generates lower triangle of a matrix of all comparisons
write.table(pHA,paste0("./t_test/paired_t_test_compareall_",i,".txt"),quote = T,sep = "\t",row.names = F)

#subset to only compared to NT control - primary interest
pHA <- as.data.frame(pHA)
pHA2 <- pHA #store copy to get lower part of matrix in a moment - provides results in lower triangle of a  matrix
pHA$names <- row.names(pHA)
pHA <- pHA[,c("names","NT ctrl")]
pHA2 <- t(pHA2)
pHA2 <- pHA2[,"NT ctrl"]
pHA2 <- as.data.frame(pHA2)
pHA2$names <- row.names(pHA2)

names(pHA2)[1]<- "NT ctrl"
pHA2 <- pHA2[,c(2,1)]
pHA <- pHA[!is.na(pHA$`NT ctrl`),]
pHA2 <- pHA2[!is.na(pHA2$`NT ctrl`),]

pHA_res <-rbind(pHA,pHA2) 
pHA_res <-rbind(pHA,pHA2)
pHA_res$P_FDR <- p.adjust(pHA_res$`NT ctrl`,method = "fdr")
write.table(pHA_res,paste0("./t_test/paired_t_test_compareNT_",i,".txt"),quote = T,sep = "\t",row.names = F)

}
