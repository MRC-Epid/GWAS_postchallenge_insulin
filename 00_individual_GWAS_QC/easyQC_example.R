args = commandArgs(trailingOnly=TRUE)

setwd("/path_to_projecteasyQC/Scripts/METSIM")

#load EasyQC R package
if(!require(EasyQC)){
install.packages("EasyQC")
library(EasyQC)
}

#run easy QC for study
EasyQC('EasyQC_V17_METSIM_FCI_AW.ecf')
