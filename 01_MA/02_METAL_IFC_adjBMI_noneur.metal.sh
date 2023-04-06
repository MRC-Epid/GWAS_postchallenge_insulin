###  metal script for meta-analysis of Insulin Fold Change cohorts.


#metal/2011-03-25  

######################## Filters and Metanalysis parameters ####################################################
SCHEME   STDERR
SEPARATOR   TAB
MARKERLABEL   cpaid
ALLELELABELS  EFFECT_ALLELE NON_EFFECT_ALLELE
PVALUELABEL   P_VAL
EFFECTLABEL   BETA
STDERR   SE
GENOMICCONTROL ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N
FREQLABEL EAF
AVERAGEFREQ ON
MINMAXFREQ ON
ADDFILTER EAF > 0.005
ADDFILTER EAF < 0.995

############################# read in the files to be meta-analysed #####################################################

### HTN-IR
GENOMICCONTROL ON
PROCESS /path_to_projectCLEANED_Studies/CLEAN_FCI/CLEANED.IFC_adjBMI_HTN-IR_AMR_rvtests_MINIMAC_20200226_JT.gz
##MACAD 
GENOMICCONTROL ON
PROCESS /path_to_projectCLEANED_Studies/CLEAN_FCI/CLEANED.IFC_adjBMI_MACAD_AMR_rvtests_MINIMAC_20200226_JT.gz

############################ OUTPUT ##########################################################
## name the outfile; saved in the current working directory  ###

OUTFILE METAANALYSIS_IFC_adjBMI_AMR .tbl


### Analyse heterogeneity - calculates the heterogeneity 
ANALYZE HETEROGENEITY

### clear metal

CLEAR


SCHEME   STDERR
SEPARATOR   TAB
MARKERLABEL   cpaid
ALLELELABELS  EFFECT_ALLELE NON_EFFECT_ALLELE
PVALUELABEL   P_VAL
EFFECTLABEL   BETA
STDERR   SE
GENOMICCONTROL ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N
FREQLABEL EAF
AVERAGEFREQ ON
MINMAXFREQ ON
ADDFILTER EAF > 0.005
ADDFILTER EAF < 0.995

## TAICHI
GENOMICCONTROL ON
PROCESS /path_to_projectCLEANED_Studies/CLEAN_FCI/CLEANED.IFC_adjBMI_TAICHI_EAS_rvtests_MINIMAC_20200324_JT.txt



OUTFILE METAANALYSIS_IFC_adjBMI_EAS .tbl


### Analyse heterogeneity - calculates the heterogeneity 
ANALYZE HETEROGENEITY

### clear metal

CLEAR



##### Double genomic control
######################## Filters and Metanalysis parameters ####################################################

### DGC
SCHEME STDERR
SEPARATOR TAB
GENOMICCONTROL ON
AVERAGEFREQ ON
MINMAXFREQ ON



# column names
MARKER MarkerName
ALLELE Allele1 Allele2
EFFECT Effect
STDERR StdErr
PVALUE P-value
WEIGHT TotalSampleSize
FREQLABEL Freq1
WEIGHTLABEL TotalSampleSize
 CUSTOMVARIABLE TotalSampleSizeAll
 LABEL TotalSampleSizeAll as TotalSampleSize


############################# read in the files to be meta-analysed #####################################################


### IFC_adjBMI_NON-EUR
PROCESS METAANALYSIS_IFC_adjBMI_AMR1.tbl

############################ OUTPUT ##########################################################
## name the outfile; saved in the current working directory  ###

OUTFILE METAANALYSIS_IFC_adjBMI_AMR_DGC .tbl
ANALYZE
### exit metal 
QUIT


### DGC
SCHEME STDERR
SEPARATOR TAB
GENOMICCONTROL ON
AVERAGEFREQ ON
MINMAXFREQ ON



# column names
MARKER MarkerName
ALLELE Allele1 Allele2
EFFECT Effect
STDERR StdErr
PVALUE P-value
WEIGHT TotalSampleSize
FREQLABEL Freq1
WEIGHTLABEL TotalSampleSize
 CUSTOMVARIABLE TotalSampleSizeAll
 LABEL TotalSampleSizeAll as TotalSampleSize


############################# read in the files to be meta-analysed #####################################################


### IFC_adjBMI_NON-EUR
PROCESS METAANALYSIS_IFC_adjBMI_EAS1.tbl

############################ OUTPUT ##########################################################
## name the outfile; saved in the current working directory  ###

OUTFILE METAANALYSIS_IFC_adjBMI_EAS_DGC .tbl
ANALYZE
### exit metal 
QUIT


