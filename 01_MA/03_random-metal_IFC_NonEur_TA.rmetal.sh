###  metal script for random effects meta-analysis of all Insulin Fold Change cohorts.
## Trans-ancestry - EUR, EAS, HIS AMR 
## NON-EUR - EAS, HIS AMR

## output columns: randomeffects "\tEffectARE\tStdErrARE\tPvalueARE\ttausq\tStdErrMRE\tPvalueMRE
# ARE - dlstats - DerSimonian-Laird random effects model
# MRE - phi

#random-metal : /path_to/resources/random-metal/executables/metal 

								##############################
								#       IFC adjBMI ALL       #
								##############################

SCHEME STDERR
SEPARATOR TAB
MARKERLABEL MarkerName
ALLELELABELS Allele1 Allele2
PVALUELABEL P-value
EFFECTLABEL Effect
STDERR StdErr

CUSTOMVARIABLE TotalSampleSizeAll
LABEL TotalSampleSizeAll as TotalSampleSize
FREQLABEL Freq1
AVERAGEFREQ ON
MINMAXFREQ ON
ADDFILTER Freq1 > 0.005
ADDFILTER Freq1 < 0.995
### EUR
PROCESS /path_to_results/METAANALYSIS_IFC_adjBMI_EUR1.tbl
### HISAMR
PROCESS /path_to_results/METAANALYSIS_IFC_adjBMI_AMR1.tbl
### EAS
PROCESS /path_to_results/METAANALYSIS_IFC_adjBMI_EAS1.tbl

OUTFILE METAANALYSIS_IFC_adjBMI_ALL_RE .tbl
#run RANDOM EFFECTS
ANALYZE RANDOM
CLEAR


									### genomic control of MA ###
SCHEME STDERR
SEPARATOR TAB
GENOMICCONTROL ON
AVERAGEFREQ ON
MINMAXFREQ ON

# column names
MARKER MarkerName
ALLELE Allele1 Allele2
EFFECT EffectARE
STDERR StdErrARE
PVALUE PvalueARE
WEIGHT TotalSampleSize
FREQLABEL Freq1
WEIGHTLABEL TotalSampleSize
 CUSTOMVARIABLE TotalSampleSize
 LABEL TotalSampleSize as TotalSampleSizeAll

PROCESS METAANALYSIS_IFC_adjBMI_ALL_RE1.tbl

OUTFILE METAANALYSIS_IFC_adjBMI_ALL_RE_DGC .tbl
ANALYZE
CLEAR


								##############################
								#       IFC noBMI ALL        #
								##############################

SCHEME STDERR
SEPARATOR TAB
MARKERLABEL MarkerName
ALLELELABELS Allele1 Allele2
PVALUELABEL P-value
EFFECTLABEL Effect
STDERR StdErr

CUSTOMVARIABLE TotalSampleSizeAll
LABEL TotalSampleSizeAll as TotalSampleSize
FREQLABEL Freq1
AVERAGEFREQ ON
MINMAXFREQ ON
ADDFILTER Freq1 > 0.005
ADDFILTER Freq1 < 0.995

### EUR
PROCESS /path_to_results/TEmeta/METAANALYSIS_IFC_noBMI_EUR1.tbl
### HISAMR
PROCESS /path_to_results/TEmeta/METAANALYSIS_IFC_noBMI_AMR1.tbl
### EAS
PROCESS /path_to_results/TEmeta/METAANALYSIS_IFC_noBMI_EAS1.tbl

OUTFILE METAANALYSIS_IFC_noBMI_ALL_RE .tbl
ANALYZE RANDOM
CLEAR

									### genomic control of MA ###
### DGC
SCHEME STDERR
SEPARATOR TAB
GENOMICCONTROL ON
AVERAGEFREQ ON
MINMAXFREQ ON
# column names
MARKER MarkerName
ALLELE Allele1 Allele2
EFFECT EffectARE
STDERR StdErrARE
PVALUE PvalueARE
WEIGHT TotalSampleSize
FREQLABEL Freq1
WEIGHTLABEL TotalSampleSizeAll
 CUSTOMVARIABLE TotalSampleSize
 LABEL TotalSampleSize as TotalSampleSizeAll

PROCESS METAANALYSIS_IFC_noBMI_ALL_RE1.tbl

OUTFILE METAANALYSIS_IFC_noBMI_ALL_RE_DGC .tbl
ANALYZE
CLEAR




								##################################
								#       IFC adjBMI NONEUR        #
								##################################

SCHEME STDERR
SEPARATOR TAB
MARKERLABEL MarkerName
ALLELELABELS Allele1 Allele2
PVALUELABEL P-value
EFFECTLABEL Effect
STDERR StdErr

CUSTOMVARIABLE TotalSampleSizeAll
LABEL TotalSampleSizeAll as TotalSampleSize
FREQLABEL Freq1
AVERAGEFREQ ON
MINMAXFREQ ON
ADDFILTER Freq1 > 0.005
ADDFILTER Freq1 < 0.995
### HISAMR
PROCESS /path_to_results/METAANALYSIS_IFC_adjBMI_AMR1.tbl
### EAS
PROCESS /path_to_results/METAANALYSIS_IFC_adjBMI_EAS1.tbl

OUTFILE METAANALYSIS_IFC_adjBMI_NONEUR_RE .tbl
#run RANDOM EFFECTS
ANALYZE RANDOM
CLEAR


									### now genomic control of MA ###
SCHEME STDERR
SEPARATOR TAB
GENOMICCONTROL ON
AVERAGEFREQ ON
MINMAXFREQ ON

# column names
MARKER MarkerName
ALLELE Allele1 Allele2
EFFECT EffectARE
STDERR StdErrARE
PVALUE PvalueARE
WEIGHT TotalSampleSize
FREQLABEL Freq1
WEIGHTLABEL TotalSampleSizeAll
 CUSTOMVARIABLE TotalSampleSize
 LABEL TotalSampleSize as TotalSampleSizeAll

PROCESS METAANALYSIS_IFC_adjBMI_NONEUR_RE1.tbl

OUTFILE METAANALYSIS_IFC_adjBMI_NONEUR_RE_DGC .tbl
ANALYZE
CLEAR


								#################################
								#       IFC noBMI NONEUR        #
								#################################

SCHEME STDERR
SEPARATOR TAB
MARKERLABEL MarkerName
ALLELELABELS Allele1 Allele2
PVALUELABEL P-value
EFFECTLABEL Effect
STDERR StdErr

CUSTOMVARIABLE TotalSampleSizeAll
LABEL TotalSampleSizeAll as TotalSampleSize
FREQLABEL Freq1
AVERAGEFREQ ON
MINMAXFREQ ON
ADDFILTER Freq1 > 0.005
ADDFILTER Freq1 < 0.995

### HISAMR
PROCESS /path_to_results/TEmeta/METAANALYSIS_IFC_noBMI_AMR1.tbl
### EAS
PROCESS /path_to_results/TEmeta/METAANALYSIS_IFC_noBMI_EAS1.tbl

OUTFILE METAANALYSIS_IFC_noBMI_NONEUR_RE .tbl
ANALYZE RANDOM
CLEAR

									### genomic control of MA ###
SCHEME STDERR
SEPARATOR TAB
GENOMICCONTROL ON
AVERAGEFREQ ON
MINMAXFREQ ON
# column names
MARKER MarkerName
ALLELE Allele1 Allele2
EFFECT EffectARE
STDERR StdErrARE
PVALUE PvalueARE
WEIGHT TotalSampleSizeAll
FREQLABEL Freq1
WEIGHTLABEL TotalSampleSizeAll
 CUSTOMVARIABLE TotalSampleSize
 LABEL TotalSampleSize as TotalSampleSizeAll

PROCESS METAANALYSIS_IFC_noBMI_NONEUR_RE1.tbl

OUTFILE METAANALYSIS_IFC_noBMI_NONEUR_RE_DGC .tbl
ANALYZE
CLEAR




