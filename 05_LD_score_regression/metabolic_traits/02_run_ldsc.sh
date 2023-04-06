
module load python-2.7.13-gcc-5.4.0-yubmrmn
#conda needs to be python 2 too
module load miniconda/2

# file directory /input_dir
source activate ldsc

## run LD score regression - following guidelines at https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation

## LD score reference was downloaded from here: https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/eur_w_ld_chr.tar.bz2

#run ldscore regression

##note does pairwise corr with the first trait listed and all subsequent traits therefore need to ensure the 1st trait is the primary trait of interest

## MODIFIED STUMVOLL ISI adjBMI
ldsc.py \
--rg /input_dir/formatted_summ_stats/ISI_adjBMI_autosomes.sumstats.gz,/input_dir/summ_stats/PASS_BMI1.sumstats,/input_dir/summ_stats/PASS_Coronary_Artery_Disease.sumstats,/input_dir/summ_stats/PASS_Type_1_Diabetes.sumstats,/input_dir/summ_stats/PASS_Type_2_Diabetes.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_AlanineAminotransferase.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Albumin.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_AlkalinePhosphatase.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_ApolipoproteinA.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_ApolipoproteinB.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_AspartateAminotransferase.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Calcium.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Cholesterol.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_CreactiveProtein.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Creatinine.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_CystatinC.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_DirectBilirubin.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_GammaGlutamyltransferase.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Glucose.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_HbA1c.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_HDLcholesterol.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_IGF1.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_LDLdirect.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_LipoproteinA.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Phosphate.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_SHBG.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_TotalBilirubin.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_TotalProtein.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Triglycerides.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Urate.sumstats,/input_dir/summ_stats/UKB_460K.body_BMIz.sumstats,/input_dir/summ_stats/UKB_460K.body_WHRadjBMIz.sumstats,/input_dir/summ_stats/UKB_460K.bp_DIASTOLICadjMEDz.sumstats,/input_dir/summ_stats/UKB_460K.bp_SYSTOLICadjMEDz.sumstats,/input_dir/summ_stats/UKB_460K.disease_CARDIOVASCULAR.sumstats,/input_dir/summ_stats/UKB_460K.disease_DIABETES_ANY_DIAGNOSED.sumstats,/input_dir/summ_stats/UKB_460K.disease_ENDOCRINE_DIABETES.sumstats,/input_dir/summ_stats/UKB_460K.disease_HYPERTENSION_DIAGNOSED.sumstats,/input_dir/summ_stats/UKB_460K.disease_T2D.sumstats,/input_dir/formatted_summ_stats/IFC_adjBMI_autosomes.sumstats.gz,/input_dir/formatted_summ_stats/FI_adjBMI.sumstats.gz,/input_dir/formatted_summ_stats/Glu120_adjBMI.sumstats.gz,/input_dir/formatted_summ_stats/ISI_adjBMI_autosomes.sumstats.gz,/input_dir/formatted_summ_stats/T2D_MVP.sumstats.gz,/input_dir/formatted_summ_stats/WHRadjBMI.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out /input_dir/results/LD_score_cardiometabolic_glycaemic_results_ISIadjBMI



## INSULIN FOLD CHANGE adjBMI
ldsc.py \
--rg /input_dir/formatted_summ_stats/IFC_adjBMI_autosomes.sumstats.gz,/input_dir/summ_stats/PASS_BMI1.sumstats,/input_dir/summ_stats/PASS_Coronary_Artery_Disease.sumstats,/input_dir/summ_stats/PASS_Type_1_Diabetes.sumstats,/input_dir/summ_stats/PASS_Type_2_Diabetes.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_AlanineAminotransferase.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Albumin.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_AlkalinePhosphatase.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_ApolipoproteinA.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_ApolipoproteinB.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_AspartateAminotransferase.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Calcium.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Cholesterol.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_CreactiveProtein.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Creatinine.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_CystatinC.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_DirectBilirubin.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_GammaGlutamyltransferase.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Glucose.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_HbA1c.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_HDLcholesterol.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_IGF1.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_LDLdirect.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_LipoproteinA.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Phosphate.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_SHBG.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_TotalBilirubin.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_TotalProtein.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Triglycerides.sumstats,/input_dir/summ_stats/UKB_460K.biochemistry_Urate.sumstats,/input_dir/summ_stats/UKB_460K.body_BMIz.sumstats,/input_dir/summ_stats/UKB_460K.body_WHRadjBMIz.sumstats,/input_dir/summ_stats/UKB_460K.bp_DIASTOLICadjMEDz.sumstats,/input_dir/summ_stats/UKB_460K.bp_SYSTOLICadjMEDz.sumstats,/input_dir/summ_stats/UKB_460K.disease_CARDIOVASCULAR.sumstats,/input_dir/summ_stats/UKB_460K.disease_DIABETES_ANY_DIAGNOSED.sumstats,/input_dir/summ_stats/UKB_460K.disease_ENDOCRINE_DIABETES.sumstats,/input_dir/summ_stats/UKB_460K.disease_HYPERTENSION_DIAGNOSED.sumstats,/input_dir/summ_stats/UKB_460K.disease_T2D.sumstats,/input_dir/formatted_summ_stats/IFC_adjBMI_autosomes.sumstats.gz,/input_dir/formatted_summ_stats/FI_adjBMI.sumstats.gz,/input_dir/formatted_summ_stats/Glu120_adjBMI.sumstats.gz,/input_dir/formatted_summ_stats/ISI_adjBMI_autosomes.sumstats.gz,/input_dir/formatted_summ_stats/T2D_MVP.sumstats.gz,/input_dir/formatted_summ_stats/WHRadjBMI.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out /input_dir/results/LD_score_cardiometabolic_glycaemic_results_IFCadjBMI



#exit conda environment
conda deactivate

