*STATA script to perform Wakefield finemapping at IFC adjBMI EUR loci. 
*use conditional stats for PPP1R3B to fine map primary and secondary signals
cd "/path_to_analysis_dir/finemap_wakefield"

insheet using "../GCTA_COJO/IFC_adjBMI_EUR_C2CD4A.ma"
gen z_stat = b/se
gen V = se^2
local W = 0.04
gen r = `W'/(V+`W')
gen double ABF = 1/(exp(-r*z_stat*z_stat/2)/sqrt(1-r))
egen double total_ABF  = total(ABF)
gen double post_prob = ABF/total_ABF
gsort - ABF
gen cum_post_prob = sum(post_prob)
*get PP of the lead SNP at this locus
list post_prob if snp == "rs7167878"

export delimited IFC_adjBMI_EUR_C2CD4A_finemap.txt,replace
clear

insheet using "../GCTA_COJO/IFC_adjBMI_EUR_MTN1RB.ma"
gen z_stat = b/se
gen V = se^2
local W = 0.04
gen r = `W'/(V+`W')
gen double ABF = 1/(exp(-r*z_stat*z_stat/2)/sqrt(1-r))
egen double total_ABF  = total(ABF)
gen double post_prob = ABF/total_ABF
gsort - ABF
gen cum_post_prob = sum(post_prob)
list post_prob if snp == "rs1387153"

export delimited IFC_adjBMI_EUR_MTN1RB_finemap.txt,replace

clear
insheet using "../GCTA_COJO/IFC_adjBMI_EUR_PPP1R3B.ma"
gen z_stat = b/se
gen V = se^2
local W = 0.04
gen r = `W'/(V+`W')
gen double ABF = 1/(exp(-r*z_stat*z_stat/2)/sqrt(1-r))
egen double total_ABF  = total(ABF)
gen double post_prob = ABF/total_ABF
gsort - ABF
gen cum_post_prob = sum(post_prob)
list post_prob if snp == "rs7012814"

export delimited IFC_adjBMI_EUR_PPP1R3B_finemap.txt,replace

clear
insheet using "../GCTA_COJO/secondary_signals/conditional_Fenland-OMICS_IFC_adjBMI_PPP1R3B_chr8.cma.cojo"
encode bc, generate(bc2)
encode bc_se, generate(bc_se2)
gen z_stat = bc2/bc2_se
gen V = bc2_se^2
local W = 0.04
gen r = `W'/(V+`W')
gen double ABF = 1/(exp(-r*z_stat*z_stat/2)/sqrt(1-r))
egen double total_ABF  = total(ABF)
gen double post_prob = ABF/total_ABF
gsort - ABF
gen cum_post_prob = sum(post_prob)
list post_prob if snp == "rs7012814"

export delimited IFC_adjBMI_EUR_PPP1R3B_conditional_finemap.txt,replace

clear
insheet using "../GCTA_COJO/IFC_adjBMI_EUR_SLC2A4.ma"
gen z_stat = b/se
gen V = se^2
local W = 0.04
gen r = `W'/(V+`W')
gen double ABF = 1/(exp(-r*z_stat*z_stat/2)/sqrt(1-r))
egen double total_ABF  = total(ABF)
gen double post_prob = ABF/total_ABF
gsort - ABF
gen cum_post_prob = sum(post_prob)
list post_prob if snp == "rs117643180"

export delimited IFC_adjBMI_EUR_SLC2A4_finemap.txt,replace
clear

import delim using "../GCTA_COJO/secondary_signals/conditional_Fenland-OMICS_IFC_adjBMI_PPP1R3B_cond_primary_chr8.cma.cojo", numericcols(10,11,12)
gen z_stat = bc/bc_se
gen V = bc_se^2
local W = 0.04
gen r = `W'/(V+`W')
gen double ABF = 1/(exp(-r*z_stat*z_stat/2)/sqrt(1-r))
egen double total_ABF  = total(ABF)
gen double post_prob = ABF/total_ABF
gsort - ABF
gen cum_post_prob = sum(post_prob)

list post_prob if snp == "rs7012814"

export delimited IFC_adjBMI_EUR_PPP1R3B_condition_on_primary_finemap.txt,replace
clear
import delim using "../GCTA_COJO/secondary_signals/conditional_Fenland-OMICS_IFC_adjBMI_PPP1R3B_cond_secondary_chr8.cma.cojo", numericcols(10,11,12)
gen z_stat = bc/bc_se
gen V = bc_se^2
local W = 0.04
gen r = `W'/(V+`W')
gen double ABF = 1/(exp(-r*z_stat*z_stat/2)/sqrt(1-r))
egen double total_ABF  = total(ABF)
gen double post_prob = ABF/total_ABF
gsort - ABF
gen cum_post_prob = sum(post_prob)

list post_prob if snp == "rs7012814"

export delimited IFC_adjBMI_EUR_PPP1R3B_condition_on_secondary_finemap.txt,replace
