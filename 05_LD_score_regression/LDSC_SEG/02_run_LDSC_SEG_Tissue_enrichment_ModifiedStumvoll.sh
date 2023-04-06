#!/bin/bash
#!
#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! Name of the job:
#SBATCH -J ldsc_seg_insulin

#! Which project should be charged:
#SBATCH -A MRC-EPID-SL0-CPU

#! How many whole nodes should be allocated?
#SBATCH --nodes=1

#! Specify required run time
#SBATCH --time=8:00:00

#SBATCH --ntasks=3

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL

#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! Do not change:
#SBATCH -p epid

#SBATCH --output=slurm-%x-%j.out
#! ############################################################
#! Modify the settings below to specify the application's environment, location 
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

module list
#! Insert additional module load commands after this line if needed:

module load python-2.7.13-gcc-5.4.0-yubmrmn
module list

#conda needs to be python 2 too
module load miniconda/2

cd /path_to_analysis_dir/LDSC_SEGG/ldsc

source activate ldsc

#source /path_to_resources/ldsc/python_ve/bin/activate


#Insulin fold change
#note copied the ./ldcts files to this directory and modified file paths based on what downloaded so they contain the full mapping needed to run
# modified for example: Multi_tissue_chromatin_1000Gv3_ldscores/ENTEX.1.    to     /path_to_analysis_dir/Multi_tissue_chromatin_1000Gv3_ldscores/ENTEX.1.
python ldsc.py \
	--h2-cts /path_to_analysis_dir/LDSC_SEGG/IFC_adjBMI_autosomes.sumstats.gz \
	--ref-ld-chr /path_to_resources/ldsc/1000G_EUR_Phase3_baseline/baseline. \
	--ref-ld-chr-cts Multi_tissue_chromatin.ldcts \
	--w-ld-chr /path_to_resources/ldsc/weights_hm3_no_hla/weights. \
  --out /path_to_analysis_dir/LDSC_SEGG/IFC_adjBMI_multi_tissue_chromatin 

#did same for the gene_exp.ldcts as outlined above

python ldsc.py \
--h2-cts /path_to_analysis_dir/LDSC_SEGG/IFC_adjBMI_autosomes.sumstats.gz \
--ref-ld-chr /path_to_resources/ldsc/1000G_EUR_Phase3_baseline/baseline. \
--ref-ld-chr-cts Multi_tissue_gene_expr.ldcts \
--w-ld-chr /path_to_resources/ldsc/weights_hm3_no_hla/weights. \
--out /path_to_analysis_dir/LDSC_SEGG/IFC_adjBMI_multi_tissue_gene_exp



#same for modified stumvoll ISI
python ldsc.py \
	--h2-cts /path_to_analysis_dir/LDSC_SEGG/ISI_adjBMI_autosomes.sumstats.gz \
	--ref-ld-chr /path_to_resources/ldsc/1000G_EUR_Phase3_baseline/baseline. \
	--ref-ld-chr-cts Multi_tissue_chromatin.ldcts \
	--w-ld-chr /path_to_resources/ldsc/weights_hm3_no_hla/weights. \
	--out /path_to_analysis_dir/LDSC_SEGG/ISI_adjBMI_multi_tissue_chromatin 



python ldsc.py \
--h2-cts /path_to_analysis_dir/LDSC_SEGG/ISI_adjBMI_autosomes.sumstats.gz \
--ref-ld-chr /path_to_resources/ldsc/1000G_EUR_Phase3_baseline/baseline. \
--ref-ld-chr-cts Multi_tissue_gene_expr.ldcts \
--w-ld-chr /path_to_resources/ldsc/weights_hm3_no_hla/weights. \
--out /path_to_analysis_dir/LDSC_SEGG/ISI_adjBMI_multi_tissue_gene_exp

#conda deactivate






