#!/bin/bash
#!


#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! Name of the job:
#SBATCH -J coloc

#! Which project should be charged:
#SBATCH -A MRC-EPID-SL0-CPU


#! How many whole nodes should be allocated?
# In your case you should leave this at 1
#SBATCH --nodes=1

#SBATCH --ntasks=3



#SBATCH --output=slurm-%x-%j.out


#! Specify required run time
#SBATCH --time=15:00:00

## restrict the number of nodes to be used (maximum of 30)


#! What types of email messages do you wish to receive?
#SBATCH --mail-type=BEGIN,FAIL,END

#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! Do not change:
#SBATCH -p epid



#! ############################################################
#! Modify the settings below to specify the application's environment, location 
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load r-3.6.0-gcc-5.4.0-bzuuksv
module load gcc/5
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

#! Insert additional module load commands after this line if needed:

## assign directories used for the analysis
DIR=/analysis_dir/

cd ${DIR}
echo "beginning coloc Insulin Fold Change"

#extract positions etc to loop over in the script
tail -n +2 ./input/MAGIC_IFC_adjBMI_regions.txt | while IFS='\t' read line;
do
echo $line 
grep "$line" ./input/MAGIC_IFC_adjBMI_regions.txt > ./input/VOI.txt
REG="$(awk '{print $1}' ./input/VOI.txt)"
echo "region is $REG"
CHR="$(awk '{print $9}' ./input/VOI.txt)" #lead variant chr and pos
echo "chromosome is $CHR"
POS="$(awk '{print $10}' ./input/VOI.txt)"
echo "lead position is $POS"

echo "Run Coloc for locus ${REG}"

## run the script
Rscript --no-save --no-restore --verbose 02_run_coloc_MVP_IFC.R $CHR $POS ;
done 


