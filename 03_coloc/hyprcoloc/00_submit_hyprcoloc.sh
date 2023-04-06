#!/bin/bash
#!


#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! Name of the job:
#SBATCH -J hypr_coloc_harmonised

#! Which project should be charged:
#SBATCH -A MRC-EPID-SL0-CPU


#! How many whole nodes should be allocated?
#SBATCH --nodes=1

#SBATCH --ntasks=5



#SBATCH --output=slurm-%x-%j.out


#! Specify required run time
#SBATCH --time=12:00:00

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
cd /analysis_dir/
tail -n +2 ./input/MAGIC_ISI_IFC_regions.txt | while IFS='\t' read line;
do
echo $line 
grep "$line" ./input/MAGIC_ISI_IFC_regions.txt > ./input/VOI.txt
REG="$(awk '{print $1}' ./input/VOI.txt)"
echo "region is $REG"
CHR="$(awk '{print $9}' ./input/VOI.txt)"
echo "chromosome is $CHR"
POSS="$(awk '{print $12}' ./input/VOI.txt)"
echo "start is $POSS"
POSE="$(awk '{print $13}' ./input/VOI.txt)"
echo "end is $POSE"



 
echo "Run HyprColoc for region ${REG} on chromosome ${CHR} between ${POSS} and ${POSE}"

## run the script
01_run_hyprcoloc.R $REG $CHR $POSS $POSE
done 

