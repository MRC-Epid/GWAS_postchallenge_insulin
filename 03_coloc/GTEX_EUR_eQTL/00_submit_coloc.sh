#!/bin/bash
#!

#! Slurm script to run jobs on UoC HPC system

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! Name of the job:
#SBATCH -J coloc_gtex_EUR

#! Which project should be charged:
#SBATCH -A MRC-EPID-SL0-CPU

#! How many whole nodes should be allocated?
#SBATCH --nodes=1

#SBATCH --ntasks=4

#! Specify required run time
#SBATCH --time=10:00:00

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL

#SBATCH --output=slurm-%x-%j.out

#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
#SBATCH --no-requeue

#! Partition to run on:
#SBATCH -p epid

#! ############################################################
#! Modify the settings below to specify the application's environment, location
#! and launch method:
#! uncomment this to exclude certain nodes e.g. 1-25 SBATCH --exclude=cpu-d-[1-25]
#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment



#! Insert additional module load commands after this line if needed:
module load r-3.6.0-gcc-5.4.0-bzuuksv
## list all modules loaded 
module list


cd /path_to_analysis_dir/

#loop over file containing the key information needed to run the coloc for each locus
#contains locus label, chromosome, position of lead variant in b37 and b38 (needed for GTEx subsetting), rsid of lead variant at locus
cat loci_coloc.txt| while IFS='\t' read line;
do
echo $line 

echo $line > VOI.txt
locus="$(awk '{print $1}' VOI.txt)"
echo "locus is $locus"
CHR="$(awk '{print $3}' VOI.txt)"
echo "chromosome is $CHR"
RSID="$(awk '{print $2}' VOI.txt)"
echo "start is $RSID"
POS37="$(awk '{print $4}' VOI.txt)"
echo "end is $POS38"
POS38="$(awk '{print $5}' VOI.txt)"
echo "end is $POS38"
echo "Run Coloc for locus ${locus} on chromosome ${CHR} at b37 ${POS37} and b38 ${POS38}"

#eQTL files all tissues +/- 200kb window 
echo "running all tissues of interest coloc with 200kb window"
Rscript --no-save --no-restore --verbose 01_colocalisation_eQTL_EUR.R ${locus} ${RSID} ${CHR} ${POS37} ${POS38}

done
