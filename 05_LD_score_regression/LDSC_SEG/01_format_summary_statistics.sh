#MarkerName	Allele1	Allele2	Freq1	Effect	StdErr	Pvalue	TotalSampleSize	chr	pos
cd /path_to_analysis_dir/LDSC_SEGG/ldsc


module load python-2.7.13-gcc-5.4.0-yubmrmn
module load anaconda/3.2019-10
#conda env create --force --file environment.yml python=2.7.13
source activate ldsc


#run for IFC. --sumstats file has been formatted to meet input requirements. autosomes only
python munge_sumstats.py \
--sumstats ../IFC_MA_adjBMI_autosomes.txt \
--chunksize 500000 \
--snp rsid \
--a1 Allele1 \
--a2 Allele2 \
--p P-value \
--frq Freq1 \
--N-col TotalSampleSizeAll \
--out /path_to_analysis_dir/LDSC_SEGG/IFC_adjBMI_autosomes


#run for ISI
python munge_sumstats.py \
--sumstats ../ISI_MA_adjBMI_autosomes.txt \
--chunksize 500000 \
--snp rsid \
--a1 Allele1 \
--a2 Allele2 \
--p Pvalue \
--frq Freq1 \
--N-col N \
--out /path_to_analysis_dir/LDSC_SEGG/ISI_adjBMI_autosomes

