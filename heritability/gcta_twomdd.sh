#!/bin/bash -l
#SBATCH --job-name=gcta
#SBATCH --output="/users/k21157612/ukb_rap/temp/ad_switch/job_logs/%j_gcta_SSRI_switch_twomdd"
#SBATCH --partition=cpu
#SBATCH --mem=50G
#SBATCH --time=16:00:00
#SBATCH --mail-user=chris.lowh@kcl.ac.uk
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --nodes=1
#SBATCH --array=1-10

#paths
OUT="/scratch/prj/ukbiobank/chrislo/gcta"
PFILE="/scratch/prj/ukbiobank/recovered/ukb82087/imputed/"
SOFTWARE="/scratch/users/k21157612/recovered/software/gcta/gcta-1.94.1-linux-kernel-3-x86_64"
BFILE="/scratch/prj/ukbiobank/chrislo/genetic_qc"
SGDP_QC="/scratch/prj/ukbiobank/recovered/ukb82087/qc/"
PHENO="/scratch/prj/ukbiobank/recovered/Edinburgh_Data/usr/chrislo/ukb_gwas/data/phenotype/pgx/switch"
SNPS="/scratch/prj/ukbiobank/recovered/Edinburgh_Data/usr/dale_handley/tools/imputation_scores"
LOG="/scratch/prj/ukbiobank/recovered/Edinburgh_Data/usr/chrislo/ad_dosage/job_logs"
COV="/users/k21157612/ukb_rap/temp/ad_switch/input/cov"
HM3="/scratch/users/k21157612/recovered/genome_ref/hm3"


# create the genetic relationship matrix using 10 subsets of the genetic data in plink binary format:
# bfile
# binary file using QC parameters from SGDP (in 82087 folder)
# --keep: phenotype file
# --extract: variants passing standard SGDP QC (should be done already in post-qc file, but can check in log)

grm="$SOFTWARE"/gcta-1.94.1" \
--make-grm-part 10 $SLURM_ARRAY_TASK_ID \
--bfile $BFILE"/ukb82087_post_qc" \
--keep $PHENO"/id_switch_SSRI_twomdd_noheader.txt" \
--extract $SGDP_QC"ukb82087_post_qc_snp_list_relatives_included.txt" \
--out $OUT"/SSRI_twomdd_gcta"" 

eval $grm

cd $OUT

# concatenate the different subsets:
cat SSRI_twomdd_gcta.part_10_*.grm.id > SSRI_twomdd_gcta.grm.id
cat SSRI_twomdd_gcta.part_10_*.grm.bin > SSRI_twomdd_gcta.grm.bin
cat SSRI_twomdd_gcta.part_10_*.grm.N.bin > SSRI_twomdd_gcta.grm.N.bin

# remove the subsets:
rm SSRI_twomdd_gcta.part_10*

cd /

# From here, run interactively (arrays unncessary)
# adjust GRM for incomplete tagging of causal SNPs:
adj="$SOFTWARE"/gcta-1.94.1" \
--grm $OUT"/SSRI_twomdd_gcta" \
--grm-adj 0 \
--make-grm \
--out $OUT"/SSRI_twomdd_gcta.adjusted""

# Remove related subjects using grm-cutoff 0.05:
remove_rel="$SOFTWARE"/gcta-1.94.1" \
--grm $OUT"/SSRI_twomdd_gcta.adjusted" \
--grm-cutoff 0.05 \
--make-grm \
--out $OUT"/SSRI_twomdd_gcta.adjusted.unrel""

# Run GREML-SC:
greml="$SOFTWARE"/gcta-1.94.1" \
--grm $OUT"/SSRI_twomdd_gcta.adjusted.unrel" \
--reml \
--reml-no-constrain \
--pheno $PHENO"/id_switch_SSRI_twomdd_noheader.txt" \
--covar $COV"/cov_switch_SSRI_twomdd_cat_noheader.txt" \
--qcovar $COV"/cov_switch_SSRI_twomdd_num_noheader.txt" \
--out $OUT"/SSRI_twomdd_gcta.adjust.unrel.GREML_SC""

#eval $adj
#eval $remove_rel
#eval $greml