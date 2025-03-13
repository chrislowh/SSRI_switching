#!/bin/bash
#SBATCH --job-name=step_2_regenie_parallel
#SBATCH --partition=cpu
#SBATCH --mem=50G
#SBATCH --array=1-22
#SBATCH --time=16:00:00
#SBATCH --mail-user=chris.lowh@kcl.ac.uk
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --nodes=1

#Please feel free to go through my directories and have a look at files, take a copy of my SNPlists from the directories below.
#You will need to submit this script as sbatch after running

# Source bash script and conda executable, otherwise conda environment installed in own shell cannot be activated
source /users/k21157612/.bashrc
source /users/k21157612/miniconda_temp/etc/profile.d/conda.sh

conda activate regenie
#Which can be installed using: 
#conda create -n regenie -c conda-forge -c bioconda regenie

#The pred files are generated from step 1. 
#The included individuals are just taken from the qc folder, took the first 2 cols from the fam, 
#And converted to text format with headers.


PFILE="/scratch/prj/ukbiobank/recovered/ukb82087/imputed/"
VARS="/scratch/prj/ukbiobank/recovered/Edinburgh_Data/usr/chrislo/ukb_gwas/data/imputed_snplist"
OUT="/scratch/prj/ukbiobank/chrislo/regenie/step_2/"
PRED="/scratch/prj/ukbiobank/chrislo/regenie/step_1/"
PHENOS="/scratch/prj/ukbiobank/recovered/Edinburgh_Data/usr/chrislo/ukb_gwas/data/phenotype/pgx/switch/"
PHENO_TEMP="/users/k21157612/ukb_rap/temp/ad_switch/input/pheno/"
COVS="/users/k21157612/ukb_rap/temp/ad_switch/input/cov/"
TOOLS="/scratch/prj/ukbiobank/recovered/Edinburgh_Data/usr/dale_handley/tools/"
LOG="/users/k21157612/ukb_rap/temp/ad_switch/job_logs/"

ssri_all="regenie \
--step 2 \
--pgen $PFILE"ukb82087_imp_chr"$SLURM_ARRAY_TASK_ID"_MAF1_INFO4_v1" \
--extract $VARS"/SNPs_imputed_qc_chr"$SLURM_ARRAY_TASK_ID".snplist" \
--keep $TOOLS"ukb82087_post_QC_individuals.txt" \
--phenoFile $PHENO_TEMP"id_switch_SSRI_all.txt" \
--covarFile $COVS"cov_switch_SSRI_all.txt" \
--phenoCol "cc" \
--bt \
--bsize 200 \
--threads 32 \
--catCovarList "Sex,Centre" \
--covarColList "index_date,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" \
--maxCatLevels 200 \
--firth \
--approx \
--pred $PRED"switch_SSRI_all_step1_pred.list" \
--out $OUT"switch_SSRI_all_GWAS_chr"$SLURM_ARRAY_TASK_ID > $LOG"SSRI_switch_all_regenie_firth_step_2""

ssri_onemdd="regenie \
--step 2 \
--pgen $PFILE"ukb82087_imp_chr"$SLURM_ARRAY_TASK_ID"_MAF1_INFO4_v1" \
--extract $VARS"/SNPs_imputed_qc_chr"$SLURM_ARRAY_TASK_ID".snplist" \
--keep $TOOLS"ukb82087_post_QC_individuals.txt" \
--phenoFile $PHENOS"id_switch_SSRI_onemdd.txt" \
--covarFile $COVS"cov_switch_SSRI_onemdd.txt" \
--phenoCol "cc" \
--bt \
--bsize 200 \
--threads 32 \
--catCovarList "Sex,Centre" \
--covarColList "index_date,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" \
--maxCatLevels 200 \
--firth \
--approx \
--pred $PRED"switch_SSRI_onemdd_step1_pred.list" \
--out $OUT"switch_SSRI_onemdd_GWAS_chr"$SLURM_ARRAY_TASK_ID > $LOG"SSRI_switch_onemdd_regenie_firth_step_2""

ssri_twomdd="regenie \
--step 2 \
--pgen $PFILE"ukb82087_imp_chr"$SLURM_ARRAY_TASK_ID"_MAF1_INFO4_v1" \
--extract $VARS"/SNPs_imputed_qc_chr"$SLURM_ARRAY_TASK_ID".snplist" \
--keep $TOOLS"ukb82087_post_QC_individuals.txt" \
--phenoFile $PHENOS"id_switch_SSRI_twomdd.txt" \
--covarFile $COVS"cov_switch_SSRI_twomdd.txt" \
--phenoCol "cc" \
--bt \
--bsize 200 \
--threads 32 \
--catCovarList "Sex,Centre" \
--covarColList "index_date,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" \
--maxCatLevels 200 \
--firth \
--approx \
--pred $PRED"switch_SSRI_twomdd_step1_pred.list" \
--out $OUT"switch_SSRI_twomdd_GWAS_chr"$SLURM_ARRAY_TASK_ID > $LOG"SSRI_switch_twomdd_regenie_firth_step_2""


#eval $ssri_all
#eval $ssri_onemdd
#eval $ssri_twomdd