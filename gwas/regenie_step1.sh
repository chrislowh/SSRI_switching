#!/bin/bash -l
#SBATCH --job-name=step_1_regenie_parallel
#SBATCH --partition=cpu
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --mail-user=chris.lowh@kcl.ac.uk
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --nodes=4

#Please feel free to go through my directories and have a look at files, take a copy of my SNPlists from the directories below.
#You will need to submit this script as sbatch after running

# Source bash script and conda executable, otherwise conda environment installed in own shell cannot be activated
source /users/k21157612/.bashrc
source /users/k21157612/miniconda_temp/etc/profile.d/conda.sh

conda activate regenie
#Which can be installed using: 
#conda create -n regenie -c conda-forge -c bioconda regenie

#The phenotypes can be submitted as a single file, but if you use the --strict flag, this removes all individuals who are missing
#Any data at any phenotype, which is undesirable. 
#However, if this flag isn't used, imputation is performed, which might be untrustworthy given the specific missingness/nature of 
#Your data. 


#Please use my bfiles that I've generated by extracting the pruned SNPs.
#If you run on the full genotyped data, even with autosome only flags, regenie will throw an error, as 
#it does not understand plink's coding of sex chromosomes.
#This bfile contains ~500k pruned SNPs (window 1000, size 100, r2 = 0.9, as per regenie documentation) in 
#the 458k individuals of European ancestry who pass all pre-defined QC steps. At some stage, I'll update this
#To run my own QCed list. 

OUT="/scratch/prj/ukbiobank/chrislo/regenie/step_1/"
BFILE="/scratch/prj/ukbiobank/chrislo/genetic_regenie/"
PHENOS="/scratch/prj/ukbiobank/recovered/Edinburgh_Data/usr/chrislo/ukb_gwas/data/phenotype/pgx/switch/"
PHENO_TEMP="/users/k21157612/ukb_rap/temp/ad_switch/input/pheno/"
COVS="/users/k21157612/ukb_rap/temp/ad_switch/input/cov/"
VARS="/scratch/prj/ukbiobank/recovered/Edinburgh_Data/usr/dale_handley/idk/variants/"
LOG="/users/k21157612/ukb_rap/temp/ad_switch/job_logs/"
SGDP_QC="/scratch/prj/ukbiobank/recovered/ukb82087/qc/"


# Flags:
# --bed: Genotype (.bed/.bim/.fam) files from UKB82087 (QC'ed), then extracted with independent SNPs pruned for LD (as the one below), and only in QCed individuals (done in 82087) 
# --extract: independent SNPs pruned for LD (from Dale Handley), as recommended in original publication for REGENIE
# --keep: post qc individual list
## bed files should have already filtered for both --extract / --keep flag, but just put them here to ensure, in case upstream genotype file has been updated. 
# --phenoFile: phenotype file
# --covarFile: covariates files
# bt - use if trait is standard binary (IE: 0/1/NA)
# No need for lowmem as only 3 phenotypes. If more than 10, use --lowmem. 
# Strict - doesn't allow imputation

# SSRI
ssri_all="regenie \
--step 1 \
--bed $BFILE"ukb82087_genotyped_independent_regenie" \
--extract $VARS"pruned_independent_regenie_SNPs_step1_all.txt" \
--keep $SGDP_QC"ukb82087_post_qc_id_list_relatives_included.fam" \
--phenoFile $PHENO_TEMP"id_switch_SSRI_all.txt" \
--covarFile $COVS"cov_switch_SSRI_all.txt" \
--phenoCol "cc" \
--bt \
--bsize 1000 \
--threads 120 \
--catCovarList "Sex,Centre" \
--covarColList "index_date,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" \
--maxCatLevels 200 \
--strict \
--lowmem \
--lowmem-prefix $OUT"regenie_tmp_SSRI_all" \
--out $OUT"switch_SSRI_all_step1" > $LOG"switch_SSRI_all_regenie_step_1""

# SSRI
ssri_onemdd="regenie \
--step 1 \
--bed $BFILE"ukb82087_genotyped_independent_regenie" \
--extract $VARS"pruned_independent_regenie_SNPs_step1_all.txt" \
--keep $SGDP_QC"ukb82087_post_qc_id_list_relatives_included.fam" \
--phenoFile $PHENOS"id_switch_SSRI_onemdd.txt" \
--covarFile $COVS"cov_switch_SSRI_onemdd.txt" \
--phenoCol "cc" \
--bt \
--bsize 1000 \
--threads 120 \
--catCovarList "Sex,Centre" \
--covarColList "index_date,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" \
--maxCatLevels 200 \
--strict \
--lowmem \
--lowmem-prefix $OUT"regenie_tmp_SSRI_onemdd" \
--out $OUT"switch_SSRI_onemdd_step1" > $LOG"switch_SSRI_onemdd_regenie_step_1""

ssri_twomdd="regenie \
--step 1 \
--bed $BFILE"ukb82087_genotyped_independent_regenie" \
--extract $VARS"pruned_independent_regenie_SNPs_step1_all.txt" \
--keep $SGDP_QC"ukb82087_post_qc_id_list_relatives_included.fam" \
--phenoFile $PHENOS"id_switch_SSRI_twomdd.txt" \
--covarFile $COVS"cov_switch_SSRI_twomdd.txt" \
--phenoCol "cc" \
--bt \
--bsize 1000 \
--threads 120 \
--catCovarList "Sex,Centre" \
--covarColList "index_date,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" \
--maxCatLevels 200 \
--strict \
--lowmem \
--lowmem-prefix $OUT"regenie_tmp_SSRI_twomdd" \
--out $OUT"switch_SSRI_twomdd_step1" > $LOG"switch_SSRI_twomdd_regenie_step_1""

#eval $ssri_all
#eval $ssri_onemdd
#eval $ssri_twomdd