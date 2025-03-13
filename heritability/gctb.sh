#!/bin/bash -l
#SBATCH --job-name=gctb
#SBATCH --output="/users/k21157612/ukb_rap/temp/ad_switch/job_logs/%j_gctb_ssri_switch"
#SBATCH --partition=cpu
#SBATCH --mem=50G
#SBATCH --time=42:00:00
#SBATCH --mail-user=chris.lowh@kcl.ac.uk
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --nodes=1

#paths
OUT="/scratch/prj/ukbiobank/chrislo/gctb/"
BFILE="/scratch/prj/ukbiobank/chrislo/genetic_qc"
SGDP_QC="/scratch/prj/ukbiobank/recovered/ukb82087/qc/"
PHENO="/scratch/prj/ukbiobank/recovered/Edinburgh_Data/usr/chrislo/ukb_gwas/data/phenotype/pgx/switch"
PHENO_TEMP="/users/k21157612/ukb_rap/temp/ad_switch/input/pheno"
SNPS="/scratch/prj/ukbiobank/recovered/Edinburgh_Data/usr/dale_handley/tools/imputation_scores"
LOG="/users/k21157612/ukb_rap/temp/ad_switch/job_logs/"
COV="/users/k21157612/ukb_rap/temp/ad_switch/input/cov/"

module add gctb/2.04.3-gcc-13.2.0

ssri_all="gctb \
--bfile $BFILE"/ukb82087_post_qc" \
--pheno $PHENO_TEMP"/id_switch_SSRI_all_noheader.txt" \
--extract $SGDP_QC"ukb82087_post_qc_snp_list_relatives_included.txt" \
--covar $COV"cov_switch_SSRI_all_gctb_noheader.txt" \
--bayes S \
--pi 0.05 \
--hsq 0.5 \
--S 0 \
--chain-length 21000 \
--burn-in 1000 \
--out-freq 500 \
--seed 123 \
--out $OUT"gctb_switch_SSRI_all" > $LOG"gctb_switch_SSRI_all""

ssri_onemdd="gctb \
--bfile $BFILE"/ukb82087_post_qc" \
--pheno $PHENO"/id_switch_SSRI_onemdd_noheader.txt" \
--extract $SGDP_QC"ukb82087_post_qc_snp_list_relatives_included.txt" \
--covar $COV"cov_switch_SSRI_onemdd_gctb_noheader.txt" \
--bayes S \
--pi 0.05 \
--hsq 0.5 \
--S 0 \
--chain-length 21000 \
--burn-in 1000 \
--out-freq 500 \
--seed 123 \
--out $OUT"gctb_switch_SSRI_onemdd" > $LOG"gctb_switch_SSRI_onemdd""

ssri_twomdd="gctb \
--bfile $BFILE"/ukb82087_post_qc" \
--pheno $PHENO"/id_switch_SSRI_twomdd_noheader.txt" \
--extract $SGDP_QC"ukb82087_post_qc_snp_list_relatives_included.txt" \
--covar $COV"cov_switch_SSRI_twomdd_gctb_noheader.txt" \
--bayes S \
--pi 0.05 \
--hsq 0.5 \
--S 0 \
--chain-length 21000 \
--burn-in 1000 \
--out-freq 500 \
--seed 123 \
--out $OUT"gctb_switch_SSRI_twomdd" > $LOG"gctb_switch_SSRI_twomdd""

eval $ssri_all
#eval $ssri_onemdd
#eval $ssri_twomdd
