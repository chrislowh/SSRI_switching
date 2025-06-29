# Antidepressant Switching as a Proxy Phenotype for Drug Nonresponse: Investigating Clinical, Demographic, and Genetic Characteristics

### This repository includes the R scripts to create switching phenotype and perform genetic analyses in the study: 
### Lo et al, 2025. Antidepressant Switching as a Proxy Phenotype for Drug Nonresponse: Investigating Clinical, Demographic, and Genetic Characteristics. Biol Psychiatry Glob Open Sci. 
#### available at doi: [https://doi.org/10.1016/j.bpsgos.2025.100502](https://doi.org/10.1016/j.bpsgos.2025.100502)

<br>

## Description of contents:

### gwas

	  # regenie_step1.sh: Bash script to run step 1 of REGENIE (https://rgcgithub.github.io/regenie/) with LD-pruned independent SNPs on HPC.

	  # regenie_step2.sh: Bash script to run step 2 of REGENIE (https://rgcgithub.github.io/regenie/) with imputed SNPs on HPC, parallelized by chromosome (1-22).

### heritability

	  # gcta_all.sh: Bash script to run GCTA on full sample for SSRI switchers and non-switchers.

 	  # gcta_onemdd.sh: Bash script to run GCTA on sample for SSRI switchers and non-switchers (with >= 1 depression diagnosis record in primary care).

    # gcta_twomdd.sh: Bash script to run GCTA on sample for SSRI switchers and non-switchers (with >= 2 depression diagnoses record in primary care).

    # gctb.sh: Bash script to run GCTB on all three samples described in the manuscript.

### pgs

    # genopred_pgs.sh: Bash script to initiate GenoPred pipeline (creating conda environments), and running polygenic scores with the pipeline on HPC.

### phenotype

    # switch_phenotype.R: R script to create switchers and non-switchers by running switch_Lo2024().

For phenotyping, we are currently developing the T-Rx package to simplify phenotyping processes and data cleaning of prescription records. Please contact Chris Lo (chris.lowh@kcl.ac.uk) for details.

