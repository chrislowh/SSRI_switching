#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --time=30:00:00
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=chris.lowh@kcl.ac.uk
#SBATCH --job-name=genopred

################################################
# Please kindly:
# 1. change the paths for batch job options as above

# 2. If you are installing your own miniconda3 (as per Ollie's standard instructions)
    # on CREATE, the default conda will be called when executing bash scripts
    # therefore, your codes would only work in command line
    # if you want to submit as bash script, the errors would look like as you have never installed GenoPred conda environment
    # So please source your own conda executable (as below)

# 3. 1st time using GenoPred
    # these steps can be done using Ollie's instruction on the website
        # create conda environment to host GenoPred
        # create ANOTHER conda environment to host packages used to run steps in GenoPred
            # the dependencies are included in the .yaml file in github by Ollie

    # if you run into problems with conda channels (which somehow happens for me on CREATE)
        # please try disabling channel priorities, should install without conflicts
        # re-set back channel priorities to strict afterwards
            # avoid dependencies conflicts in the future if you want to installing many packages into conda with multiple channels
            # last time I tried to install GenoPred, it was specified with multiple channels, but the packages were not in conflict
###############################################

### --- Activate conda environment --- ###
# Source bash script and conda executable
# otherwise conda environment installed in own shell cannot be activated
# To : your home directory should have ./bashrc and you will need to have a look where you have installed conda executable
source /users/k21157612/.bashrc
source /scratch/users/k21157612/miniconda3/etc/profile.d/conda.sh


### --- 1st Time running GenoPred --- ###
# create conda environment (genopred) to host snakemake running pipeline
cd /scratch/prj/ukbiobank/usr/chrislo/multiPGS
conda env create -f GenoPred/pipeline/envs/pipeline.yaml

# activate environment to host snakemake (genopred)
conda activate genopred
# Navigate to working directory
cd /scratch/prj/ukbiobank/usr/chrislo/multiPGS/GenoPred/pipeline

# create another conda environment for the pipeline under .yaml file in github

### --- Useful if running into problems with channel conflicts --- ###
# some conflicts with channels when installing
# disable channel priority, and set it back to strict after creating the conda environment (to avoid future conflicts)
conda config --set channel_priority disabled

# download software and dependencies
# I am copying this line from Ollie's configurations
# But if you run into problems with mamba, you might want to remove that flag
# It is not absolutely necessary to run GenoPred with mamba, but it will make the installation quicker
snakemake -j 1 --use-conda --conda-frontend mamba install_r_packages

# re-enable channel priority
conda config --set channel_priority strict

# Note:
# Comment the above lines after running it for first time
# the conda environment is created already so you don't need to do it again


### --- Actual running GenoPred --- ###
# one line command
# --profile slurm: allows to specify resources and number of jobs to be run, very useful on HPC if you want to speed up
snakemake --configfile=config_ukb.yaml target_pgs --profile slurm --latency-wait 150 --rerun-incomplete
