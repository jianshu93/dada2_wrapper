#!/bin/bash

#SBATCH --partition=ieg_64g,ieg_128g,ieg_lm  ### Partition (like a queue in PBS)
#SBATCH --job-name=Haibei_wp_megahit  ### Job Name
#SBATCH -o ~/jarray.%j.%N.out  ### File in which to store job output
#SBATCH -e ~/jarray.%j.%N.err  ### File in which to store job error
#SBATCH --time=08-00:00:00          ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                   ### Node count required for the job
#SBATCH --ntasks=1                  ### Nuber of tasks to be launched per Node
#SBATCH --cpus-per-task=24          ### Number of threads per task (OMP threads)
#SBATCH --mem=120G                   ### memory for each job
#SBATCH --mail-type=FAIL            ### When to send mail
#SBATCH --mail-user=jianshuzhao@yahoo.com. ### mail to send
#SBATCH --get-user-env              ### Import your user environment setup
#SBATCH --requeue                   ### On failure, requeue for another try
#SBATCH --verbose                   ### Increase informational messages

###
module purge

### activate conda environment, install miniconda3 first

source ~/.bash_profile
conda init bash
conda activate base
which Rscript


cd /storage/home/hcoda1/4/jzhao399/p-ktk3-0/rich_project_bio-konstantinidis/dada2/scripts
Rscript ./dada2_cli.r --input_dir=input --output_dir=output --pool
Rscript ./dada2_cli.r --input_dir=input --output_dir=output --pool