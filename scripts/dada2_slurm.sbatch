#!/bin/bash

#SBATCH --partition=ieg_lm, large_mem  ### Partition (like a queue in PBS)
#SBATCH --job-name=dada2_duolun  ### Job Name
#SBATCH -o /condo/ieg/jianshu/log/jarray.%j.%N.out  ### File in which to store job output
#SBATCH -e /condo/ieg/jianshu/log/jarray.%j.%N.err  ### File in which to store job error
#SBATCH --time=15-00:00:00          ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                   ### Node count required for the job
#SBATCH --ntasks=1                  ### Nuber of tasks to be launched per Node
#SBATCH --cpus-per-task=24          ### Number of threads per task (OMP threads)
#SBATCH --mem=200G                   ### memory for each job, this is large dataset
#SBATCH --mail-type=FAIL            ### When to send mail
#SBATCH --mail-user=jianshuzhao@yahoo.com. ### mail to send
#SBATCH --get-user-env              ### Import your user environment setup
#SBATCH --requeue                   ### On failure, requeue for another try
#SBATCH --verbose                   ### Increase informational messages

module purge
### activate conda environment with installed R, install miniconda3 first
source ~/.bashrc
conda init bash
conda activate base
which Rscript

git clone https://github.com/jianshu93/dada2_wrapper
cd dada2_wrapper/scripts
time Rscript ./dada2_wrapper.r --input_dir=../demo_input --output_dir=output_paired --pool --threads 24