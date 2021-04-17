#!/bin/bash

#SBATCH --partition=ieg_lm  ### Partition (like a queue in PBS)
#SBATCH --job-name=dada2_duolun  ### Job Name
#SBATCH -o /condo/ieg/jianshu/log/jarray.%j.%N.out  ### File in which to store job output
#SBATCH -e /condo/ieg/jianshu/log/jarray.%j.%N.err  ### File in which to store job error
#SBATCH --time=15-00:00:00          ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                   ### Node count required for the job
#SBATCH --ntasks=1                  ### Nuber of tasks to be launched per Node
#SBATCH --cpus-per-task=16          ### Number of threads per task (OMP threads)
#SBATCH --mem=1200G                   ### memory for each job
#SBATCH --mail-type=FAIL            ### When to send mail
#SBATCH --mail-user=jianshuzhao@yahoo.com. ### mail to send
#SBATCH --get-user-env              ### Import your user environment setup
#SBATCH --requeue                   ### On failure, requeue for another try
#SBATCH --verbose                   ### Increase informational messages

module purge
### activate conda environment, install miniconda3 first
source ~/.bash_profile
conda init bash
conda activate base
which Rscript

cd /condo/ieg/jianshu/app/dada2_wrapper/scripts
time Rscript ./dada2_wrapper.r --input_dir=/condo/ieg/jianshu/app/duolun --output_dir=duolun_output --pool --threads 16