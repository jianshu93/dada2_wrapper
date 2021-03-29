DADA2 R workflow for profiling 16S sequence reads
======


* This repo contains workflows for analysing exact sequence variants from 16S sequencing reads using the DADA2 algorithm (https://github.com/benjjneb/dada2).

* Primer removal

Before you start, make sure that you remove all the primers in the forwared and 
reverse reads. There is a bash script based on cutadapt software in the demo_input 
directory (run in this directory only), run it will remoe primers for the forward and reverse reads. The primer is now set to universal 16S V3-V4 primer set 515F and 806R. Change the -a and -A option in the cut_primer.sh to your primer used during the amplification experiment. The demo_input reads contains no primer

* IMPORTANT
this pipeline is extensively tested under conda R version 4.0.2 on both MacOS and Linux (Ubuntu 18.0.4 and RHEL 7). I strongly suggest you reinstall a new R 4.0.2 from scratch but not update R in you conda. INSTALL A COMPLETELY NEW R 4.0.2! This will save you lot of trouble (updating R packages from an old version to a new one is annoying) 

# DADA2 can now be run with this command (*_R1.fastq, *_R2.fastq or gzipped formatshould be in the input directory):
```
git clone https://github.com/jianshu93/dada2_wrapper
cd dada2_wrapper/scripts
./dada2_cli.r --input_dir=../demo_input --output_dir=output --pool
```


DADA2 R workflow for profiling 16S sequence reads consists of these two scripts:
--------------------------------------------------------------------------------

#### 1. dada2_cli.r - This wrapper script parses arguments from the command line and passes them to and invokes the DADA2 workflow (dada2_16S_paired-end.Rmd).

#### 2. dada2_16S_paired-end.Rmd - DADA2 workflow for resolving sequence variants from 16S amplicon reads. In addition to OTU picking and taxonomy assignment, it will also produce a QC folder with visuals and read counts (before/after). 
  NOTE: This workflow requires paired-end 16S rDNA gene amplicon reads.
  NOTE: This workflow will analyze each sample independently to be able to process any sample size.  The option to pool samples for resolving sequence variants will be added.
  NOTE: Currently, the wrapper script dada2_cli.r and the workflow dada2_16S_paired-end.Rmd need to be in the same directory (and not in the input directory). By default the database is at the parent directory of the script dada2_cli.r, which is: ../reference_dbs_16S

It is suggested that the pool mode is more sensitive than sample-based mode for detecting rare ASVs (http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/)


#### Installation of DADA2 from bioconductor


#### Installation of FastTree if not available in your environment, please using the OpenMP version, which is much faster and have identical results with regular version
http://www.microbesonline.org/fasttree/#Install

#### I strongly suggest that you use miniconda3 to manage your environment. Both on clusters and locally. See how to install miniconda3 here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

#### Installation of dada2 stable version. Go to R commandline interface first and then (R version >= 3.6.3 is preferred, it is tested for both r3.6.3 and r4.0.2). For MacOS, there are some problems with R conda version 3.6.3 and thus I suggest you use the newest R conda version 4.0.2. It works nicely:
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager",repos='http://cran.us.r-project.org')
library(BiocManager)

requiredPackages = c('RCurl','ggplot2','GenomicRanges','SummarizedExperiment',
                    'BiocParallel','Rsamtools','dada2','msa','phangorn','gridExtra','knitr')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) BiocManager::install(p,update = FALSE)
}
```
#### Time needed and platform dependent questions
for the testing dataset, 6 samples, a few million reads per sample, it takes ~20 minutes on a 4-core (8 threads) Macbook Pro. For a 72 sample dataset with more than 10 million reads per sample, it takes about 5 hours on a 24 threads Linux system on a cluster node.

#### Running DADA2 on Phoenix cluster/TORQUE based or SLURM based cluster


```


# # If logged in on a cluster (by default, this workflow will use all the processes/threads you have on your cluster or the processes you requested), install miniconda3 first and then activate conda environment and install R packages as mentioned above and then clone the repo. DADA2 can then be run on Phoenix using the Torque queue manager with the example script in the repo:
qsub dada2_torque.pbs
## or on a slurm queue manager with the sample
sbatch dada2_wrapper.sbatch


## for interactive mode on TORQUE and SLURM system:

qsub -I -l nodes=1:ppn=24 -l walltime=10:00:00 -l mem=600gb -q inferno -A GT-ktk3-CODA20 -o /storage/home/hcoda1/4/jzhao399/p-ktk3-0/rich_project_bio-konstantinidis/scripts/log/${PBS_JOBNAME}_${PBS_JOBID}.out -e /storage/home/hcoda1/4/jzhao399/p-ktk3-0/rich_project_bio-konstantinidis/scripts/log/${PBS_JOBNAME}_${PBS_JOBID}.err -N dada2

srun --partition=ieg_plus --nodes=1 --ntasks-per-node=1 --cpus-per-task=24 --time=24:00:00 --mem=120G --pty bash -i

## then run as it is a regular terminal

git clone https://github.com/jianshu93/dada2_wrapper
cd dada2_wrapper/scripts
./dada2_cli.r --input_dir=../demo_input --output=output --pool

# NOTE: Prior to running make sure to install DADA2 R package dependencies. I have added a few lines to install all packages in the dada2_cli.r script after you have successfully install miniconda3 and activate base environment.

```
you will need pandoc and texlive to generate pdf files.

```
conda install pandoc

```

if you do not have textlive installed in your system, the programm will stop and only a .tex file will be generated. No worries! You can transform this tex file to pdf file here: https://cloudconvert.com/tex-to-pdf


For you own dataset just change the input directory to where your pair-end fastq file is. Do remember to remove primers first before running this pipeline. Check the script in the demo_input directory on how to do that. 

Please contact jianshu.zhao@gatech.edu if you have any questions.

