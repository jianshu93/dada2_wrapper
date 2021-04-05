DADA2 R workflow for profiling 16S sequence reads
======


* This repo contains workflows for analysing exact sequence variants from 16S sequencing reads using the DADA2 algorithm (https://github.com/benjjneb/dada2).

YOU MUST HAVE R INSTALLED BEFORE DOING ANYTHING

* Primer removal

Before you start, I suggest you remove all the primers in the forwared and 
reverse reads. There is a bash script based on cutadapt software in the demo_input 
directory (run in this directory only), run it will remove primers for the forward and reverse reads. The primer is now set to universal 16S V3-V4 primer set 515F and 806R. Change the -a and -A option in the cut_primer.sh to your primer used during the amplification experiment. The demo_input reads contains no primer. You can also remove primer in this script by providing to it via --forward= and --reverse=

# IMPORTANT
This pipeline is extensively tested under conda R version 4.0.2 on both MacOS and Linux (Ubuntu 18.0.4 and RHEL 7). There are known bugs with old R version on conda in terms of compiling R packages. I strongly suggest you reinstall a new R 4.0.2 from scratch but not update R in you conda. INSTALL A COMPLETELY NEW R 4.0.2! This will save you lot of trouble (updating R packages from an old version to a new one is annoying). For MacOS, you can installed it here: https://cran.r-project.org/bin/macosx/ or by running conda install r-base=4.0.2 -c conda-forge after installing and activating conda.

# DADA2 can now be run with those command (*_R1.fastq, *_R2.fastq or gzipped formatshould be in the input directory):
```
git clone https://github.com/jianshu93/dada2_wrapper
cd dada2_wrapper/scripts

#### Pair end model, default
time Rscript ./dada2_wrapper.r --input_dir=../demo_input --output_dir=output_paired --pool --threads=4
#### Forward reads only, single mode, if you have primer you only need to provide the forward primer on forward reads with --forward="..."
time Rscript ./dada2_wrapper.r --input_dir=../demo_input --output_dir=output_single --pool --single --threads=4

### if you provider primers directly to the wrapper, we can also work it out. But it's slower than cutadapt. Therefore, I still suggest that you use cutadapt first. If you are lasy, you can just do it here. It takes about another 15 minutes to cut primers, very slow compare to cutadapt. 

time Rscript ./dada2_wrapper.r --input_dir=output_paired_primer --output_dir=output1 --pool --forward=GTGCCAGCMGCCGCGGTAA --reverse=GGACTACHVGGGTWTCTAAT --threads=4

### merged mode, merge reads first and then infer based on merged reads using single mode (vsearch with gzip support is need here, script will install vserach from conda automatically if you do not have it. Use the latest 15.2 version)
time Rscript ./dada2_wrapper.r --input_dir=../demo_input --output_dir=output_merged --pool --merge --threads=4

### Or merged mode with primer

time Rscript ./dada2_wrapper.r --input_dir=../demo_input --output_dir=output_merged --pool --merge --forward=GTGCCAGCMGCCGCGGTAA --reverse=GGACTACHVGGGTWTCTAAT --threads 4





```


DADA2 R workflow for profiling 16S sequence reads consists of these 4 scripts:
--------------------------------------------------------------------------------

#### 1. dada2_wrapper.r - This wrapper script parses arguments from the command line and passes them to and invokes the DADA2 workflow (dada2_16S_paired-end.Rmd).

#### 2. dada2_16S_paired-end.Rmd - DADA2 workflow for resolving sequence variants from 16S amplicon reads. In addition to OTU picking and taxonomy assignment, it will also produce a QC folder with visuals and read counts (before/after). 
  NOTE: This workflow requires paired-end 16S rDNA gene amplicon reads.
  NOTE: This workflow will analyze each sample independently to be able to process any sample size.  The option to pool samples for resolving sequence variants will be added.
  NOTE: Currently, the wrapper script dada2_wrapper.r and the workflow dada2_16S_paired-end.Rmd need to be in the same directory (and not in the input directory). By default the database is at the parent directory of the script dada2_wrapper.r, which is: ../reference_dbs_16S

#### 3.dada2_16S_single_end.Rmd - DADA2 workflow for single end 16S amplicon reads. It is very similar to the pair end version but only use forwared reads in the directory (*R1*.fastq). Please name you reads file with *R1* if you have only single end reads data 
#### 4.dada2_16S_merged.Rmd - DADA2 workflow for merged reads. Merged reads first using vsearch -fastq_mergepairs and then run single mode.
#### 5. Diversity.Rmd - Diversity module that calculate both alpha and beta diversity using various metric such as chao1, bray-curtis, unifrac

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

requiredPackages = c('parallel','dplyr','RCurl','ggplot2','GenomicRanges','SummarizedExperiment',
                    'BiocParallel','Rsamtools','dada2','msa','phangorn','gridExtra','knitr')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) BiocManager::install(p,update = FALSE)
}
```
The whole process of install packages will take a few minutes on MacOS (.dmg installed R from R website) but will take much longer on Linux, like ~25 minutes (also ~25 minutes if you are using conda R on MacOS, because you need compiling using conda R). I woud strongly suggest that you use MacOS for small dataset (less than 30 samples) and linux cluster for large dataset.

#### Time needed and platform dependent questions
for the testing dataset, 6 samples, fifty thousand reads per sample, it takes ~8 minutes on a 4-core (8 threads) Macbook Pro but ~15 minutes on a dual-core (4 threads) iMac. For a 72 sample dataset with 200 thousand reads per sample, it takes about 5 hours on a 24 threads Linux system on a cluster node.

#### Running DADA2 on TORQUE based or SLURM based cluster


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
./dada2_wrapper.r --input_dir=../demo_input --output=output --pool

# NOTE: Prior to running make sure to install DADA2 R package dependencies. I have added a few lines to install all packages in the dada2_wrapper.r script after you have successfully install miniconda3 and activate base environment.

```
you will need pandoc and texlive to generate pdf files.

```
conda install pandoc

```

if you do not have textlive installed in your system, the programm will stop and only a .tex file will be generated. No worries! You can transform this tex file to pdf file here: https://cloudconvert.com/tex-to-pdf


For you own dataset just change the input directory to where your pair-end fastq file is. Do remember to remove primers first before running this pipeline. Check the script in the demo_input directory on how to do that. 

Please contact jianshu.zhao@gatech.edu if you have any questions.




# Reference

Bodenhofer, U., Bonatesta, E., & Hochreiter, S. (2015). msa: an R package for multiple sequence alignment. Bioinformatics, 1–3. http://doi.org/10.1093/bioinformatics/btv494

Murali, A., Bhargava, A., & Wright, E. S. (2018). IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences, 1–14. http://doi.org/10.1186/s40168-018-0521-5

Price, M. N., Deha, P. S., & Arkin, A. P. (2010). FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments, 5(3), e9490. http://doi.org/10.1371/journal.pone.0009490

Callahan, B. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13(7), 581–583. http://doi.org/10.1038/nmeth.3869

phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. (2013). phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data, 8(4), e61217–11. http://doi.org/10.1371/journal.pone.0061217

Lozupone, C., & Knight, R. (2005). UniFrac: a new phylogenetic method for comparing microbial communities., 71(12), 8228–8235. http://doi.org/10.1128/AEM.71.12.8228-8235.2005

Chen, J., Bittinger, K., Charlson, E. S., Hoffmann, C., Lewis, J., Wu, G. D., et al. (2012). Associating microbiome composition with environmental covariates using generalized UniFrac distances. Bioinformatics, 28(16), 2106–2113. http://doi.org/10.1093/bioinformatics/bts342

Chao, A., GOTELLI, N. J., Hsieh, T. C., SANDER, E. L., Ma, K. H., COLWELL, R. K., & Ellison, A. M. (2014). Rarefaction and extrapolation with Hill numbers: a framework for sampling and estimation in species diversity studies. Ecological Monographs, 84(1), 45–67. http://doi.org/https://doi.org/10.1890/13-0133.1

Hsieh, T. C., Ma, K. H., & Chao, A. (2016). iNEXT: an R package for rarefaction and extrapolation of species diversity (Hill numbers). Methods in Ecology and Evolution, 7(12), 1451–1456. http://doi.org/10.1111/2041-210X.12613

Shenhav, L., Thompson, M., Joseph, T. A., Briscoe, L., Furman, O., Bogumil, D., et al. (2019). FEAST: fast expectation-maximization for microbial source tracking. Nature Methods, 1–10. http://doi.org/10.1038/s41592-019-0431-x 

Edgar, R. C., & Flyvbjerg, H. (2015). Error filtering, pair assembly and error correction for next-generation sequencing reads. Bioinformatics, 1–7. http://doi.org/10.1093/bioinformatics/btv401

Rognes, T., Flouri, T., Ben Nichols, Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4(17), e2584–22. http://doi.org/10.7717/peerj.2584










