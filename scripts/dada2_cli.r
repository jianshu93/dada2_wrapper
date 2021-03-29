#!/usr/bin/env Rscript
### Jianshu Zhao (jianshu.zhao@gatech.edu). Wrapper for dada2 and IDTAXA for amplicon sequencing analysis
### requireed R packages are dada2, phangorn, msa, ggplot2 and vegan for diveristy analysis. This script will automatically
### install all the required packages and check required external softwares.
### fasttree or FastTreeMP (rename to fasttree) is required for building phylogenetic trees.
### all barcodes and primers shoud be removed for all R1 and R2 reads in the input directory.
### you can use either catadaper or usearch --search_oligo command to remove primers before using this script

## Collect arguments
args <- commandArgs(TRUE)

## Print help message when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      dada2_cli.r
      
      Arguments:
      --input_dir=path/to/dir     - (relative) path to directory containing raw input '.fastq' or '.fastq.gz' files
      --output_dir=path/to/dir    - (relative) path to output directory [default='output']
      --pool                      - logical, if provided then samples are pooled for analysis. It is suggested that the pool mode is
                                    more sensitive than sample-based mode for detecting rare ASVs (http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/)
                                    Default setting is analysis of each sample independently [pool=FALSE].
      --help                      - print this text
      
      Example:
      ./dada2_cli.r --input_dir=input --output_dir=output --pool \n\n")
  
  q(save="no")
}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager",repos='http://cran.us.r-project.org')
library(BiocManager)
  ### check if packages exists and then install packages that does not exist
requiredPackages = c('RCurl','vegan','GUniFrac','ggplot2','phyloseq','GenomicRanges','SummarizedExperiment',
                    'BiocParallel','Rsamtools','dada2','msa','phangorn','gridExtra','rmarkdown','knitr')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) BiocManager::install(p,update = FALSE)
}
### install dependencies using conda
if(Sys.which("pandoc") == "") {
  system("conda install pandoc -y")
}


### load libraries
library("knitr")
# If sample pooling logical is not provided, then default to FALSE, i.e. for independent analysis of each sample.
# Otherwise, if set then TRUE for pooling samples.
if( any( grepl("--pool", args) ) ) {
  args <- gsub( "--pool", "--pool=TRUE", args)
} else {
  args <- c(args, "--pool=FALSE")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args.df <- as.data.frame(do.call("rbind", parseArgs(args)))

args.list <- as.list(as.character(args.df$V2))
names(args.list) <- args.df$V1

## Arg1 default
if(is.null(args.list$input_dir)) {
  stop("input directory must be supplied (input folder).\n", call.=FALSE)
}
if(is.null(args.list$output_dir)) {
  stop("Output directory must be supplied (output folder).\n", call.=FALSE)
}
# Print args list to STDOUT
if(length(args) > 1) {
  for( i in names(args.list) ) {
    cat( i, "\t", args.list[[i]], "\n")
  }
}

# print contents of folder
cat( grep( "*\\.fastq", list.files(args.list$input_dir), value=T ), sep = "\n" )

# these variables are passed to the workflow
input.path <- normalizePath( args.list$input_dir )

output.dir <- ifelse( is.null(args.list$output_dir), "output", args.list$output_dir )

if (!dir.exists(output.dir)){
dir.create(output.dir)
} else {
    print("Dir already exists! If you do not want to overwrite it, please use a new directory")
}
pool.samples <- args.list$pool
# Run dada2 Rmarkdown workflow and output report using Rmarkdown
rmarkdown::render("dada2_16S_paired-end.Rmd",
                  output_file = paste( output.dir, "/16Sreport_dada2_", Sys.Date(), ".pdf", sep='')
)


