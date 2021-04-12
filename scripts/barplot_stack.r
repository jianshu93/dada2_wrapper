#!/usr/bin/env Rscript

# Summarize sequence variant (OTU) count table, from dada2 workflow, through taxonomic ranks.  
# Note that taxonomic ranks must be separate tab-delimited columns at the right end of the table, as output from the dada2 workflow. 

## Collect arguments
args <- commandArgs(TRUE)

## Print help message when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      barplot_stack.r
      
      Arguments:
      --input=path/to/file        - Path to input abundance/relative abundance table; tab-delimited.
                                    Samples as columns, feature (phylum, class et.al., output from summarize_taxa.r) as rows.
                                    Taxonomy ranks at right end of table as separate columns.
      --output_file=path/to/file    - Path to output file [default='barplot_stack.pdf']
      --help                      - print this text
      
      Example:
      ./summarize_taxa.r --input=input --output_file=barplot_stack.pdf --top=25 \n\n")
  
  q(save="no")
}

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager",repos='http://cran.us.r-project.org')
library(BiocManager,quietly=TRUE)
  ### check if packages exists and then install packages that does not exist
requiredPackages = c('reshape2','dplyr',,'ggplot2','cowplot','ggalluvial')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE, quietly=TRUE)) BiocManager::install(p,update = FALSE,quietly=TRUE)
}
}

library("reshape2")
library("ggplot")
library("ggalluvial")

# Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args.df <- as.data.frame(do.call("rbind", parseArgs(args)))

args.list <- as.list(as.character(args.df$V2))
names(args.list) <- args.df$V1

# Arg1 default
if(is.null(args.list$input)) {
  stop("At least one argument must be supplied (input folder).\n", call.=FALSE)
}

# If output directory is not provided then make default
output.dir <- ifelse( is.null(args.list$output_file), "barplot_stack.pdf", args.list$output_file )
