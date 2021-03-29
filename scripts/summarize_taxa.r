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
      summarize_taxa.r
      
      Arguments:
      --input=path/to/file        - Path to input sequence variant count table; tab-delimited.
                                    Samples as columns, sequence variants (OTUs) as rows.
                                    Taxonomy ranks at right end of table as separate columns.
      --output_dir=path/to/dir    - Path to output directory [default='taxa_summary']
      --help                      - print this text
      
      Example:
      ./summarize_taxa.r --input=input \n\n")
  
  q(save="no")
}

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
output.dir <- ifelse( is.null(args.list$output_dir), "taxa_summary", args.list$output_dir )

# Print args list to STDOUT
for( i in names(args.list) ) {
  cat( i, "\t", args.list[[i]], "\n")
}


# Read in counts table, samples as columns and sequence variants (or OTUs) as rows
otu.dat <- read.table( args.list$input, sep = "\t", header = T, row.names = 1, stringsAsFactors = FALSE)

# column index of rank kingdom 
k_index <- which( colnames( otu.dat ) == "Kingdom" )
# integer lowest taxonomic rank; e.g. Species is 7
lowest_rank <- length( otu.dat[ k_index : ncol( otu.dat ) ] )

# Create directory to store output files
taxa.folder <- file.path( getwd(), output.dir )
ifelse(!dir.exists(taxa.folder), dir.create(taxa.folder, recursive = TRUE), FALSE)

# loop through ranks, summarize taxa at each rank, sum normalize, and save to file
for( i in 1 : (lowest_rank-1) ) {  # substract 1 because it already starts at kingdom, so kingdom + 6 = species
  tmp.counts <- otu.dat
  tmp.counts$taxonomy <- apply( tmp.counts[ , k_index : (k_index + i) ], 1, function(x) paste( x, collapse = ';' ) )
  tmp.counts <- tmp.counts[ , -c(k_index : (ncol(tmp.counts) - 1) ) ]
  tmp.counts <- aggregate( x = tmp.counts[,-ncol(tmp.counts)], by = list(tmp.counts$taxonomy), FUN = 'sum' )
  rownames( tmp.counts ) <- tmp.counts[,1]
  tmp.counts <- tmp.counts[,-1]
  write.table( tmp.counts, paste0( taxa.folder, "/tax.table_L", i+1, "_counts.txt" ), sep = "\t", quote = F, eol = "\n", col.names = NA )
  # sum normalize
  tmp.relab <- apply( tmp.counts, 2, function(x) x/sum(x) )
  write.table( tmp.relab, paste0( taxa.folder, "/tax.table_L", i+1, "_relab.txt" ), sep = "\t", quote = F, eol = "\n", col.names = NA )
}
