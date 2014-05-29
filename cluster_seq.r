#!/usr/bin/env Rscript

# Load and if needed install libraries
if(suppressMessages(!require(optparse))) {
  print("Trying to install optparse.")
  install.packages(optparse)
  if(require(optparse)){
    print("Optparse installed and loaded.")
  } else {
    stop("Could not install optparse")
  }
}

if(suppressMessages(!require(Biostrings))) {
  print("Trying to install Biostrings.")
  install.packages(Biostrings)
  if(require(Biostrings)){
    print("Biostrings installed and loaded.")
  } else {
    stop("Could not install Biostrings.")
  }
}

# Define options
option_list <- list(
  make_option(c("-m", "--method"),
              type="character",
              default="complete",
              help="Method for hclust")
)

# Define parser
parser = OptionParser(usage = "%prog [options] FASTA OUTPUT", option_list=option_list)

# Parse arguments and options
arguments = parse_args(parser, positional_arguments = 2)
opt = arguments$options
args = arguments$args

# Loads the FASTA file
reads = readDNAStringSet(args[1], format="fasta")

# Hierarchical Clustering using Levensthein distance
clusters <- hclust(stringDist(reads), method=opt$method)

# Save
clusters_file_name <- paste("pansen", ".clusters.Rdata", sep="")
save(file=clusters_file_name, clusters)

png_file_name <- paste(clusters_file_name, ".png", sep="")
png(file=png_file_name, heigh=800, width=1600, type="cairo")
plot(clusters)

cat(png_file_name, "\n")