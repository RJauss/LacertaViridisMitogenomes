#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("requires 2 parameters: mafft aligned fasta filename and output filename\n", call.=FALSE)
}

# using genetics library
# install.packages("ape") # only uncomment, if not already installed
library('ape')

# load fasta file
alignment <- read.dna(args[1], format = 'fasta')

# calculate bp size of one sequence
sequence_size <- length(alignment) / length(labels(alignment))

# create nex file
write.nexus.data(alignment, file=args[2], format="dna", missing="?", gap="-",charsperline=sequence_size)
