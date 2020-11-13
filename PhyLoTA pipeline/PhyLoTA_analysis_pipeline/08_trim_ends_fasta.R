#!/usr/bin/env Rscript
args_trim = commandArgs(trailingOnly=TRUE)

## This script is just a quick little wrapper around the ips trimEnds function
##    to trim and ragged ends from a FASTA file
##    takes 3 command line arguments: input file, output file, and a minimum numer of sequences
##    or fraction of sequences that must have non missing/gapped/ambiguous data
## run this from command line using Rscript --vanilla Trim_ends_fasta.R <input_path> <output_path> <threshold_for_trimming>

## Set up the arguments from command line
infile<-args_trim[1] # input fasta file
outfile<-args_trim[2] # output for trimmed file
thresh<-as.numeric(args_trim[3]) # threshold for missing data at ends

library(ips)

input<-read.FASTA(infile) # read in the fasta file
input<-as.matrix(input) # convert it into a matrix for trimEnds
out <- trimEnds(input, min.n.seq = thresh) # trim the ends to the desired fraction
write.FASTA(out, outfile)# write the output file


