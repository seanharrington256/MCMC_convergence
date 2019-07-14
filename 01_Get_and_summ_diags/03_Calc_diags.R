## Script to use compressed folders - uncompresses reads in chains, then deletes uncompressed folder
## loops convergence diagnostics over chains that are read ins
## Includes checkpointing to save the output every x number of chains that have been analyzed to track progress better

## Set the batch number, to determine which folders to analyze and propoerly set the output file name
## Do this using an argument passed in through from the terminal (or shell script)
batch <- commandArgs(trailingOnly = TRUE)

rwty.processors<-1
batch_size<-10  # set the number of folders to iterate over in each batch

# Set up the main directory with the scripts and the directory with the folders to be iterated over
main_dir<-getwd()
storage<-"~/apps/1_All_orig_chains"  ## Set this to the location of the chains we are working with


trim<-1    ## Trim is the nth tree to be kept for thinning out MCMC - do this to prevent overly large lists of trees when necessary
target_samp<-1000  # we want to retain 1000 samples
max_gens<-10000000  # we want all chains to be read in up to the 10,000,000th generation so that all are the same


library(rwty)
library(coda)
library(boa)
library(gridExtra)
library(adegenet)
library(adegraphics)
library(phytools)
library(parallel)
library(tracerer)
library(gsubfn)


setwd(main_dir)

source("01_R_functions_calculate_diagnostics_V2.R")
source("02_Wrapper_func_calc_diags.R")

setwd(storage)


# List the gz folders
folders_all<-as.list(list.files(pattern=".tar.gz"))
if(length(folders_all)>batch_size){
  batches<-split(folders_all, ceiling(seq_along(folders_all)/batch_size)) # Split into batches of whatever size specified by batch size above
  folders<-batches[[batch]]#select the specific folders to analyze in each
}else{
  folders<-folders_all
}

# Split the folders into batches of 10 to be read in, analyzed, and removed to prevent storing too many trees in memory at once
split_folders<-split(folders, ceiling(seq_along(folders)/10))

# Function to read in some trees then remove them to prevent using up huge amounts of memory    
summ.and.rm.chains<-function(folders, trim=1, target_samp=1000, max_gens=10000000){
  # read in the chains from the compressed folders
  all_info<-lapply(folders, read.comp.info.catch, trim=trim, target_samp=target_samp, max_gens=max_gens)
  chains<-lapply(all_info, function(x) x[[1]])
  info<-lapply(all_info, function(x) x[2:length(x)])
  # Get convergence summary from these
  summ<-lapply(chains, conv.summ.catch, topo_ESS_reps=20)
  res<-list()
  for(i in 1:length(summ)){
    res[[i]]<-c(summ[[i]], info[[i]])
  }
  rm(list=c("chains", "all_info"))
  return(res)
}


summ_pre<-list()
for(i in 1:length(split_folders)){
  summ_pre[[i]]<-summ.and.rm.chains(split_folders[[i]], trim=trim, target_samp=target_samp, max_gens=max_gens)
  
  # To track progress and have results if there is a failure, write a checkpoint each
  # time that one batch of chains is finished
  setwd(main_dir)
  ckpt<-summ_pre[[i]]
  save(ckpt, file=paste("batch_", batch, "_ckpt_", i, ".RData", sep=""))
  if(i>1){
    unlink(paste("batch_", batch, "_ckpt_", i-1, ".RData", sep=""), recursive=TRUE)  # Remove the old checkpoint to prevent file ballooning
  }
  setwd(storage)
}

summ<-unlist(summ_pre, recursive=FALSE)

setwd(main_dir)

save(summ, file=paste("summary_", batch, ".RData", sep=""))


