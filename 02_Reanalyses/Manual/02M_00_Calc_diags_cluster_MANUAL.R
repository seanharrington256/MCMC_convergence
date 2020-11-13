### This is a modified version of 03_Calc_diags.R to get diagnostics for chains that were reanalyzed after manual examination

# primary modification is that instead of the second command line argument specifying the directory of chains
#    to be iterated over, it specifies an analysis set, which then selects among several directories and 
#    changes setting depending. The analyses to run longer require different settings than other analyses do

# Otherwise, pretty much the same.


# as is that case for that script, this requires that the scripts "01_R_functions_calculate_diagnostics_V2.R" and "02_Wrapper_func_calc_diags.R" are
##   in the directory that this script is run from

## The directories on lines 49-65 will need to be changed to reflect the locations of each set of analyses to be run on your setup


# print out the command line args used
print(paste0("batches = ", commandArgs(trailingOnly = TRUE)[[1]]))
print(paste0("analysis_set = ", commandArgs(trailingOnly = TRUE)[[2]]))

## Set the batch number, to determine which folders to analyze and properly set the output file name
## Do this using an argument passed in through from the terminal (or shell script)
batch<-commandArgs(trailingOnly = TRUE)[[1]]

rwty.processors<-1
batch_size<-10  # set the number of folders to iterate over in each batch

# Set up the main directory with the scripts and the directory with the folders to be iterated over
main_dir<-getwd()
analysis_set<-commandArgs(trailingOnly = TRUE)[[2]]  ## Set this to the location of the chains we are working with


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

if(analysis_set=="HKY"){
  storage<-"MCMC_reanalyses/manual/HKY"
}
if(analysis_set=="EmpBF"){
  storage<-"MCMC_reanalyses/manual/empBF"
}
if(analysis_set=="RunLonger"){
  storage<-"MCMC_reanalyses/manual/RunLonger"
}
if(analysis_set=="gamma"){
  storage<-"MCMC_reanalyses/manual/gamma"
}
if(analysis_set=="NstMixed"){
  storage<-"MCMC_reanalyses/manual/NstMixed"
}
if(analysis_set=="EmpBFLong"){
  storage<-"MCMC_reanalyses/manual/EmpBFLong"
}



if(analysis_set %in% c("HKY", "EmpBF", "gamma", "NstMixed")){
  trim<-1    ## Trim is the nth tree to be kept for thinning out MCMC - do this to prevent overly large lists of trees when necessary
  target_samp<-1000  # we want to retain 1000 samples
  max_gens<-10000000  # we want all chains to be read in up to the 10,000,000th generation so that all are the same
}
if(analysis_set %in% c("RunLonger", "EmpBFLong")){
  trim<-1    ## Trim is the nth tree to be kept for thinning out MCMC - do this to prevent overly large lists of trees when necessary
  target_samp<-3000  # we want to retain 3,000 samples
  max_gens<-30000000  # Here we read in chains up to 30,000,000 M instead of just 10 M because these chains needed to be run longer
}


#####
##### Get out the diags
#####
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
  summ<-lapply(chains, conv.summ.catch, topo_ESS_reps=20, burnin=25)
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
    if(analysis_set=="HKY"){
      save(ckpt, file=paste("batch_", batch, "_ckpt_HKY_", i, ".RData", sep=""))
      if(i>1){
        unlink(paste("batch_", batch, "_ckpt_HKY_", i-1, ".RData", sep=""), recursive=TRUE)  # Remove the old checkpoint to prevent file ballooning
      }
    }
    if(analysis_set=="EmpBF"){
      save(ckpt, file=paste("batch_", batch, "_ckpt_EmpBF_", i, ".RData", sep=""))
      if(i>1){
        unlink(paste("batch_", batch, "_ckpt_EmpBF_", i-1, ".RData", sep=""), recursive=TRUE)  # Remove the old checkpoint to prevent file ballooning
      }
    }
    if(analysis_set=="RunLonger"){
      save(ckpt, file=paste("batch_", batch, "_ckpt_RunLonger_", i, ".RData", sep=""))
      if(i>1){
        unlink(paste("batch_", batch, "_ckpt_RunLonger_", i-1, ".RData", sep=""), recursive=TRUE)  # Remove the old checkpoint to prevent file ballooning
      }
    }
    if(analysis_set=="gamma"){
      save(ckpt, file=paste("batch_", batch, "_ckpt_gamma_", i, ".RData", sep=""))
      if(i>1){
       unlink(paste("batch_", batch, "_ckpt_gamma_", i-1, ".RData", sep=""), recursive=TRUE)  # Remove the old checkpoint to prevent file ballooning
      }
    }
    if(analysis_set=="NstMixed"){
      save(ckpt, file=paste("batch_", batch, "_ckpt_NstMixed_", i, ".RData", sep=""))
      if(i>1){
        unlink(paste("batch_", batch, "_ckpt_NstMixed_", i-1, ".RData", sep=""), recursive=TRUE)  # Remove the old checkpoint to prevent file ballooning
      }
    }
    if(analysis_set=="EmpBFLong"){
      save(ckpt, file=paste("batch_", batch, "_ckpt_EmpBFLong_", i, ".RData", sep=""))
      if(i>1){
        unlink(paste("batch_", batch, "_ckpt_EmpBFLong_", i-1, ".RData", sep=""), recursive=TRUE)  # Remove the old checkpoint to prevent file ballooning
      }  
    }
    
    
setwd(storage)
}

summ<-unlist(summ_pre, recursive=FALSE)

setwd(main_dir)

if(analysis_set=="HKY"){
  save(summ, file=paste("summary_HKY_", batch, ".RData", sep=""))
}
if(analysis_set=="EmpBF"){
  save(summ, file=paste("summary_EmpBF_", batch, ".RData", sep=""))
}
if(analysis_set=="RunLonger"){
  save(summ, file=paste("summary_RunLonger_", batch, ".RData", sep=""))
}
if(analysis_set=="gamma"){
  save(summ, file=paste("summary_gamma_", batch, ".RData", sep=""))
}
if(analysis_set=="NstMixed"){
  save(summ, file=paste("summary_NstMixed_", batch, ".RData", sep=""))
}
if(analysis_set=="EmpBFLong"){
  save(summ, file=paste("summary_EmpBFLong_", batch, ".RData", sep=""))
}




