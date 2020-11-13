## Script to compare the tree topologies usnig RF distances for the datasets were reanalyzed following manual inspection

library(phangorn)
library(rwty)

### Set up a bunch of directories
##
## this will need to be changed to reflect where your files are stored
##
main_dir<-"02_Reanalyses/Manual" # directory containing this script and 01_R_functions_calculate_diagnostics_V2.R
storage_old<-"02_Reanalyses/Manual/chains/orig_chains"  # Directory with the original chains for comparison
storage_HKY<-"02_Reanalyses/Manual/chains/HKY"         # Directory for chains reanalyzed using HKY
storage_empbfreq<-"02_Reanalyses/Manual/chains/EmpBF"   # Directory for chains reanalyzed using empirical base frequencies
storage_run_longer<-"02_Reanalyses/Manual/chains/RunLonger" # Directory for chains reanalyzed by running longer and retaining more samples
storage_NstMixed<-"02_Reanalyses/Manual/chains/NstMixed"   # Directory for chains reanalyzed using Nst=mixed model averaging over subsitution models
out_dir<-"02_Reanalyses/Manual/manual_output"             # Directory to write output files to

setwd(main_dir)

trim<-1    ## Trim is the nth tree to be kept for thinning out MCMC - do this to prevent overly large lists of trees when necessary
target_samp<-1000  # we want to retain 1000 samples
max_gens<-10000000  # we want all chains to be read in up to the 10,000,000th generation so that all are the same
        ### ^ however, we will change this for the analyses that we run longer to allow there to be more samples in these analyses

source("01_R_functions_calculate_diagnostics_V2.R")


# Read in the new chains HKY
setwd(storage_HKY)
folders_HKY<-as.list(list.files(pattern=".tar.gz"))
chains_HKY<-lapply(folders_HKY, read.comp.catch, trim=trim, target_samp=target_samp, max_gens=max_gens)

# Read in the new chains with empirical base frequencies
setwd(storage_empbfreq)
folders_empbfreq<-as.list(list.files(pattern=".tar.gz"))
chains_empbfreq<-lapply(folders_empbfreq, read.comp.catch, trim=trim, target_samp=target_samp, max_gens=max_gens)

# Read in the new chains that were run longer - multiple max gens by 3, since they were run 3x longher
setwd(storage_run_longer)
folders_run_longer<-as.list(list.files(pattern=".tar.gz"))
chains_run_longer<-lapply(folders_run_longer, read.comp.catch, trim=trim, target_samp=target_samp, max_gens=(3*max_gens))
## also read in a version that retains 3x as many samples since these were run for 3x longer
chains_run_longer_more<-lapply(folders_run_longer, read.comp.catch, trim=trim, target_samp=(3*target_samp), max_gens=(3*max_gens))

# Read in the chains with gamma
setwd(storage_gamma)
folders_gamma<-as.list(list.files(pattern=".tar.gz"))
chains_gamma<-lapply(folders_gamma, read.comp.catch, trim=trim, target_samp=target_samp, max_gens=max_gens)

# Read in the new chains with Nst=mixed
setwd(storage_NstMixed)
folders_NstMixed<-as.list(list.files(pattern=".tar.gz"))
chains_NstMixed<-lapply(folders_NstMixed, read.comp.catch, trim=trim, target_samp=target_samp, max_gens=max_gens)

## Read in the old chains
setwd(storage_old)
folders_old<-as.list(list.files(pattern=".tar.gz"))
chains_old<-lapply(folders_old, read.comp.catch, trim=trim, target_samp=target_samp, max_gens=max_gens) # read in the chains

# Save a checkpoint here because reading chains takes a while
setwd(main_dir)
# save.image(file="chains_for_tree_compare.RData")
# load("chains_for_tree_compare.RData")


#############################################################################
#############################################################################
#############     Compare tree distances    #################################
#############################################################################
#############################################################################


### Make a function to compare the distances of trees among runs
### This function computes the normalized RF distance between trees from pre and post reanalysis
### calculates the RF distances for 95% consensus trees, 50% consensus trees, and MCC trees
### trees from multiple chains of the same analysis are lumped together
compare_rf<-function(reanalyses, originals, out_file, new_burn, old_burn){   
 ## Match names of reanalyses to original chains
  names(reanalyses)<-sapply(reanalyses, function(x) gsub("\\.run.*$|\\.nex.*$|_processed_mbReady.*$", "",names(x)[[1]])) # get the names of the reanalyses
  names(originals)<-sapply(originals, function(x) gsub("\\.run.*$|\\.nex.*$", "",names(x)[[1]])) # get the names of the analyses
  names_match<-names(reanalyses)[which(names(reanalyses) %in% names(originals))] # get the names that match up - 1 analysis failed to reanalyze at some point
  # combine chains
  old_new_combined<-lapply(names_match, function(x) c(unlist(originals[x], recursive=FALSE), unlist(reanalyses[x], recursive=FALSE)))
  res<-list()
  for(i in 1:length(old_new_combined)){
      ind_old<-length(old_new_combined[[i]])/2  # get the index of the last of the old analyses in the old_new_combined object
      old_only<-old_new_combined[[i]][1:ind_old]  # get the first half of the analyese - these are the old ones
      new_only<-old_new_combined[[i]][(ind_old+1):length(old_new_combined[[1]])]  # Get the second half of the anlyses - these are the new ones
      old_start_tree<-(old_burn/100)*(length(old_only[[1]][[1]])-1)+1 ## Get the number of the first post-burnin tree 
      new_start_tree<-(new_burn/100)*(length(new_only[[1]][[1]])-1)+1 ## Get the number of the first post-burnin tree 
      # Get a consensus tree for the old analysis across all chains
      old_trees_pre<-lapply(old_only, function(x) x$trees[old_start_tree:length(x$trees)])  # Get all of the post-burnin old trees across chains
      old_trees<-Reduce(c, old_trees_pre) # Each chain is a separate element of the list - combine these together into a single multiphylo
      old_con_95<-consensus(old_trees, p=0.95)  # Get a consensus tree out for comparison
      # Same for the new trees
      new_trees_pre<-lapply(new_only, function(x) x$trees[new_start_tree:length(x$trees)])  # Get all of the post-burnin old trees across chains
      new_trees<-Reduce(c, new_trees_pre)
      new_con_95<-consensus(new_trees, p=0.95)
      #Calculate the RF distance for 95% consensus trees
      norm_RF_95<-RF.dist(old_con_95, new_con_95, normalize = TRUE)
      #Calculate the RF distance for 50% consensus trees
      old_con_50<-consensus(old_trees, p=0.5) 
      new_con_50<-consensus(new_trees, p=0.5)
      norm_RF_50<-RF.dist(old_con_50, new_con_50, normalize = TRUE)
      # Calculate distance among MCC trees
      old_con_MCC<-maxCladeCred(old_trees, rooted = FALSE) 
      new_con_MCC<-maxCladeCred(new_trees, rooted = FALSE)
      norm_RF_MCC<-RF.dist(old_con_MCC, new_con_MCC, normalize = TRUE)
      ## List out all of the output we want to retain
      old_con_trees<-list(old_con_95, old_con_50, old_con_MCC)
      names(old_con_trees)<-c("Old_95_con", "Old_50_con", "Old_MCC")
      new_con_trees<-list(new_con_95, new_con_50, new_con_MCC)
      names(new_con_trees)<-c("New_95_con", "New_50_con", "New_MCC")
      trees_res<-list(old_con_trees, new_con_trees)
      names(trees_res)<-c("Old_summ_trees", "New_summ_trees")
      rf_dists<-c(norm_RF_95, norm_RF_50, norm_RF_MCC)
      names(rf_dists)<-c("norm_RF_95", "norm_RF_50", "norm_RF_MCC")
      res[[i]]<-list(rf_dists, trees_res)
  }  
  return(res)
}  

## Get normalized RF distances for pre and post HKY
HKY_RF<-compare_rf(reanalyses = chains_HKY, originals = chains_old, new_burn = 25, old_burn = 25)

## RF distances pre and post empirical base frequency
empbfreq_RF<-compare_rf(reanalyses = chains_empbfreq, originals = chains_old, new_burn = 25, old_burn = 25)

### RF distances comparing run longer with more samples
run_longerMore_RF<-compare_rf(reanalyses = chains_run_longer_more, originals = chains_old, new_burn = 25, old_burn = 25)

## RF pre and post changing gamma
gamma_RF<-compare_rf(reanalyses = chains_gamma, originals = chains_old, new_burn = 25, old_burn = 25)

## RF pre and post changing gamma
NstMixed_RF<-compare_rf(reanalyses = chains_NstMixed, originals = chains_old, new_burn = 25, old_burn = 25)



rf_comparisons<-list(HKY_RF, empbfreq_RF, run_longerMore_RF, gamma_RF, NstMixed_RF)
names(rf_comparisons)<-c("HKY_RF", "empbfreq_RF", "run_longerMore_RF", "gamma_RF", "NstMixed_RF")

save(rf_comparisons, file="rf_comp_Manual.RData")


################################################################
#### Summarize the RF distances
################################################################
load("rf_comp_Manual.RData")

HKY_RF<-rf_comparisons[[1]]
empbfreq_RF<-rf_comparisons[[2]]
run_longerMore_RF<-rf_comparisons[[3]]
gamma_RF<-rf_comparisons[[4]]
NstMixed_RF<-rf_comparisons[[5]]

## Get out the RF scores only from each of these objects and plot them out
HKY_RF_only<-do.call(rbind, lapply(HKY_RF, function(x) x[[1]]))
empbfreq_RF_only<-do.call(rbind, lapply(empbfreq_RF, function(x) x[[1]]))
run_longerMore_RF_only<-do.call(rbind, lapply(run_longerMore_RF, function(x) x[[1]]))
gamma_RF_only<-do.call(rbind, lapply(gamma_RF, function(x) x[[1]]))
NstMixed_RF_only<-do.call(rbind, lapply(NstMixed_RF, function(x) x[[1]]))



### Plots of the change between new and old heating for 95% and 50% consensus trees and Maximum Clade Credibility (MCC) trees
setwd(out_dir)
pdf(file="RF manual.pdf", height=10, width=11)
par(mfrow=c(4,3))
# HKY
hist(HKY_RF_only[,1], col="gray", main=paste("HKY 95% con \nMedian = ", round(median(HKY_RF_only[,1]), 3),sep=""), xlab=NULL, xlim=c(0,1), ylab=NULL, breaks=6)
hist(HKY_RF_only[,2], col="gray", main=paste("HKY 50% con \nMedian = ", round(median(HKY_RF_only[,2]), 3),sep=""), xlab=NULL, xlim=c(0,1), ylab=NULL, breaks=6)
hist(HKY_RF_only[,3], col="gray", main=paste("HKY MCC \nMedian = ", round(median(HKY_RF_only[,3]), 3),sep=""), xlab=NULL, xlim=c(0,1), ylab=NULL, breaks=10)
# EmpBFreq
hist(empbfreq_RF_only[,1], col="gray", main=paste("Empirical Base Freq 95% con \nMedian = ", round(median(empbfreq_RF_only[,1]), 3), sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL, breaks=6)
hist(empbfreq_RF_only[,2], col="gray", main=paste("Empirical Base Freq 50% con \nMedian = ", round(median(empbfreq_RF_only[,2]), 3), sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL, breaks=6)
hist(empbfreq_RF_only[,3], col="gray", main=paste("Empirical Base Freq MCC \nMedian = ", round(median(empbfreq_RF_only[,3]), 3), sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL, breaks=6)
# Run Longer more samples
hist(run_longerMore_RF_only[,1], col="gray", main=paste("Run longer 95% con \nMedian = ", round(median(run_longerMore_RF_only[,1]), 3), sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL, breaks=1)
hist(run_longerMore_RF_only[,2], col="gray", main=paste("Run longer 50% con \nMedian = ", round(median(run_longerMore_RF_only[,2]), 3), sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL, breaks=1)
hist(run_longerMore_RF_only[,3], col="gray", main=paste("Run longer MCC \nMedian = ", round(median(run_longerMore_RF_only[,3]), 3), sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL, breaks=2)
# NstMixed
hist(NstMixed_RF_only[,1], col="gray", main=paste("Nst=Mixed 95% con \nMedian = ", round(median(NstMixed_RF_only[,1]), 3), sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL, breaks=3)
hist(NstMixed_RF_only[,2], col="gray", main=paste("Nst=Mixed 50% con \nMedian = ", round(median(NstMixed_RF_only[,2]), 3), sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL, breaks=3)
hist(NstMixed_RF_only[,3], col="gray", main=paste("Nst=Mixed MCC \nMedian = ", round(median(NstMixed_RF_only[,3]), 3), sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL, breaks=8)
dev.off()# turn off the pdf plotting device

