## Script to compare the tree topologies for the datasets were reanalyzed with different heating or with nst=mixed or with IG switched to just G

## Set the batch number, to determine which folders to analyze and propoerly set the output file name
## Do this using an argument passed in through from the terminal (or shell script)
batch <- commandArgs(trailingOnly = TRUE)[[1]]
### Set up here if we're analyzing IG vs G ("G"), Nst=mixed ("mixed"), or New Heating ("NewHeat")
analysis_set<-commandArgs(trailingOnly = TRUE)[[2]]


rwty.processors<-1
batch_size<-10  # set the number of folders to iterate over in each batch


library(phangorn)
library(rwty)

# Set up the main directory with the scripts and the directory with the folders to be iterated over
main_dir<-getwd()
storage_old<-"1_All_orig_chains"
storage_new_heat<-"Nex_new_heat_chains_complete/RESULTS_tar_gz"
storage_nst_mixed<-"Nex_NstMixed_complete/RESULTS_tar_gz"
storage_G_only<-"Nex_IG_reanalyses_complete/RESULTS_tar_gz"
# Out directories
out_IG<-paste(main_dir, "IG_rf_out", sep="/")
out_NewHeat<-paste(main_dir, "NewHeat_rf_out", sep="/")
out_nstMixed<-paste(main_dir, "nstMixed_out", sep="/")

setwd(main_dir)

trim<-1    ## Trim is the nth tree to be kept for thinning out MCMC - do this to prevent overly large lists of trees when necessary
target_samp<-1000  # we want to retain 1000 samples
max_gens<-10000000  # we want all chains to be read in up to the 10,000,000th generation so that all are the same


source("01_R_functions_calculate_diagnostics_V2.R")

## originally, we included analyses that had failed Geweke's test but no other diagnostics
# Decided to not deal with analyses that failed Geweke's, so remove those here
# but because of potential overlap among analyses in the different subsets, can't just remove the Geweke's ones, instead keep everything else
bad_ESS<-read.csv("fail_oneormore_ESS_LARGE.csv", stringsAsFactors = FALSE, row.names=1)[[1]]
bad_topo_bad_LnL<-read.csv("fail_topo_bad_LnL_ESS_LARGE.csv", stringsAsFactors = FALSE, row.names=1)[[1]]
bad_topo_good_LnL<-read.csv("fail_topo_good_LnL_ESS_LARGE.csv", stringsAsFactors = FALSE, row.names=1)[[1]]
to_keep<-c(bad_ESS, bad_topo_bad_LnL, bad_topo_good_LnL)



if(analysis_set=="NewHeat"){
  # Read in the new chains new heat
  setwd(storage_new_heat)
  folders_new_heat_all<-as.list(list.files(pattern=".tar.gz"))
  folders_new_heat_all<-folders_new_heat_all[which(gsub(".tar.gz", "", unlist(folders_new_heat_all)) %in% to_keep)] ## Get just the non-Geweke ones to keep
  batches<-split(folders_new_heat_all, ceiling(seq_along(folders_new_heat_all)/batch_size)) # Split into batches of whatever size specified by batch size above
  folders_new_heat<-batches[[batch]]#select the specific folders to analyze in each
  chains_new_heat<-lapply(folders_new_heat, read.comp.catch, trim=trim, target_samp=target_samp, max_gens=max_gens)
  failed_read<-chains_new_heat[grep("FAILED", chains_new_heat)]
  if(length(failed_read)>0){
    chains_new_heat<-chains_new_heat[-grep("FAILED", chains_new_heat)]
  }
  folders<-folders_new_heat
}


if(analysis_set=="mixed"){
  # Read in the new chains nst mixed
  setwd(storage_nst_mixed)
  folders_nst_mixed_all<-as.list(list.files(pattern=".tar.gz"))
  folders_nst_mixed_all<-folders_nst_mixed_all[which(gsub(".tar.gz", "", unlist(folders_nst_mixed_all)) %in% to_keep)] ## Get just the non-Geweke ones to keep
  batches<-split(folders_nst_mixed_all, ceiling(seq_along(folders_nst_mixed_all)/batch_size)) # Split into batches of whatever size specified by batch size above
  folders_nst_mixed<-batches[[batch]]#select the specific folders to analyze in each
  chains_nst_mixed<-lapply(folders_nst_mixed, read.comp.catch, trim=trim, target_samp=target_samp, max_gens=max_gens)
  failed_read<-chains_nst_mixed[grep("FAILED", chains_nst_mixed)]
  if(length(failed_read)>0){
    chains_nst_mixed<-chains_nst_mixed[-grep("FAILED", chains_nst_mixed)]
  }
  folders<-folders_nst_mixed
}



if(analysis_set=="G"){
  # Read in the new chains G only
  setwd(storage_G_only)
  folders_G_only_all<-as.list(list.files(pattern=".tar.gz"))
  batches<-split(folders_G_only_all, ceiling(seq_along(folders_G_only_all)/batch_size)) # Split into batches of whatever size specified by batch size above
  folders_G_only<-batches[[batch]]#select the specific folders to analyze in each
  
  chains_G_only<-lapply(folders_G_only, read.comp.catch, trim=trim, target_samp=target_samp, max_gens=max_gens)
  failed_read<-chains_G_only[grep("FAILED", chains_G_only)]
  if(length(failed_read)>0){
    chains_G_only<-chains_G_only[-grep("FAILED", chains_G_only)]
  }
  folders<-folders_G_only
}


## Make a list of all of the original chains that need to be read in
setwd(storage_old) # set working directory
old_to_read<-list.files()[which(gsub("_READY", "", list.files()) %in% unique(folders))]

## Read in the old chains
chains_old<-lapply(old_to_read, read.comp.catch, trim=trim, target_samp=target_samp, max_gens=max_gens) # read in the chains


setwd(main_dir)

if(analysis_set %in% c("mixed", "NewHeat")){
  ## Load in the burnin for each analysis - use the maximum of all chains (new and old)
  load("example_chains_burnin_LARGE.RData")
  if(length(which(duplicated(names(example_chains_burn_LARGE))))>0){
    example_chains_burn_LARGE<-example_chains_burn_LARGE[-which(duplicated(names(example_chains_burn_LARGE)))] # Remove any duplicates
  }
}


if(analysis_set=="G"){
  ## Load in G vs IG burns
  load("example_IG_chains_burn_LARGE.RData")  
  IG_ex_chains_burn_LARGE<-unlist(IG_ex_chains_burn_LARGE, recursive=FALSE)## Combine these into a more usable format for farther down
  names(IG_ex_chains_burn_LARGE)<-gsub("high_IG_cor_bad_ESS_LARGE_burns.|high_IG_cor_good_conv_LARGE_burns.", "",  names(IG_ex_chains_burn_LARGE)) # clip off pieces of names specifying which subset each goes to
  if(length(which(duplicated(names(IG_ex_chains_burn_LARGE))))>0){
    IG_ex_chains_burn_LARGE<-IG_ex_chains_burn_LARGE[-which(duplicated(names(IG_ex_chains_burn_LARGE)))] # Remove any duplicates
  }
} 


## New burn ins
if(analysis_set=="NewHeat"){
  load("NewHeat_burns.RData")
}
if(analysis_set=="mixed"){
  load("nstMixed_burns.RData")
}
if(analysis_set=="G"){
  load("G_only_burns.RData")
}


#############################################################################
#############################################################################
#############     Compare tree distances    #################################
#############################################################################
#############################################################################
#############################################################################



### Make a function to compare the distances of trees among runs
### This function computes the normalized RF distance between trees from pre and post reanalysis
### calculates the RF distances for 95% consensus trees, 50% consensus trees, and MCC trees
### trees from multiple chains of the same analysis are lumped together
###   if the burnin for a chain was not successfully calculated for a reanalyzed chain, return NA
compare_rf<-function(reanalyses, originals, out_file, new_burn, old_burn){   
  ## Match names of reanalyses to original chains
  names(reanalyses)<-sapply(reanalyses, function(x) gsub("\\.run.*$|\\.nex.*$|_processed_mbReady.*$", "",names(x)[[1]])) # get the names of the reanalyses
  names(originals)<-sapply(originals, function(x) gsub("\\.run.*$|\\.nex.*$|_processed_mbReady.*$", "",names(x)[[1]])) # get the names of the analyses
  names_match<-names(reanalyses)[which(names(reanalyses) %in% names(originals))] # get the names that match up - analyses may have failed to reanalyze at some point
  # combine chains
  old_new_combined<-lapply(names_match, function(x) c(unlist(originals[x], recursive=FALSE), unlist(reanalyses[x], recursive=FALSE)))
  res<-list()
  for(i in 1:length(old_new_combined)){
    if(names_match[[i]] %in% names(new_burn)){
      ind_old<-length(old_new_combined[[i]])/2  # get the index of the last of the old analyses in the old_new_combined object
      old_only<-old_new_combined[[i]][1:ind_old]  # get the first half of the analyese - these are the old ones
      new_only<-old_new_combined[[i]][(ind_old+1):length(old_new_combined[[1]])]  # Get the second half of the anlyses - these are the new ones
      old_start_tree<-(old_burn[names_match[[i]]][[1]]*10)+1 ## Get the number of the first post-burnin tree - mulitply by 10 because all chains are of length 1000 here - (divide by 100 and then multiply by 1,000 = multiply by 10) - add 1 to get to post-burnin
      new_start_tree<-(new_burn[names_match[[i]]][[1]]*10)+1 ## Get the number of the first post-burnin tree - mulitply by 10 because all chains are of length 1000 here - (divide by 100 and then multiply by 1,000 = multiply by 10) - add 1 to get to post-burnin
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
    }else{
      res[[i]]<-list(NA, "New burn not found")
    }
  }  
  return(res)
}  


if(analysis_set=="NewHeat"){
  ## Get normalized RF distances for pre and post new heating
  rf_comparisons_New_heat<-compare_rf(reanalyses = chains_new_heat, originals = chains_old, new_burn = NewHeat_burns, old_burn = example_chains_burn_LARGE)
  setwd(out_NewHeat)
  save(rf_comparisons_New_heat, file=paste("rf_comp_New_heat_", batch, ".RData", sep=""))
  save(failed_read, file=paste("rf_New_heat_FAILED", batch, ".RData", sep=""))
}

if(analysis_set=="mixed"){
  ## RF distances pre and post nst_mixed
  rf_comparisons_NstMixed<-compare_rf(reanalyses = chains_nst_mixed, originals = chains_old, new_burn = nstMixed_burns, old_burn = example_chains_burn_LARGE)
  setwd(out_nstMixed)
  save(rf_comparisons_NstMixed, file=paste("rf_comp_NstMixed_", batch, ".RData", sep=""))
  save(failed_read, file=paste("rf_NstMixed_FAILED", batch, ".RData", sep=""))
}

if(analysis_set=="G"){
  ### RF distances comparing IG and G
  rf_comparisons_IG<-compare_rf(reanalyses = chains_G_only, originals = chains_old, new_burn = G_only_burns, old_burn = IG_ex_chains_burn_LARGE)
  setwd(out_IG)
  save(rf_comparisons_IG, file=paste("rf_comp_IG_", batch, ".RData", sep=""))
  save(failed_read, file=paste("rf_IG_FAILED", batch, ".RData", sep=""))
}


