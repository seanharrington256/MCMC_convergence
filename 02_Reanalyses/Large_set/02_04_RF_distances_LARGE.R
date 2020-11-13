## Script to compare the tree topologies for the datasets were reanalyzed with different heating or with nst=mixed or with IG switched to just G
##   the script saves an Rdata file with  objects containing RF distances that can then be further analyzed using
##   script 02_05_Summarize_RF_LARGE.R


## Requires that the scripts "01_R_functions_calculate_diagnostics_V2.R" and "02_Wrapper_func_calc_diags.R" are
##   in the directory that this script is run from

## Like script 03_Calc_diags.R, this script is designed to be run on a cluster to split the analysis of thousands of chains up into
##     many smaller jobs that can be run in parallel, the batch_size object determines how many chains are
##     analyzed in a single sequential set - it is set to 10, but can be set to whatever you want
##  note also that the number of samples of each chain that are read in and used to calculate convergence
##    diagnostics can be adjusted by changing the values of trim, target_samp, and max_gens

## it is therefore set up to run from the command line. 2 command line arguments are required to be passed into R
##    the first determines which batch will be run - if the batch size is set at 10, then batch 1 will 
##    be the first 10 analyses, batch 2 will be 11-20, etc.
##   The second argument specifies which set of analyses are being compared - this is really only to specify the directories that contain the relevant files
##  Usage of this script from the command line: >Rscript 02_04_RF_distances_LARGE.R <batch_number> <analysis_set>

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
###### **** Note that these directories were how we had things stored - change as needed
### storage_old_heat_nst and storage_old_IG are the directories containing the original chains for analyses that change nst=mixed/new heating or IG to G respectively
###   storage_new_heat, storage_nst_mixed, and storage_G_only all contain the new, reanalyzed chains to compare topologies to the old ones

main_dir<-getwd()
storage_old_heat_nst<-paste0(main_dir, "/all_original_analyses")
storage_old_IG<-paste0(main_dir, "/IG_original_analyses")
storage_new_heat<-paste0(main_dir, "/Nex_new_heat_chains_READY")
storage_nst_mixed<-paste0(main_dir, "/Nex_nstMixed_chains_READY")
storage_G_only<-paste0(main_dir, "/Nex_IG_chains_READY")
# Out directories - these are simply where output RData filw will end up
out_IG<-paste0(main_dir, "/IG_rf_out")
out_NewHeat<-paste0(main_dir, "/NewHeat_rf_out")
out_nstMixed<-paste0(main_dir, "/nstMixed_rf_out")

setwd(main_dir)

trim<-1    ## Trim is the nth tree to be kept for thinning out MCMC - do this to prevent overly large lists of trees when necessary
target_samp<-1000  # we want to retain 1000 samples
max_gens<-10000000  # we want all chains to be read in up to the 10,000,000th generation so that all are the same


source("01_R_functions_calculate_diagnostics_V2.R")


if(analysis_set=="NewHeat"){
  # Read in the new chains new heat
  setwd(storage_new_heat)
  folders_new_heat_all<-as.list(list.files(pattern=".tar.gz"))
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
if(analysis_set %in% c("mixed", "NewHeat")){
  setwd(storage_old_heat_nst) # set working directory
}
if(analysis_set=="G"){
  setwd(storage_old_IG) # set working directory
}
old_to_read<-list.files()[which(gsub("_READY", "", list.files()) %in% unique(folders))]

## Read in the old chains
chains_old<-lapply(old_to_read, read.comp.catch, trim=trim, target_samp=target_samp, max_gens=max_gens) # read in the chains


setwd(main_dir)


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
compare_rf<-function(reanalyses, originals, out_file, new_burn, old_burn){   
  ## Match names of reanalyses to original chains
  names(reanalyses)<-sapply(reanalyses, function(x) gsub("\\.run.*$|\\.nex.*$|_processed_mbReady.*$", "",names(x)[[1]])) # get the names of the reanalyses
  names(originals)<-sapply(originals, function(x) gsub("\\.run.*$|\\.nex.*$|_processed_mbReady.*$", "",names(x)[[1]])) # get the names of the analyses
  names_match<-names(reanalyses)[which(names(reanalyses) %in% names(originals))] # get the names that match up - analyses may have failed to reanalyze at some point
  # combine chains
  old_new_combined<-lapply(names_match, function(x) c(unlist(originals[x], recursive=FALSE), unlist(reanalyses[x], recursive=FALSE)))
  res<-list()
  for(i in 1:length(old_new_combined)){
      ind_old<-length(old_new_combined[[i]])/2  # get the index of the last of the old analyses in the old_new_combined object
      old_only<-old_new_combined[[i]][1:ind_old]  # get the first half of the analyese - these are the old ones
      new_only<-old_new_combined[[i]][(ind_old+1):length(old_new_combined[[1]])]  # Get the second half of the anlyses - these are the new ones
      old_start_tree<-(old_burn*10)+1 ## Get the number of the first post-burnin tree - mulitply by 10 because all chains are of length 1000 here - (divide by 100 and then multiply by 1,000 = multiply by 10) - add 1 to get to post-burnin
      new_start_tree<-(new_burn*10)+1 ## Get the number of the first post-burnin tree - mulitply by 10 because all chains are of length 1000 here - (divide by 100 and then multiply by 1,000 = multiply by 10) - add 1 to get to post-burnin
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


if(analysis_set=="NewHeat"){
  ## Get normalized RF distances for pre and post new heating
  rf_comparisons_New_heat<-compare_rf(reanalyses = chains_new_heat, originals = chains_old, new_burn = 25, old_burn = 25)
  setwd(out_NewHeat)
  save(rf_comparisons_New_heat, file=paste("rf_comp_New_heat_", batch, ".RData", sep=""))
  save(failed_read, file=paste("rf_New_heat_FAILED", batch, ".RData", sep=""))
}

if(analysis_set=="mixed"){
  ## RF distances pre and post nst_mixed
  rf_comparisons_NstMixed<-compare_rf(reanalyses = chains_nst_mixed, originals = chains_old, new_burn = 25, old_burn = 25)
  setwd(out_nstMixed)
  save(rf_comparisons_NstMixed, file=paste("rf_comp_NstMixed_", batch, ".RData", sep=""))
  save(failed_read, file=paste("rf_NstMixed_FAILED", batch, ".RData", sep=""))
}

if(analysis_set=="G"){
  ### RF distances comparing IG and G
  rf_comparisons_IG<-compare_rf(reanalyses = chains_G_only, originals = chains_old, new_burn = 25, old_burn = 25)
  setwd(out_IG)
  save(rf_comparisons_IG, file=paste("rf_comp_IG_", batch, ".RData", sep=""))
  save(failed_read, file=paste("rf_IG_FAILED", batch, ".RData", sep=""))
}


