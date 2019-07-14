# Load up the packages
library(ggplot2)
library(RMySQL)
library(plyr)
library(phangorn)
library(rwty)

## Set up the locations of the diagnostic summaries of each of the datasets
storage<-"~/Manual_Reanalyses/02_Diagnostics"
# Set up the location of the original chains and RData files with the original burnin and diagnostics
orig_storage<-"~/Manual_Reanalyses/00_Orig_chains"
# Set up a directory to write output to
write_to<-"~/Manual_Reanalyses/04_Plots_etc"
setwd(write_to)

### Read in the original diags for comparison below
setwd(orig_storage)
load("orig_manual_diags.RData")

## Bring in summary
  ## set the working directory to where summaries of the chains are stored
  setwd(storage)
  load("Long_summary_1.RData")
  long_out<-summ

  ## Now that analyses are all laoded up, we can set the working directory to where we want to output figures, etc.
  setwd(write_to)
  
  # Pull out the analyses that failed
  long_failed_summs<-long_out[which(lapply(long_out, function(x) x[[1]])=="FAILED")]  # Get the analyses that failed
  
  # Then get the analyses that were successfully summarized
  diagnostics<-long_out[which(!lapply(long_out, function(x) x[[1]])=="FAILED")]
  ## Name the diagnostics by analysis
  analysis_names<-lapply(diagnostics, function(x) gsub("\\.nex|\\.run1.t.*", "", names(x$Burnin[1])))
  # Set these as the names for the elements of the diagnostics list
  names(diagnostics)<-analysis_names
  
  ## If any analyses are duplicated, remove them
  dup_analysis<-analysis_names[which(duplicated(analysis_names))]
  if(length(dup_analysis)>0){
    diagnostics<-diagnostics[!names(diagnostics)==dup_analysis]
  }
  

# For each analysis, get out a list of which chains passed and failed the diagnostic tests
pass_fail<-lapply(diagnostics, function(x) x[[1]])
pass_fail_combined<-do.call(rbind, pass_fail) #Combine these into a single matrix where each row is a chain
rownames(pass_fail_combined)<-unlist(lapply(diagnostics, function(x) names(x$Burnin)))

## Get the proportion of anlyses that pass each threshold for convergence
prop_pass<-apply(pass_fail_combined, 2, function(x) length(which(x==TRUE))/length(x))
prop_pass

# Split things out based on which set of reanalyses they were in - remembering that one analysis is in more than one category
setwd(storage) # set working directory to where the summaries are - this also has csv files saying which analyses are in each category

####### Identify how many diagnostics are failed for each chain, how many 
  ### parameters fail ESS for each chain
  # Identify parameters that failed ESS in each chain & number
  failed_ESS_pre<-unlist(lapply(diagnostics, function(x) x$Params_ESS_less_200), recursive=FALSE)
  failed_ESS<-lapply(failed_ESS_pre, names)   # get the actual names of these parameters as the value in the object
  failed_ESS_num<-sapply(failed_ESS, length)
  # Identify how many diags are failed for each chain
  pass_fail_reduced<-pass_fail_combined[,!colnames(pass_fail_combined) %in% c("Pass topo pseudo ESS", "Pass Approx topo ESS", "Pass H&W", "Pass Geweke")]  #cut out redundant diags (pseudo and approx topo ESS) and H&W
  num_failed_diags<-apply(pass_fail_reduced, 1, function(x) length(which(x==FALSE)))
  

  # Pass/fail diags
  pass_fail_orig<-orig_manual_diags[[1]] #pass/fail table for original diagnostics
  pass_fail_reduced_orig<-pass_fail_orig[,!colnames(pass_fail_orig) %in% c("Pass topo pseudo ESS", "Pass Approx topo ESS", "Pass H&W", "Pass Geweke")]  #cut out redundant diags (pseudo and approx topo ESS) and H&W
  num_failed_diags_orig<-apply(pass_fail_reduced_orig, 1, function(x) length(which(x==FALSE)))
  ## Pars fail ESS
  failed_ESS_orig<-orig_manual_diags[[3]] #pars fail ESS for original diags
  failed_ESS_num_orig<-sapply(failed_ESS_orig, length)# get numnber of pars that originally failed ESS
  ## Raw diags
  raw_diags_orig<-orig_manual_diags[[4]] #all raw original diags
  
## Let's compare these numbers between original and reestimated analyses
## Start with the overall number of failed diagnostics
  # first, need to sum across the multiple chains in a single analysis to compare whole analyses, because chain 1 or 2 has no real meaning
  analysis_base_names<-unique(gsub("\\.nex.*$|\\.run.*$", "", names(num_failed_diags)))
  analysis_base_names<-gsub("_processed_mbReady", "", analysis_base_names)  #get rid of any _processed_mbReady still on these names
  fail_diags_per_ana<-lapply(analysis_base_names, function(x) sum(num_failed_diags[grep(x, names(num_failed_diags))])) # Use an lapply to go through this list of names and get the total number of diags failed for these chains
  names(fail_diags_per_ana)<-analysis_base_names
  # Then, let's do the same for the original diagnostics
  fail_diags_per_ana_orig<-lapply(analysis_base_names, function(x) sum(num_failed_diags_orig [grep(x, names(num_failed_diags_orig))])) # Use an lapply to go through this list of names and get the total number of diags failed for these chains
  names(fail_diags_per_ana_orig)<-analysis_base_names
  
  # Let's subtract these and see how it shakes out
  diff_num_failed_diags<-unlist(fail_diags_per_ana_orig)-unlist(fail_diags_per_ana)
  mean(diff_num_failed_diags)
  median(diff_num_failed_diags)
  sd(diff_num_failed_diags)
  hist(diff_num_failed_diags)

  ## Let's compare how many parameters fail the ESS between original and reanalyses - need to compare the proportion of parameters failing, because down lower for HKY and empirical base frequency reanalyses, the number of paramters is lower in the reanalyses
  fail_ESS_per_ana<-lapply(analysis_base_names, function(x) sum(failed_ESS_num[grep(x, names(failed_ESS_num))]))
  names(fail_ESS_per_ana)<-analysis_base_names
  fail_ESS_per_ana_orig<-lapply(analysis_base_names, function(x) sum(failed_ESS_num_orig[grep(x, names(failed_ESS_num_orig))]))
  names(fail_ESS_per_ana_orig)<-analysis_base_names
  num_params<-lapply(diagnostics, function(x) 2*(length(x$ESS_values[[1]])-1))  ## Get the total number of parameters that ESS is calculated for (subtract 1 because topo us not used in calculation of number of params that fail ESS) - multiply by 2 because we have 2 chains
  ratio_ESS_failed<-unlist(fail_ESS_per_ana)/unlist(num_params) # ratio of failed ESS to total # pars
  ## Get the ratio for the original analyses
  num_params_orig<-lapply(raw_diags_orig[names(fail_ESS_per_ana_orig)], function(x) 2*(length(x$ESS_values[[1]])-1))
  ratio_ESS_failed_orig<-unlist(fail_ESS_per_ana_orig)/unlist(num_params_orig)
  # Difference in the ratios
  diff_num_failed_ESS<-unlist(ratio_ESS_failed_orig)-unlist(ratio_ESS_failed)
  mean(diff_num_failed_ESS)
  median(diff_num_failed_ESS)
  sd(diff_num_failed_ESS)
  hist(diff_num_failed_ESS)
  

  

## combine these all into a table
  diffs_reanalysis_long<-cbind(analysis_base_names, diff_num_failed_diags, diff_num_failed_ESS)  ## Add empty placeholder column that will be replaced with which subset each analysis belongs to
  
## Write out a csv with this table in it
  setwd(write_to)
  write.csv(diffs_reanalysis_long, file="Diffs_Run_longer.csv")
  
## Get the mean, median, and standard deviation for each of the subsets
  # Mean
  summ_diffs_long<-apply(diffs_reanalysis_long[, 2:3], 2, function(x) mean(as.numeric(x)))
  names(summ_diffs_long)<-c("Failed diags", "failed ESS")
  
  # Median
  summ_diffs_med_long<-apply(diffs_reanalysis_long[, 2:3], 2, function(x) median(as.numeric(x)))
  names(summ_diffs_med_long)<-c("Failed diags", "failed ESS")
  
  # Std dev
  summ_diffs_sd_long<-apply(diffs_reanalysis_long[, 2:3], 2, function(x) sd(as.numeric(x)))
  names(summ_diffs_sd_long)<-c("Failed diags", "failed ESS")
  
  summ_table_manual_long<-rbind(summ_diffs_long, summ_diffs_med_long, summ_diffs_sd_long)
  rownames(summ_table_manual_long)<-c("mean", "median", "sd")
  
  write.csv(summ_table_manual_long, file="Diffs_summ_manual_Run_longer.csv")
  
  
# Get the burnin for each chain for later
  Run_longer_burns<-lapply(diagnostics, function(x) max(unlist(x$Burnin)))
  names(Run_longer_burns)<-lapply(diagnostics, function(x) gsub("\\.run.*$|\\.nex.*$|_processed_mbReady.*$", "", names(x$Burnin)[[1]]))
  save(Run_longer_burns, file="Run_longer_burnin.RData")
  
  
########################################################################################################################
########################################################################################################################
###########  Same as above but for the set of analyses that were run longer and retaining more samples 
########################################################################################################################
########################################################################################################################
  ## Bring in summary
  ## set the working directory to where summaries of the chains are stored
  setwd(storage)
  load("LongMoreSample_summary_1.RData")
  long_more_samp_out<-summ
  
  ## Now that analyses are all laoded up, we can set the working directory to where we want to output figures, etc.
  setwd(write_to)
  
  # Pull out the analyses that failed
  long_more_samp_failed_summs<-long_more_samp_out[which(lapply(long_more_samp_out, function(x) x[[1]])=="FAILED")]  # Get the analyses that failed
  
  
  # Then get the analyses that were successfully summarized
  diagnostics<-long_more_samp_out[which(!lapply(long_more_samp_out, function(x) x[[1]])=="FAILED")]
  ## Name the diagnostics by analysis
  analysis_names<-lapply(diagnostics, function(x) gsub("\\.nex|\\.run1.t.*", "", names(x$Burnin[1])))
  # Set these as the names for the elements of the diagnostics list
  names(diagnostics)<-analysis_names
  
  ## If any analyses are duplicated, remove them
  dup_analysis<-analysis_names[which(duplicated(analysis_names))]
  if(length(dup_analysis)>0){
    diagnostics<-diagnostics[!names(diagnostics)==dup_analysis]
  }
  
  
  # For each analysis, get out a list of which chains passed and failed the diagnostic tests
  pass_fail<-lapply(diagnostics, function(x) x[[1]])
  pass_fail_combined<-do.call(rbind, pass_fail) #Combine these into a single matrix where each row is a chain
  rownames(pass_fail_combined)<-unlist(lapply(diagnostics, function(x) names(x$Burnin)))
  
  ## Get the proportion of anlyses that pass each threshold for convergence
  prop_pass<-apply(pass_fail_combined, 2, function(x) length(which(x==TRUE))/length(x))
  prop_pass
  
  # Split things out based on which set of reanalyses they were in - remembering that one analysis is in more than one category
  setwd(storage) # set working directory to where the summaries are - this also has csv files saying which analyses are in each category
  
  ####### Identify how many diagnostics are failed for each chain, how many 
  ### parameters fail ESS for each chain
  # Identify parameters that failed ESS in each chain & number
  failed_ESS_pre<-unlist(lapply(diagnostics, function(x) x$Params_ESS_less_200), recursive=FALSE)
  failed_ESS<-lapply(failed_ESS_pre, names)   # get the actual names of these parameters as the value in the object
  failed_ESS_num<-sapply(failed_ESS, length)
  # Identify how many diags are failed for each chain
  pass_fail_reduced<-pass_fail_combined[,!colnames(pass_fail_combined) %in% c("Pass topo pseudo ESS", "Pass Approx topo ESS", "Pass H&W", "Pass Geweke")]  #cut out redundant diags (pseudo and approx topo ESS) and H&W
  num_failed_diags<-apply(pass_fail_reduced, 1, function(x) length(which(x==FALSE)))
  
  ## Compare numbers between original and reestimated analyses
  ## Start with the overall number of failed diagnostics
  # first, need to sum across the multiple chains in a single analysis to compare whole analyses, because chain 1 or 2 has no real meaning
  analysis_base_names<-unique(gsub("\\.nex.*$|\\.run.*$", "", names(num_failed_diags)))
  analysis_base_names<-gsub("_processed_mbReady", "", analysis_base_names)  #get rid of any _processed_mbReady still on these names
  fail_diags_per_ana<-lapply(analysis_base_names, function(x) sum(num_failed_diags[grep(x, names(num_failed_diags))])) # Use an lapply to go through this list of names and get the total number of diags failed for these chains
  names(fail_diags_per_ana)<-analysis_base_names
  # Then, let's do the same for the original diagnostics
  fail_diags_per_ana_orig<-lapply(analysis_base_names, function(x) sum(num_failed_diags_orig [grep(x, names(num_failed_diags_orig))])) # Use an lapply to go through this list of names and get the total number of diags failed for these chains
  names(fail_diags_per_ana_orig)<-analysis_base_names
  
  # Let's subtract these and see how it shakes out
  diff_num_failed_diags<-unlist(fail_diags_per_ana_orig)-unlist(fail_diags_per_ana)
  mean(diff_num_failed_diags)
  median(diff_num_failed_diags)
  sd(diff_num_failed_diags)
  hist(diff_num_failed_diags)
  
  ## Let's compare how many parameters fail the ESS between original and reanalyses - need to compare the proportion of parameters failing, because down lower for HKY and empirical base frequency reanalyses, the number of paramters is lower in the reanalyses
  fail_ESS_per_ana<-lapply(analysis_base_names, function(x) sum(failed_ESS_num[grep(x, names(failed_ESS_num))]))
  names(fail_ESS_per_ana)<-analysis_base_names
  fail_ESS_per_ana_orig<-lapply(analysis_base_names, function(x) sum(failed_ESS_num_orig[grep(x, names(failed_ESS_num_orig))]))
  names(fail_ESS_per_ana_orig)<-analysis_base_names
  num_params<-lapply(diagnostics, function(x) 2*(length(x$ESS_values[[1]])-1))  ## Get the total number of parameters that ESS is calculated for (subtract 1 because topo us not used in calculation of number of params that fail ESS) - multiply by 2 because we have 2 chains
  ratio_ESS_failed<-unlist(fail_ESS_per_ana)/unlist(num_params) # ratio of failed ESS to total # pars
  ## Get the ratio for the original analyses
  num_params_orig<-lapply(raw_diags_orig[names(fail_ESS_per_ana_orig)], function(x) 2*(length(x$ESS_values[[1]])-1))
  ratio_ESS_failed_orig<-unlist(fail_ESS_per_ana_orig)/unlist(num_params_orig)
  # Difference in the ratios
  diff_num_failed_ESS<-unlist(ratio_ESS_failed_orig)-unlist(ratio_ESS_failed)
  mean(diff_num_failed_ESS)
  median(diff_num_failed_ESS)
  sd(diff_num_failed_ESS)
  hist(diff_num_failed_ESS)
  

  
  
  ## combine these all into a table
  diffs_reanalysis_long_more<-cbind(analysis_base_names, diff_num_failed_diags, diff_num_failed_ESS)  ## Add empty placeholder column that will be replaced with which subset each analysis belongs to
  
  ## Write out a csv with this table in it
  setwd(write_to)
  write.csv(diffs_reanalysis_long_more, file="Diffs_Run_longerMoreSamples.csv")
  
  ## Get the mean, median, and standard deviation for each of the subsets
  # Mean
  summ_diffs_long_more<-apply(diffs_reanalysis_long_more[, 2:3], 2, function(x) mean(as.numeric(x)))
  names(summ_diffs_long_more)<-c("Failed diags", "failed ESS")
  
  # Median
  summ_diffs_med_long_more<-apply(diffs_reanalysis_long_more[, 2:3], 2, function(x) median(as.numeric(x)))
  names(summ_diffs_med_long_more)<-c("Failed diags", "failed ESS")
  
  # Std dev
  summ_diffs_sd_long_more<-apply(diffs_reanalysis_long_more[, 2:3], 2, function(x) sd(as.numeric(x)))
  names(summ_diffs_sd_long_more)<-c("Failed diags", "failed ESS")
  
  summ_table_manual_long_more<-rbind(summ_diffs_long_more, summ_diffs_med_long_more, summ_diffs_sd_long_more)
  rownames(summ_table_manual_long_more)<-c("mean", "median", "sd")
  
  write.csv(summ_table_manual_long_more, file="Diffs_summ_manual_Run_longerMoreSample.csv")
  
  
  # Get the burnin for each chain for later
  Run_longerMore_burns<-lapply(diagnostics, function(x) max(unlist(x$Burnin)))
  names(Run_longerMore_burns)<-lapply(diagnostics, function(x) gsub("\\.run.*$|\\.nex.*$|_processed_mbReady.*$", "", names(x$Burnin)[[1]]))
  save(Run_longerMore_burns, file="Run_longerMoreSample_burnin.RData")
  
  
########################################################################################################################
########################################################################################################################
###########  Same as above but for the set of HKY analyses
########################################################################################################################
########################################################################################################################
  ## Bring in summary
  ## set the working directory to where summaries of the chains are stored
  setwd(storage)
  load("HKY_summary_1.RData")
  HKY_out<-summ
  
  ## Now that analyses are all laoded up, we can set the working directory to where we want to output figures, etc.
  setwd(write_to)
  
  # Pull out the analyses that failed
  HKY_failed_summs<-HKY_out[which(lapply(HKY_out, function(x) x[[1]])=="FAILED")]  # Get the analyses that failed
  
  # Then get the analyses that were successfully summarized
  diagnostics<-HKY_out[which(!lapply(HKY_out, function(x) x[[1]])=="FAILED")]
  ## Name the diagnostics by analysis
  analysis_names<-lapply(diagnostics, function(x) gsub("\\.nex|\\.run1.t.*", "", names(x$Burnin[1])))
  # Set these as the names for the elements of the diagnostics list
  names(diagnostics)<-analysis_names
  
  ## If any analyses are duplicated, remove them
  dup_analysis<-analysis_names[which(duplicated(analysis_names))]
  if(length(dup_analysis)>0){
    diagnostics<-diagnostics[!names(diagnostics)==dup_analysis]
  }
  
  
  # For each analysis, get out a list of which chains passed and failed the diagnostic tests
  pass_fail<-lapply(diagnostics, function(x) x[[1]])
  pass_fail_combined<-do.call(rbind, pass_fail) #Combine these into a single matrix where each row is a chain
  rownames(pass_fail_combined)<-unlist(lapply(diagnostics, function(x) names(x$Burnin)))
  
  ## Get the proportion of anlyses that pass each threshold for convergence
  prop_pass<-apply(pass_fail_combined, 2, function(x) length(which(x==TRUE))/length(x))
  prop_pass
  

  ####### Identify how many diagnostics are failed for each chain, how many 
  ### parameters fail ESS for each chain
  # Identify parameters that failed ESS in each chain & number
  failed_ESS_pre<-unlist(lapply(diagnostics, function(x) x$Params_ESS_less_200), recursive=FALSE)
  failed_ESS<-lapply(failed_ESS_pre, names)   # get the actual names of these parameters as the value in the object
  failed_ESS_num<-sapply(failed_ESS, length)
  # Identify how many diags are failed for each chain
  pass_fail_reduced<-pass_fail_combined[,!colnames(pass_fail_combined) %in% c("Pass topo pseudo ESS", "Pass Approx topo ESS", "Pass H&W", "Pass Geweke")]  #cut out redundant diags (pseudo and approx topo ESS) and H&W
  num_failed_diags<-apply(pass_fail_reduced, 1, function(x) length(which(x==FALSE)))
  
  ## Compare numbers between original and reestimated analyses
  ## Start with the overall number of failed diagnostics
  # first, need to sum across the multiple chains in a single analysis to compare whole analyses, because chain 1 or 2 has no real meaning
  analysis_base_names<-unique(gsub("\\.nex.*$|\\.run.*$", "", names(num_failed_diags)))
  analysis_base_names<-gsub("_processed_mbReady", "", analysis_base_names)  #get rid of any _processed_mbReady still on these names
  fail_diags_per_ana<-lapply(analysis_base_names, function(x) sum(num_failed_diags[grep(x, names(num_failed_diags))])) # Use an lapply to go through this list of names and get the total number of diags failed for these chains
  names(fail_diags_per_ana)<-analysis_base_names
  # Then, let's do the same for the original diagnostics
  fail_diags_per_ana_orig<-lapply(analysis_base_names, function(x) sum(num_failed_diags_orig [grep(x, names(num_failed_diags_orig))])) # Use an lapply to go through this list of names and get the total number of diags failed for these chains
  names(fail_diags_per_ana_orig)<-analysis_base_names
  
  # Let's subtract these and see how it shakes out
  diff_num_failed_diags<-unlist(fail_diags_per_ana_orig)-unlist(fail_diags_per_ana)
  mean(diff_num_failed_diags)
  median(diff_num_failed_diags)
  sd(diff_num_failed_diags)
  hist(diff_num_failed_diags)
  
  ## Let's compare how many parameters fail the ESS between original and reanalyses - need to compare the proportion of parameters failing, because down lower for HKY and empirical base frequency reanalyses, the number of paramters is lower in the reanalyses
  fail_ESS_per_ana<-lapply(analysis_base_names, function(x) sum(failed_ESS_num[grep(x, names(failed_ESS_num))]))
  names(fail_ESS_per_ana)<-analysis_base_names
  fail_ESS_per_ana_orig<-lapply(analysis_base_names, function(x) sum(failed_ESS_num_orig[grep(x, names(failed_ESS_num_orig))]))
  names(fail_ESS_per_ana_orig)<-analysis_base_names
  num_params<-lapply(diagnostics, function(x) 2*(length(x$ESS_values[[1]])-1))  ## Get the total number of parameters that ESS is calculated for (subtract 1 because topo us not used in calculation of number of params that fail ESS) - multiply by 2 because we have 2 chains
  ratio_ESS_failed<-unlist(fail_ESS_per_ana)/unlist(num_params) # ratio of failed ESS to total # pars
  ## Get the ratio for the original analyses
  num_params_orig<-lapply(raw_diags_orig[names(fail_ESS_per_ana_orig)], function(x) 2*(length(x$ESS_values[[1]])-1))
  ratio_ESS_failed_orig<-unlist(fail_ESS_per_ana_orig)/unlist(num_params_orig)
  # Difference in the ratios
  diff_num_failed_ESS<-unlist(ratio_ESS_failed_orig)-unlist(ratio_ESS_failed)
  mean(diff_num_failed_ESS)
  median(diff_num_failed_ESS)
  sd(diff_num_failed_ESS)
  hist(diff_num_failed_ESS)

  
  ## combine these all into a table
  diffs_reanalysis_hky<-cbind(analysis_base_names, diff_num_failed_diags, diff_num_failed_ESS)  ## Add empty placeholder column that will be replaced with which subset each analysis belongs to
  
  ## Write out a csv with this table in it
  setwd(write_to)
  write.csv(diffs_reanalysis_hky, file="Diffs_HKY.csv")
  
  ## Get the mean, median, and standard deviation for each of the subsets
  # Mean
  summ_diffs_hky<-apply(diffs_reanalysis_hky[, 2:3], 2, function(x) mean(as.numeric(x)))
  names(summ_diffs_hky)<-c("Failed diags", "failed ESS")
  
  # Median
  summ_diffs_med_hky<-apply(diffs_reanalysis_hky[, 2:3], 2, function(x) median(as.numeric(x)))
  names(summ_diffs_med_hky)<-c("Failed diags", "failed ESS")
  
  # Std dev
  summ_diffs_sd_hky<-apply(diffs_reanalysis_hky[, 2:3], 2, function(x) sd(as.numeric(x)))
  names(summ_diffs_sd_hky)<-c("Failed diags", "failed ESS")
  
  summ_table_manual_hky<-rbind(summ_diffs_hky, summ_diffs_med_hky, summ_diffs_sd_hky)
  rownames(summ_table_manual_hky)<-c("mean", "median", "sd")
  
  write.csv(summ_table_manual_hky, file="Diffs_summ_manual_HKY.csv")
  
  
  # Get the burnin for each chain so that we can plot things in RWTY later
  hky_burns<-lapply(diagnostics, function(x) max(unlist(x$Burnin)))
  names(hky_burns)<-lapply(diagnostics, function(x) gsub("\\.run.*$|\\.nex.*$|_processed_mbReady.*$", "", names(x$Burnin)[[1]]))
  save(hky_burns, file="HKY_burnin.RData")
  
  
  
########################################################################################################################
########################################################################################################################
###########  Same as above but for the set of empirical base frequency analyses
########################################################################################################################
########################################################################################################################
  ## Bring in summary
  ## set the working directory to where summaries of the chains are stored
  setwd(storage)
  load("EmpBF_summary_1.RData")
  empBF_out<-summ
  
  ## Now that analyses are all laoded up, we can set the working directory to where we want to output figures, etc.
  setwd(write_to)
  
  # Pull out the analyses that failed
  empBF_failed_summs<-empBF_out[which(lapply(empBF_out, function(x) x[[1]])=="FAILED")]  # Get the analyses that failed
  
  # Then get the analyses that were successfully summarized
  diagnostics<-empBF_out[which(!lapply(empBF_out, function(x) x[[1]])=="FAILED")]
  ## Name the diagnostics by analysis
  analysis_names<-lapply(diagnostics, function(x) gsub("\\.nex|\\.run1.t.*", "", names(x$Burnin[1])))
  # Set these as the names for the elements of the diagnostics list
  names(diagnostics)<-analysis_names
  
  ## If any analyses are duplicated, remove them
  dup_analysis<-analysis_names[which(duplicated(analysis_names))]
  if(length(dup_analysis)>0){
    diagnostics<-diagnostics[!names(diagnostics)==dup_analysis]
  }
  
  
  # For each analysis, get out a list of which chains passed and failed the diagnostic tests
  pass_fail<-lapply(diagnostics, function(x) x[[1]])
  pass_fail_combined<-do.call(rbind, pass_fail) #Combine these into a single matrix where each row is a chain
  rownames(pass_fail_combined)<-unlist(lapply(diagnostics, function(x) names(x$Burnin)))
  
  ## Get the proportion of anlyses that pass each threshold for convergence
  prop_pass<-apply(pass_fail_combined, 2, function(x) length(which(x==TRUE))/length(x))
  prop_pass
  
  # Split things out based on which set of reanalyses they were in - remembering that one analysis is in more than one category
  setwd(storage) # set working directory to where the summaries are - this also has csv files saying which analyses are in each category
  
  ####### Identify how many diagnostics are failed for each chain, how many 
  ### parameters fail ESS for each chain
  # Identify parameters that failed ESS in each chain & number
  failed_ESS_pre<-unlist(lapply(diagnostics, function(x) x$Params_ESS_less_200), recursive=FALSE)
  failed_ESS<-lapply(failed_ESS_pre, names)   # get the actual names of these parameters as the value in the object
  failed_ESS_num<-sapply(failed_ESS, length)
  # Identify how many diags are failed for each chain
  pass_fail_reduced<-pass_fail_combined[,!colnames(pass_fail_combined) %in% c("Pass topo pseudo ESS", "Pass Approx topo ESS", "Pass H&W", "Pass Geweke")]  #cut out redundant diags (pseudo and approx topo ESS) and H&W
  num_failed_diags<-apply(pass_fail_reduced, 1, function(x) length(which(x==FALSE)))
  
  ## Compare numbers between original and reestimated analyses
  ## Start with the overall number of failed diagnostics
  # first, need to sum across the multiple chains in a single analysis to compare whole analyses, because chain 1 or 2 has no real meaning
  analysis_base_names<-unique(gsub("\\.nex.*$|\\.run.*$", "", names(num_failed_diags)))
  analysis_base_names<-gsub("_processed_mbReady", "", analysis_base_names)  #get rid of any _processed_mbReady still on these names
  fail_diags_per_ana<-lapply(analysis_base_names, function(x) sum(num_failed_diags[grep(x, names(num_failed_diags))])) # Use an lapply to go through this list of names and get the total number of diags failed for these chains
  names(fail_diags_per_ana)<-analysis_base_names
  # Then, let's do the same for the original diagnostics
  fail_diags_per_ana_orig<-lapply(analysis_base_names, function(x) sum(num_failed_diags_orig [grep(x, names(num_failed_diags_orig))])) # Use an lapply to go through this list of names and get the total number of diags failed for these chains
  names(fail_diags_per_ana_orig)<-analysis_base_names
  
  # Let's subtract these and see how it shakes out
  diff_num_failed_diags<-unlist(fail_diags_per_ana_orig)-unlist(fail_diags_per_ana)
  mean(diff_num_failed_diags)
  median(diff_num_failed_diags)
  sd(diff_num_failed_diags)
  hist(diff_num_failed_diags)
  
  ## Let's compare how many parameters fail the ESS between original and reanalyses - need to compare the proportion of parameters failing, because down lower for HKY and empirical base frequency reanalyses, the number of paramters is lower in the reanalyses
  fail_ESS_per_ana<-lapply(analysis_base_names, function(x) sum(failed_ESS_num[grep(x, names(failed_ESS_num))]))
  names(fail_ESS_per_ana)<-analysis_base_names
  fail_ESS_per_ana_orig<-lapply(analysis_base_names, function(x) sum(failed_ESS_num_orig[grep(x, names(failed_ESS_num_orig))]))
  names(fail_ESS_per_ana_orig)<-analysis_base_names
  num_params<-lapply(diagnostics, function(x) 2*(length(x$ESS_values[[1]])-1))  ## Get the total number of parameters that ESS is calculated for (subtract 1 because topo us not used in calculation of number of params that fail ESS) - multiply by 2 because we have 2 chains
  ratio_ESS_failed<-unlist(fail_ESS_per_ana)/unlist(num_params) # ratio of failed ESS to total # pars
  ## Get the ratio for the original analyses
  num_params_orig<-lapply(raw_diags_orig[names(fail_ESS_per_ana_orig)], function(x) 2*(length(x$ESS_values[[1]])-1))
  ratio_ESS_failed_orig<-unlist(fail_ESS_per_ana_orig)/unlist(num_params_orig)
  # Difference in the ratios
  diff_num_failed_ESS<-unlist(ratio_ESS_failed_orig)-unlist(ratio_ESS_failed)
  mean(diff_num_failed_ESS)
  median(diff_num_failed_ESS)
  sd(diff_num_failed_ESS)
  hist(diff_num_failed_ESS)
  

  ## combine these all into a table
  diffs_reanalysis_empBF<-cbind(analysis_base_names, diff_num_failed_diags, diff_num_failed_ESS)  ## Add empty placeholder column that will be replaced with which subset each analysis belongs to
  
  ## Write out a csv with this table in it
  setwd(write_to)
  write.csv(diffs_reanalysis_empBF, file="Diffs_empBF.csv")
  
  ## Get the mean, median, and standard deviation for each of the subsets
  # Mean
  summ_diffs_empBF<-apply(diffs_reanalysis_empBF[, 2:3], 2, function(x) mean(as.numeric(x)))
  names(summ_diffs_empBF)<-c("Failed diags", "failed ESS")
  
  # Median
  summ_diffs_med_empBF<-apply(diffs_reanalysis_empBF[, 2:3], 2, function(x) median(as.numeric(x)))
  names(summ_diffs_med_empBF)<-c("Failed diags", "failed ESS")
  
  # Std dev
  summ_diffs_sd_empBF<-apply(diffs_reanalysis_empBF[, 2:3], 2, function(x) sd(as.numeric(x)))
  names(summ_diffs_sd_empBF)<-c("Failed diags", "failed ESS")
  
  summ_table_manual_empBF<-rbind(summ_diffs_empBF, summ_diffs_med_empBF, summ_diffs_sd_empBF)
  rownames(summ_table_manual_empBF)<-c("mean", "median", "sd")
  
  write.csv(summ_table_manual_empBF, file="Diffs_summ_manual_EmpBF.csv")
  
  
  # Get the burnin for each chain so that we can plot things in RWTY later
  emp_BF_burns<-lapply(diagnostics, function(x) max(unlist(x$Burnin)))
  names(emp_BF_burns)<-lapply(diagnostics, function(x) gsub("\\.run.*$|\\.nex.*$|_processed_mbReady.*$", "", names(x$Burnin)[[1]]))
  save(emp_BF_burns, file="EmpBFreq_burnin.RData")
  
  
### Combine all into a single table
  all_diffs_ESS_props<-rbind("Longer", diffs_reanalysis_long, "Longer more samples", diffs_reanalysis_long_more,
        "HKY", diffs_reanalysis_hky, "Emp Base Freq", diffs_reanalysis_empBF)
  
  write.csv(all_diffs_ESS_props, file="Diffs_manual_ESS_props.csv")
  
  
  
  
  
  