## Script to compare teh diagnostics for chains before and after rerunning chains that were manually examined to improve convergence


# Load up the packages
library(ggplot2)
library(plyr)
library(phangorn)
library(rwty)


## Change these directories as necessary
## Set up the locations of the diagnostic summaries of each of the datasets
storage<-"02_Reanalyses/Manual/manual_diags"
# Set up the location of the original chains and RData files with the original diagnostics
orig_storage<-"02_Reanalyses/Manual"
# Set up a directory to write output to
write_to<-"02_Reanalyses/Manual/manual_output"
setwd(write_to)

### Read in the original diags for comparison below
setwd(orig_storage)
load("orig_manual_diags.RData")  # this file was output by script 04_Analyze_diagnostic_output.R

# Pass/fail diags
pass_fail_orig<-orig_manual_diags[[1]] #pass/fail table for original diagnostics
pass_fail_reduced_orig<-pass_fail_orig[,!colnames(pass_fail_orig) %in% c("Pass topo pseudo ESS", "Pass ASDSF 2 chains", "Pass SF corr 2 chains", "Pass PSRF 2 chains")]
num_failed_diags_orig<-apply(pass_fail_reduced_orig, 1, function(x) length(which(x==FALSE)))
## Pars fail ESS
failed_ESS_orig<-orig_manual_diags[[2]] #pars fail ESS for original diags
failed_ESS_num_orig<-sapply(failed_ESS_orig, length)# get numnber of pars that originally failed ESS
## Raw diags
raw_diags_orig<-orig_manual_diags[[3]] #all raw original diags

## Number of diagnostics failed for all original analyses
analysis_base_names_orig<-unique(gsub("\\.nex.*$|\\.run.*$", "", names(num_failed_diags_orig)))
fail_diags_per_ana_orig_all<-lapply(analysis_base_names_orig, function(x) sum(num_failed_diags_orig [grep(x, names(num_failed_diags_orig))])) # Use an lapply to go through this list of names and get the total number of diags failed for these chains
names(fail_diags_per_ana_orig_all)<-analysis_base_names_orig

## Ratio of ESS values failed in original analyses
fail_ESS_per_ana_orig_all<-lapply(analysis_base_names_orig, function(x) sum(failed_ESS_num_orig[grep(x, names(failed_ESS_num_orig))]))
names(fail_ESS_per_ana_orig_all)<-analysis_base_names_orig
num_params_orig_all<-lapply(raw_diags_orig[names(fail_ESS_per_ana_orig_all)], function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))
ratio_ESS_failed_orig_all<-unlist(fail_ESS_per_ana_orig_all)/unlist(num_params_orig_all)

# whole analyses pass or fail asdsf
analyses_pass_asdsf_orig<-sapply(raw_diags_orig, function(x) unique((x[[1]][,"Pass ASDSF"])))




########################################################################################################################
########################################################################################################################
###########  Set of analyses that were run longer and retaining more samples 
########################################################################################################################
########################################################################################################################
## Bring in summary
## set the working directory to where summaries of the chains are stored
setwd(storage)

# Read in all output - need to ensure that don't overwrite each time, as each file has same object names for summ
to_read<-list.files(pattern="summary_RunLonger.*.RData")
out_all_pre<-list()    #make a list to dump the batches into
for(i in 1:length(to_read)){   # Dump the summ object from each batch into an element of a new list
  load(to_read[[i]])
  out_all_pre[[i]]<-summ
}
## Combine everything into a list where each element of the list is a single analysis
out_all<-unlist(out_all_pre, recursive=FALSE)

## Now that analyses are all laoded up, we can set the working directory to where we want to output figures, etc.
setwd(write_to)

# Pull out the analyses that failed
failed_summs<-out_all[which(lapply(out_all, function(x) x[[1]])=="FAILED")]  # Get the analyses that failed

# Then get the analyses that were successfully summarized
diagnosticsAll<-out_all[which(!lapply(out_all, function(x) x[[1]])=="FAILED")]
## Name the diagnostics by analysis
analysis_names<-lapply(diagnosticsAll, function(x) gsub("\\.nex|\\.run1.t.*", "", names(x$"ACT params")[[1]]))
# Set these as the names for the elements of the diagnostics list
names(diagnosticsAll)<-analysis_names

## If any analyses are duplicated, remove them
dup_analysis<-analysis_names[which(duplicated(analysis_names))]
if(length(dup_analysis)>0){
  diagnosticsAll<-diagnosticsAll[!names(diagnosticsAll)==dup_analysis]
}


# For each analysis, get out a list of which chains passed and failed the diagnostic tests
pass_fail<-lapply(diagnosticsAll, function(x) as.data.frame(x[[1]]))
pass_fail_combined<-do.call(rbind.fill, pass_fail) #Combine these into a single matrix where each row is a chain
rownames(pass_fail_combined)<-unlist(lapply(diagnosticsAll, function(x) names(x$"ACT params")))

## Get the proportion of anlyses that pass each threshold for convergence
prop_pass<-apply(pass_fail_combined, 2, function(x) length(which(x==TRUE))/length(x))
prop_pass
  
# whole analyses pass or fail asdsf
analyses_pass_asdsf_runLonger<-sapply(diagnosticsAll, function(x) unique((x[[1]][,"Pass ASDSF"])))



  
  ####### Identify how many diagnostics are failed for each chain, how many 
  ### parameters fail ESS for each chain
  # Identify parameters that failed ESS in each chain & number
  failed_ESS_pre<-unlist(lapply(diagnosticsAll, function(x) x$Params_ESS_less_200), recursive=FALSE)
  failed_ESS<-lapply(failed_ESS_pre, names)   # get the actual names of these parameters as the value in the object
  failed_ESS_num<-sapply(failed_ESS, length)
  # Identify how many diags are failed for each chain
  pass_fail_reduced<-pass_fail_combined[,!colnames(pass_fail_combined) %in% c("Pass topo pseudo ESS", "Pass ASDSF 2 chains", "Pass SF corr 2 chains", "Pass PSRF 2 chains")]  #cut out redundant diags (pseudo and approx topo ESS) and H&W
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
  hist(diff_num_failed_diags, main="RunLonger")
  
  ## Let's compare how many parameters fail the ESS between original and reanalyses - need to compare the proportion of parameters failing, because down lower for HKY and empirical base frequency reanalyses, the number of paramters is lower in the reanalyses
  fail_ESS_per_ana<-lapply(analysis_base_names, function(x) sum(failed_ESS_num[grep(x, names(failed_ESS_num))]))
  names(fail_ESS_per_ana)<-analysis_base_names
  fail_ESS_per_ana_orig<-lapply(analysis_base_names, function(x) sum(failed_ESS_num_orig[grep(x, names(failed_ESS_num_orig))]))
  names(fail_ESS_per_ana_orig)<-analysis_base_names
  num_params<-lapply(diagnosticsAll, function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))  ## Get the total number of parameters that ESS is calculated for (subtract 1 because topo us not used in calculation of number of params that fail ESS) - multiply by 2 because we have 2 chains
  ratio_ESS_failed<-unlist(fail_ESS_per_ana)/unlist(num_params) # ratio of failed ESS to total # pars
  ## Get the ratio for the original analyses
  num_params_orig<-lapply(raw_diags_orig[names(fail_ESS_per_ana_orig)], function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))
  ratio_ESS_failed_orig<-unlist(fail_ESS_per_ana_orig)/unlist(num_params_orig)
  # Difference in the ratios
  diff_num_failed_ESS<-unlist(ratio_ESS_failed_orig)-unlist(ratio_ESS_failed)
  mean(diff_num_failed_ESS)
  median(diff_num_failed_ESS)
  sd(diff_num_failed_ESS)
  hist(diff_num_failed_ESS, main="RunLonger")
  
  prop_pass_runLonger<-prop_pass
  ratio_ESS_failed_runLonger<-ratio_ESS_failed
  fail_diags_per_ana_runLonger<-fail_diags_per_ana
  
  
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
  
  
########################################################################################################################
########################################################################################################################
###########  Same as above but for the set of HKY analyses
########################################################################################################################
########################################################################################################################
  ## Bring in summary
  ## set the working directory to where summaries of the chains are stored
  setwd(storage)
  
  # Read in all output - need to ensure that don't overwrite each time, as each file has same object names for summ
  to_read<-list.files(pattern="summary_HKY.*.RData")
  out_all_pre<-list()    #make a list to dump the batches into
  for(i in 1:length(to_read)){   # Dump the summ object from each batch into an element of a new list
    load(to_read[[i]])
    out_all_pre[[i]]<-summ
  }
  ## Combine everything into a list where each element of the list is a single analysis
  out_all<-unlist(out_all_pre, recursive=FALSE)
  

  ## Now that analyses are all laoded up, we can set the working directory to where we want to output figures, etc.
  setwd(write_to)
  
  # Pull out the analyses that failed
  failed_summs<-out_all[which(lapply(out_all, function(x) x[[1]])=="FAILED")]  # Get the analyses that failed
  
  # Then get the analyses that were successfully summarized
  diagnosticsAll<-out_all[which(!lapply(out_all, function(x) x[[1]])=="FAILED")]
  ## Name the diagnostics by analysis
  analysis_names<-lapply(diagnosticsAll, function(x) gsub("\\.nex|\\.run1.t.*", "", names(x$"ACT params")[[1]]))
  # Set these as the names for the elements of the diagnostics list
  names(diagnosticsAll)<-analysis_names
  
  ## If any analyses are duplicated, remove them
  dup_analysis<-analysis_names[which(duplicated(analysis_names))]
  if(length(dup_analysis)>0){
    diagnosticsAll<-diagnosticsAll[!names(diagnosticsAll)==dup_analysis]
  }
  
  
  # For each analysis, get out a list of which chains passed and failed the diagnostic tests
  pass_fail<-lapply(diagnosticsAll, function(x) as.data.frame(x[[1]]))
  pass_fail_combined<-do.call(rbind.fill, pass_fail) #Combine these into a single matrix where each row is a chain
  rownames(pass_fail_combined)<-unlist(lapply(diagnosticsAll, function(x) names(x$"ACT params")))
  
  ## Get the proportion of anlyses that pass each threshold for convergence
  prop_pass<-apply(pass_fail_combined, 2, function(x) length(which(x==TRUE))/length(x))
  prop_pass
  
  # whole analyses pass or fail asdsf
  analyses_pass_asdsf_HKY<-sapply(diagnosticsAll, function(x) unique((x[[1]][,"Pass ASDSF"])))
  
  
  
  ####### Identify how many diagnostics are failed for each chain, how many 
  ### parameters fail ESS for each chain
  # Identify parameters that failed ESS in each chain & number
  failed_ESS_pre<-unlist(lapply(diagnosticsAll, function(x) x$Params_ESS_less_200), recursive=FALSE)
  failed_ESS<-lapply(failed_ESS_pre, names)   # get the actual names of these parameters as the value in the object
  failed_ESS_num<-sapply(failed_ESS, length)
  # Identify how many diags are failed for each chain
  pass_fail_reduced<-pass_fail_combined[,!colnames(pass_fail_combined) %in% c("Pass topo pseudo ESS", "Pass ASDSF 2 chains", "Pass SF corr 2 chains", "Pass PSRF 2 chains")]  #cut out redundant diags (pseudo and approx topo ESS) and H&W
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
  hist(diff_num_failed_diags, main="HKY")
  
  ## Let's compare how many parameters fail the ESS between original and reanalyses - need to compare the proportion of parameters failing, because down lower for HKY and empirical base frequency reanalyses, the number of paramters is lower in the reanalyses
  fail_ESS_per_ana<-lapply(analysis_base_names, function(x) sum(failed_ESS_num[grep(x, names(failed_ESS_num))]))
  names(fail_ESS_per_ana)<-analysis_base_names
  fail_ESS_per_ana_orig<-lapply(analysis_base_names, function(x) sum(failed_ESS_num_orig[grep(x, names(failed_ESS_num_orig))]))
  names(fail_ESS_per_ana_orig)<-analysis_base_names
  num_params<-lapply(diagnosticsAll, function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))  ## Get the total number of parameters that ESS is calculated for (subtract 1 because topo us not used in calculation of number of params that fail ESS) - multiply by 2 because we have 2 chains
  ratio_ESS_failed<-unlist(fail_ESS_per_ana)/unlist(num_params) # ratio of failed ESS to total # pars
  ## Get the ratio for the original analyses
  num_params_orig<-lapply(raw_diags_orig[names(fail_ESS_per_ana_orig)], function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))
  ratio_ESS_failed_orig<-unlist(fail_ESS_per_ana_orig)/unlist(num_params_orig)
  # Difference in the ratios
  diff_num_failed_ESS<-unlist(ratio_ESS_failed_orig)-unlist(ratio_ESS_failed)
  mean(diff_num_failed_ESS)
  median(diff_num_failed_ESS)
  sd(diff_num_failed_ESS)
  hist(diff_num_failed_ESS, main="HKY")

  
  prop_pass_HKY<-prop_pass
  ratio_ESS_failed_HKY<-ratio_ESS_failed
  fail_diags_per_ana_HKY<-fail_diags_per_ana
  
  
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
  
  

  
########################################################################################################################
########################################################################################################################
###########  Same as above but for the set of empirical base frequency analyses
########################################################################################################################
########################################################################################################################
  ## Bring in summary
  ## set the working directory to where summaries of the chains are stored
  setwd(storage)
  
  # Read in all output - need to ensure that don't overwrite each time, as each file has same object names for summ
  to_read<-list.files(pattern="summary_EmpBF.*.RData")
  out_all_pre<-list()    #make a list to dump the batches into
  for(i in 1:length(to_read)){   # Dump the summ object from each batch into an element of a new list
    load(to_read[[i]])
    out_all_pre[[i]]<-summ
  }
  ## Combine everything into a list where each element of the list is a single analysis
  out_all<-unlist(out_all_pre, recursive=FALSE)
  
  
  ## Now that analyses are all laoded up, we can set the working directory to where we want to output figures, etc.
  setwd(write_to)
  
  # Pull out the analyses that failed
  failed_summs<-out_all[which(lapply(out_all, function(x) x[[1]])=="FAILED")]  # Get the analyses that failed
  
  # Then get the analyses that were successfully summarized
  diagnosticsAll<-out_all[which(!lapply(out_all, function(x) x[[1]])=="FAILED")]
  ## Name the diagnostics by analysis
  analysis_names<-lapply(diagnosticsAll, function(x) gsub("\\.nex|\\.run1.t.*", "", names(x$"ACT params")[[1]]))
  # Set these as the names for the elements of the diagnostics list
  names(diagnosticsAll)<-analysis_names
  
  ## If any analyses are duplicated, remove them
  dup_analysis<-analysis_names[which(duplicated(analysis_names))]
  if(length(dup_analysis)>0){
    diagnosticsAll<-diagnosticsAll[!names(diagnosticsAll)==dup_analysis]
  }
  
  
  # For each analysis, get out a list of which chains passed and failed the diagnostic tests
  pass_fail<-lapply(diagnosticsAll, function(x) as.data.frame(x[[1]]))
  pass_fail_combined<-do.call(rbind.fill, pass_fail) #Combine these into a single matrix where each row is a chain
  rownames(pass_fail_combined)<-unlist(lapply(diagnosticsAll, function(x) names(x$"ACT params")))
  
  ## Get the proportion of anlyses that pass each threshold for convergence
  prop_pass<-apply(pass_fail_combined, 2, function(x) length(which(x==TRUE))/length(x))
  prop_pass
  
  # whole analyses pass or fail asdsf
  analyses_pass_asdsf_EmpBF<-sapply(diagnosticsAll, function(x) unique((x[[1]][,"Pass ASDSF"])))
  
  
  
  ####### Identify how many diagnostics are failed for each chain, how many 
  ### parameters fail ESS for each chain
  # Identify parameters that failed ESS in each chain & number
  failed_ESS_pre<-unlist(lapply(diagnosticsAll, function(x) x$Params_ESS_less_200), recursive=FALSE)
  failed_ESS<-lapply(failed_ESS_pre, names)   # get the actual names of these parameters as the value in the object
  failed_ESS_num<-sapply(failed_ESS, length)
  # Identify how many diags are failed for each chain
  pass_fail_reduced<-pass_fail_combined[,!colnames(pass_fail_combined) %in% c("Pass topo pseudo ESS", "Pass ASDSF 2 chains", "Pass SF corr 2 chains", "Pass PSRF 2 chains")]  #cut out redundant diags (pseudo and approx topo ESS) and H&W
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
  hist(diff_num_failed_diags, main="EmpBF")
  
  ## Let's compare how many parameters fail the ESS between original and reanalyses - need to compare the proportion of parameters failing, because down lower for HKY and empirical base frequency reanalyses, the number of paramters is lower in the reanalyses
  fail_ESS_per_ana<-lapply(analysis_base_names, function(x) sum(failed_ESS_num[grep(x, names(failed_ESS_num))]))
  names(fail_ESS_per_ana)<-analysis_base_names
  fail_ESS_per_ana_orig<-lapply(analysis_base_names, function(x) sum(failed_ESS_num_orig[grep(x, names(failed_ESS_num_orig))]))
  names(fail_ESS_per_ana_orig)<-analysis_base_names
  num_params<-lapply(diagnosticsAll, function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))  ## Get the total number of parameters that ESS is calculated for (subtract 1 because topo us not used in calculation of number of params that fail ESS) - multiply by 2 because we have 2 chains
  ratio_ESS_failed<-unlist(fail_ESS_per_ana)/unlist(num_params) # ratio of failed ESS to total # pars
  ## Get the ratio for the original analyses
  num_params_orig<-lapply(raw_diags_orig[names(fail_ESS_per_ana_orig)], function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))
  ratio_ESS_failed_orig<-unlist(fail_ESS_per_ana_orig)/unlist(num_params_orig)
  # Difference in the ratios
  diff_num_failed_ESS<-unlist(ratio_ESS_failed_orig)-unlist(ratio_ESS_failed)
  mean(diff_num_failed_ESS)
  median(diff_num_failed_ESS)
  sd(diff_num_failed_ESS)
  hist(diff_num_failed_ESS, main="EmpBF")
  

  prop_pass_EmpBF<-prop_pass
  ratio_ESS_failed_EmpBF<-ratio_ESS_failed
  fail_diags_per_ana_EmpBF<-fail_diags_per_ana
  
  
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
  
  
########################################################################################################################
########################################################################################################################
###########  Same as above but for the set of Nst=Mixed analyses
########################################################################################################################
########################################################################################################################
  ## Bring in summary
  ## set the working directory to where summaries of the chains are stored
  setwd(storage)
  
  # Read in all output - need to ensure that don't overwrite each time, as each file has same object names for summ
  to_read<-list.files(pattern="summary_NstMixed.*.RData")
  out_all_pre<-list()    #make a list to dump the batches into
  for(i in 1:length(to_read)){   # Dump the summ object from each batch into an element of a new list
    load(to_read[[i]])
    out_all_pre[[i]]<-summ
  }
  ## Combine everything into a list where each element of the list is a single analysis
  out_all<-unlist(out_all_pre, recursive=FALSE)
  
  
  ## Now that analyses are all laoded up, we can set the working directory to where we want to output figures, etc.
  setwd(write_to)
  
  # Pull out the analyses that failed
  failed_summs<-out_all[which(lapply(out_all, function(x) x[[1]])=="FAILED")]  # Get the analyses that failed
  
  # Then get the analyses that were successfully summarized
  diagnosticsAll<-out_all[which(!lapply(out_all, function(x) x[[1]])=="FAILED")]
  ## Name the diagnostics by analysis
  analysis_names<-lapply(diagnosticsAll, function(x) gsub("\\.nex|\\.run1.t.*", "", names(x$"ACT params")[[1]]))
  # Set these as the names for the elements of the diagnostics list
  names(diagnosticsAll)<-analysis_names
  
  ## If any analyses are duplicated, remove them
  dup_analysis<-analysis_names[which(duplicated(analysis_names))]
  if(length(dup_analysis)>0){
    diagnosticsAll<-diagnosticsAll[!names(diagnosticsAll)==dup_analysis]
  }
  
  
  # For each analysis, get out a list of which chains passed and failed the diagnostic tests
  pass_fail<-lapply(diagnosticsAll, function(x) as.data.frame(x[[1]]))
  pass_fail_combined<-do.call(rbind.fill, pass_fail) #Combine these into a single matrix where each row is a chain
  rownames(pass_fail_combined)<-unlist(lapply(diagnosticsAll, function(x) names(x$"ACT params")))
  
  ## Get the proportion of anlyses that pass each threshold for convergence
  prop_pass<-apply(pass_fail_combined, 2, function(x) length(which(x==TRUE))/length(x))
  prop_pass
  
  # whole analyses pass or fail asdsf
  analyses_pass_asdsf_NstMixed<-sapply(diagnosticsAll, function(x) unique((x[[1]][,"Pass ASDSF"])))
  
  

  ####### Identify how many diagnostics are failed for each chain, how many 
  ### parameters fail ESS for each chain
  # Identify parameters that failed ESS in each chain & number
  failed_ESS_pre<-unlist(lapply(diagnosticsAll, function(x) x$Params_ESS_less_200), recursive=FALSE)
  failed_ESS<-lapply(failed_ESS_pre, names)   # get the actual names of these parameters as the value in the object
  failed_ESS_num<-sapply(failed_ESS, length)
  # Identify how many diags are failed for each chain
  pass_fail_reduced<-pass_fail_combined[,!colnames(pass_fail_combined) %in% c("Pass topo pseudo ESS", "Pass ASDSF 2 chains", "Pass SF corr 2 chains", "Pass PSRF 2 chains")]  #cut out redundant diags (pseudo and approx topo ESS) and H&W
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
  hist(diff_num_failed_diags, main="NstMixed")
  
  ## Let's compare how many parameters fail the ESS between original and reanalyses - need to compare the proportion of parameters failing, because down lower for HKY and empirical base frequency reanalyses, the number of paramters is lower in the reanalyses
  fail_ESS_per_ana<-lapply(analysis_base_names, function(x) sum(failed_ESS_num[grep(x, names(failed_ESS_num))]))
  names(fail_ESS_per_ana)<-analysis_base_names
  fail_ESS_per_ana_orig<-lapply(analysis_base_names, function(x) sum(failed_ESS_num_orig[grep(x, names(failed_ESS_num_orig))]))
  names(fail_ESS_per_ana_orig)<-analysis_base_names
  num_params<-lapply(diagnosticsAll, function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))  ## Get the total number of parameters that ESS is calculated for (subtract 1 because topo us not used in calculation of number of params that fail ESS) - multiply by 2 because we have 2 chains
  ratio_ESS_failed<-unlist(fail_ESS_per_ana)/unlist(num_params) # ratio of failed ESS to total # pars
  ## Get the ratio for the original analyses
  num_params_orig<-lapply(raw_diags_orig[names(fail_ESS_per_ana_orig)], function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))
  ratio_ESS_failed_orig<-unlist(fail_ESS_per_ana_orig)/unlist(num_params_orig)
  # Difference in the ratios
  diff_num_failed_ESS<-unlist(ratio_ESS_failed_orig)-unlist(ratio_ESS_failed)
  mean(diff_num_failed_ESS)
  median(diff_num_failed_ESS)
  sd(diff_num_failed_ESS)
  hist(diff_num_failed_ESS, main="NstMixed")
  
  
  prop_pass_NstMixed<-prop_pass
  ratio_ESS_failed_NstMixed<-ratio_ESS_failed
  fail_diags_per_ana_NstMixed<-fail_diags_per_ana
  
  
  
  ### ESS ratio table
  ESS_ratio_table<-cbind(ratio_ESS_failed_orig_all, rep(NA, 20), rep(NA, 20), rep(NA, 20), rep(NA, 20))
  ESS_ratio_table[match(names(ratio_ESS_failed_runLonger), names(ratio_ESS_failed_orig_all)),2]<-ratio_ESS_failed_runLonger
  ESS_ratio_table[match(names(ratio_ESS_failed_HKY), names(ratio_ESS_failed_orig_all)),3]<-ratio_ESS_failed_HKY
  ESS_ratio_table[match(names(ratio_ESS_failed_EmpBF), names(ratio_ESS_failed_orig_all)),4]<-ratio_ESS_failed_EmpBF
  ESS_ratio_table[match(names(ratio_ESS_failed_NstMixed), names(ratio_ESS_failed_orig_all)),5]<-ratio_ESS_failed_NstMixed
  colnames(ESS_ratio_table)<-c("original", "run longer", "HKY", "Emp BF", "Nst=Mixed")
  write.csv(ESS_ratio_table, file="Table_S2_ESS_ratio_table.csv")
  
  ### Diags pass table
  Diags_pass_table<-cbind(fail_diags_per_ana_orig_all, rep(NA, 20), rep(NA, 20), rep(NA, 20), rep(NA, 20))
  Diags_pass_table[match(names(fail_diags_per_ana_runLonger), names(fail_diags_per_ana_orig_all)),2]<-fail_diags_per_ana_runLonger
  Diags_pass_table[match(names(fail_diags_per_ana_HKY), names(fail_diags_per_ana_orig_all)),3]<-fail_diags_per_ana_HKY
  Diags_pass_table[match(names(fail_diags_per_ana_EmpBF), names(fail_diags_per_ana_orig_all)),4]<-fail_diags_per_ana_EmpBF
  Diags_pass_table[match(names(fail_diags_per_ana_NstMixed), names(fail_diags_per_ana_orig_all)),5]<-fail_diags_per_ana_NstMixed
  colnames(Diags_pass_table)<-c("original", "run longer", "HKY", "Emp BF", "Nst=Mixed")
  write.csv(Diags_pass_table, file="Table_S3_Diags_pass_table.csv")
  
  
  ### Pass ASDSF table
  ASDSF_pass_table<-cbind(analyses_pass_asdsf_orig, rep(NA, 20), rep(NA, 20), rep(NA, 20), rep(NA, 20))
  ASDSF_pass_table[match(names(analyses_pass_asdsf_runLonger), names(analyses_pass_asdsf_orig)),2]<-analyses_pass_asdsf_runLonger
  ASDSF_pass_table[match(names(analyses_pass_asdsf_HKY), names(analyses_pass_asdsf_orig)),3]<-analyses_pass_asdsf_HKY
  ASDSF_pass_table[match(names(analyses_pass_asdsf_EmpBF), names(analyses_pass_asdsf_orig)),4]<-analyses_pass_asdsf_EmpBF
  ASDSF_pass_table[match(names(analyses_pass_asdsf_NstMixed), names(analyses_pass_asdsf_orig)),5]<-analyses_pass_asdsf_NstMixed
  colnames(ASDSF_pass_table)<-c("original", "run longer", "HKY", "Emp BF", "Nst=Mixed")
  write.csv(ASDSF_pass_table, file="Table_S4_ASDSF_pass_table.csv")
  
  