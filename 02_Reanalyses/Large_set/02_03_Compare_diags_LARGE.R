## This script compares the convergence diagnostics for large sets reestimated analyses to the original chains

# Load up the packages
library(ggplot2)
library(plyr)
library(phangorn)
library(rwty)




# This part of designating directories assumes that the parent directory of the directory containing this script contains directories named
#    diags_IG, diags_New_heat, diags_NstMixed, (these previous all output by running the 03_Calc_diags.R script on the )
#     sets of analyses to reanalyze) and Orig_diags (contains the "Orig_diags_for_IG_reestimated.RData" and "Orig_diags_for_reestimated_LARGE.RData"
#     files from the 04_Analyze_diagnostic_output.R script)


# only do this if the basedir object doesn't exist - prevents accidental rerunning of this step from another directory and messing up all paths
if(!exists("basedir")){
  basedir<-getwd()  # get the path to the base directory
}



## Set up the locations of the diagnostic summaries of each of the datasets
storage_IG<-paste0(basedir, "/diags_IG")
storage_NewHeat<-paste0(basedir, "/diags_New_heat")
storage_NstMixed<-paste0(basedir, "/diags_NstMixed")
# Set up the location of the original  burnin and diagnostics, etc.
orig_storage<-paste0(basedir, "/Orig_diags")
# Set up a directory to write output to
write_to<-paste0(basedir, "/output")
setwd(write_to)

### Read in the original diags for comparison below
setwd(orig_storage)
load("Orig_diags_for_IG_reestimated.RData")
load("Orig_diags_for_reestimated_LARGE.RData")

#########################################################################################################
#########################################################################################################
############# start by comparing G to IG
#########################################################################################################
#########################################################################################################
## Bring in summaries
## set the working directory to where summaries of the chains are stored
setwd(storage_IG)
# Read in all output - need to ensure that don't overwrite each time, as each file has same object names for summ
to_read<-list.files(pattern="summary_")
out_all_IG_1<-list()    #make a list to dump the batches into
for(i in 1:length(to_read)){   # Dump the summ object from each batch into an element of a new list
  load(to_read[[i]])
  out_all_IG_1[[i]]<-summ
}
## Combine everything into a single list, where each element of the list is a single chain
out_all_IG<-unlist(out_all_IG_1, recursive=FALSE)

## Now that analyses are all loaded up, we can set the working directory to where we want to output figures, etc.
setwd(write_to)

# Pull out the analyses that failed
failed_summs_IG<-out_all_IG[which(lapply(out_all_IG, function(x) x[[1]])=="FAILED")]  # Get the analyses that failed

# Then get the analyses that were successfully summarized
diagnosticsAll<-out_all_IG[which(!lapply(out_all_IG, function(x) x[[1]])=="FAILED")]
## Name the diagnostics by analysis
analysis_names<-lapply(diagnosticsAll, function(x) gsub("\\.nex|\\.run1.t.*", "", names(x$"ACT params")[[1]]))
# Set these as the names for the elements of the diagnostics list
names(diagnosticsAll)<-analysis_names
## Remove some diags that accidentally got included here
# diagnosticsAll<-diagnosticsAll[!gsub("_processed_mbReady", "", names(diagnosticsAll)) %in% c("my_ENSGALG00000012355.macse_DNA_gb", "ensGene.ENST00000464418.1.inc", "morf_452603_11codonalign.trim_interval.ctx")]



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

# Split things out based on which set of reanalyses they were in - remembering that analyses may be in more than one category
setwd(orig_storage)

## Bad ESS
bad_ESS<-read.csv("high_IG_cor_bad_ESS_LARGE.csv", stringsAsFactors = FALSE, row.names=1)[[1]]
# Pull out the analyses that match these
indices_bad_ESS<-unlist(lapply(bad_ESS, function(x) which(gsub("\\.nex.*$|\\.run.*$|_processed_mbReady", "", rownames(pass_fail_combined))==x)))
prop_pass_bad_ESS<-apply(pass_fail_combined[indices_bad_ESS,], 2, function(x) length(which(x==TRUE))/length(x))

## Good ESS
good_ESS<-read.csv("high_IG_cor_good_conv_LARGE.csv", stringsAsFactors = FALSE, row.names=1)[[1]]
# Pull out the analyses that match these
indices_good_ESS<-unlist(lapply(good_ESS,  function(x) which(gsub("\\.nex.*$|\\.run.*$|_processed_mbReady", "",rownames(pass_fail_combined))==x)))
prop_pass_good_ESS<-apply(pass_fail_combined[indices_good_ESS,], 2, function(x) length(which(x==TRUE))/length(x))



####### Identify how many parameters fail ESS for each chain
# Identify parameters that failed ESS in each chain & number
failed_ESS_pre<-unlist(lapply(diagnosticsAll, function(x) x$Params_ESS_less_200), recursive=FALSE)
failed_ESS<-lapply(failed_ESS_pre, names)   # get the actual names of these parameters as the value in the object
failed_ESS_num<-sapply(failed_ESS, length)
# Identify how many diags are failed for each chain
pass_fail_reduced<-pass_fail_combined[,!colnames(pass_fail_combined) %in% c("Pass topo pseudo ESS", "Pass ASDSF 2 chains", "Pass SF corr 2 chains", "Pass PSRF 2 chains")]  #cut out redundant diags
num_failed_diags<-apply(pass_fail_reduced, 1, function(x) length(which(x==FALSE)))


# Pass/fail diags
pass_fail_orig<-orig_IG_diags_LARGE[[1]] #pass/fail table for original diagnostics
pass_fail_reduced_orig<-pass_fail_orig[,!colnames(pass_fail_orig) %in% c("Pass topo pseudo ESS", "Pass ASDSF 2 chains", "Pass SF corr 2 chains", "Pass PSRF 2 chains")]  #cut out redundant diags
num_failed_diags_orig<-apply(pass_fail_reduced_orig, 1, function(x) length(which(x==FALSE)))
## Pars fail ESS
failed_ESS_orig<-orig_IG_diags_LARGE[[2]] #pars fail ESS for original diags
failed_ESS_num_orig<-sapply(failed_ESS_orig, length)# get numnber of pars that originally failed ESS
## Raw diags
raw_diags_orig<-orig_IG_diags_LARGE[[3]] #all raw original diags

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

## Let's compare how many parameters fail the ESS between original and reanalyses - compare the proportion of parameters failing, because this changes with some reanalyses
fail_ESS_per_ana<-lapply(analysis_base_names, function(x) sum(failed_ESS_num[grep(x, names(failed_ESS_num))]))
names(fail_ESS_per_ana)<-analysis_base_names
fail_ESS_per_ana_orig<-lapply(analysis_base_names, function(x) sum(failed_ESS_num_orig[grep(x, names(failed_ESS_num_orig))]))
names(fail_ESS_per_ana_orig)<-analysis_base_names
num_params<-lapply(diagnosticsAll, function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))  ## Get the total number of parameters that ESS is calculated for by getting the length of ESS values for first chain in an analysis (subtract 1 because topo us not used in calculation of number of params that fail ESS) - multiply by length of $ESS_values - number of chains
ratio_ESS_failed<-unlist(fail_ESS_per_ana)/unlist(num_params) # ratio of failed ESS to total # pars

## Get the ratio for the original analyses
num_params_orig<-lapply(raw_diags_orig[names(fail_ESS_per_ana_orig)], function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))
ratio_ESS_failed_orig<-unlist(fail_ESS_per_ana_orig)/unlist(num_params_orig)
# Difference in the ratios
diff_num_failed_ESS<-unlist(ratio_ESS_failed_orig)-unlist(ratio_ESS_failed)
mean(diff_num_failed_ESS)
median(diff_num_failed_ESS)
sd(diff_num_failed_ESS)
hist(diff_num_failed_ESS)

## combine these all into a table
diffs_reanalysis_IG<-cbind(analysis_base_names, diff_num_failed_diags, diff_num_failed_ESS, rep("", length(diff_num_failed_ESS)))  ## Add empty placeholder column that will be replaced with which subset each analysis belongs to
## add in a column that specifies which of the subset analyses each of these was from
diffs_reanalysis_IG[which(diffs_reanalysis_IG[,"analysis_base_names"] %in% bad_ESS), 4]<-sapply(diffs_reanalysis_IG[which(diffs_reanalysis_IG[,"analysis_base_names"] %in% bad_ESS), 4], function(x) paste(x, "Bad_ESS", sep=""))
diffs_reanalysis_IG[which(diffs_reanalysis_IG[,"analysis_base_names"] %in% good_ESS), 4]<-sapply(diffs_reanalysis_IG[which(diffs_reanalysis_IG[,"analysis_base_names"] %in% good_ESS), 4], function(x) paste(x, "Good_ESS", sep=""))

## Write out a csv with this table in it
setwd(write_to)
write.csv(diffs_reanalysis_IG, file="Diffs_IG_large.csv")

## Get the mean, median, and standard deviation for each of the subsets
## Get the mean, median, and standard deviation for each of the subsets
# Means
summ_diffs_bad_ESS<-apply(diffs_reanalysis_IG[which(diffs_reanalysis_IG[,"analysis_base_names"] %in% bad_ESS), 2:3], 2, function(x) mean(as.numeric(x)))
summ_diffs_good_ESS<-apply(diffs_reanalysis_IG[which(diffs_reanalysis_IG[,"analysis_base_names"] %in% good_ESS), 2:3], 2, function(x) mean(as.numeric(x)))
summ_diffs<-rbind(summ_diffs_bad_ESS, summ_diffs_good_ESS)
colnames(summ_diffs)<-c("Failed diags", "failed ESS")
rownames(summ_diffs)<-c("Bad ESS", "Good ESS")
write.csv(summ_diffs, file="Mean_diffs_summ_IG.csv")
# Medians
summ_diffs_bad_ESS_med<-apply(diffs_reanalysis_IG[which(diffs_reanalysis_IG[,"analysis_base_names"] %in% bad_ESS), 2:3], 2, function(x) median(as.numeric(x)))
summ_diffs_good_ESS_med<-apply(diffs_reanalysis_IG[which(diffs_reanalysis_IG[,"analysis_base_names"] %in% good_ESS), 2:3], 2, function(x) median(as.numeric(x)))
summ_diffs_med<-rbind(summ_diffs_bad_ESS_med, summ_diffs_good_ESS_med)
colnames(summ_diffs_med)<-c("Failed diags", "failed ESS")
rownames(summ_diffs_med)<-c("Bad ESS", "Good ESS")
write.csv(summ_diffs_med, file="Median_diffs_summ_IG.csv")
# Std Deviation
summ_diffs_bad_ESS_sd<-apply(diffs_reanalysis_IG[which(diffs_reanalysis_IG[,"analysis_base_names"] %in% bad_ESS), 2:3], 2, function(x) sd(as.numeric(x)))
summ_diffs_good_ESS_sd<-apply(diffs_reanalysis_IG[which(diffs_reanalysis_IG[,"analysis_base_names"] %in% good_ESS), 2:3], 2, function(x) sd(as.numeric(x)))
summ_diffs_sd<-rbind(summ_diffs_bad_ESS_sd, summ_diffs_good_ESS_sd)
colnames(summ_diffs_sd)<-c("Failed diags", "failed ESS")
rownames(summ_diffs_sd)<-c("Bad ESS", "Good ESS")
write.csv(summ_diffs_sd, file="StdDev_diffs_summ_IG.csv")


########################################################################################################################
########################################################################################################################
###########  Compare analyses with new heating to originals
########################################################################################################################
########################################################################################################################
## Bring in summaries
## set the working directory to where summaries of the chains are stored
setwd(storage_NewHeat)
# Read in all output - need to ensure that don't overwrite each time, as each file has same object names for summ
to_read<-list.files(pattern="summary_")
out_all_heat_1<-list()    #make a list to dump the batches into
for(i in 1:length(to_read)){   # Dump the summ object from each batch into an element of a new list
  load(to_read[[i]])
  out_all_heat_1[[i]]<-summ
}
## Combine everything into a single list, where each element of the list is a single chain
out_all_heat<-unlist(out_all_heat_1, recursive=FALSE)

## Now that analyses are all loaded up, we can set the working directory to where we want to output figures, etc.
setwd(write_to)

# Pull out the analyses that failed
failed_summs_heat<-out_all_heat[which(lapply(out_all_heat, function(x) x[[1]])=="FAILED")]  # Get the analyses that failed

# Then get the analyses that were successfully summarized
diagnostics<-out_all_heat[which(!lapply(out_all_heat, function(x) x[[1]])=="FAILED")]
## Name the diagnostics by analysis
analysis_names<-lapply(diagnostics, function(x) gsub("\\.nex|\\.run1.t.*|_processed_mbReady", "", names(x$"ACT params")[[1]]))
# Set these as the names for the elements of the diagnostics list
names(diagnostics)<-analysis_names
## Remove some diags that accidentally got included here - this was from an earlier, buggy run
# diagnostics<-diagnostics[!gsub("_processed_mbReady", "", names(diagnostics)) %in% c("16190", "chr5_3418", "chr6_1350", "ensGene.ENST00000458043.1", "knownGene.uc001dig.3.1")]



## If any analyses are duplicated, remove them
dup_analysis<-analysis_names[which(duplicated(analysis_names))]
if(length(dup_analysis)>0){
  diagnostics<-diagnostics[!names(diagnostics)==dup_analysis]
}


# Decided to not deal with analyses that failed Geweke's, so remove those here
# but because of potential overlap among analyses in the different subsets, can't just remove the Geweke's ones, instead keep everything else
setwd(orig_storage)
bad_ESS<-read.csv("fail_oneormore_ESS_LARGE.csv", stringsAsFactors = FALSE, row.names=1)[[1]]
bad_topo_bad_LnL<-read.csv("fail_topo_bad_LnL_ESS_LARGE.csv", stringsAsFactors = FALSE, row.names=1)[[1]]
bad_topo_good_LnL<-read.csv("fail_topo_good_LnL_ESS_LARGE.csv", stringsAsFactors = FALSE, row.names=1)[[1]]
bad_asdsf<-read.csv("fail_ASDSF_LARGE.csv", stringsAsFactors = FALSE, row.names=1)[[1]]

diagnostics<-diagnostics[which(names(diagnostics) %in% c(bad_ESS, bad_topo_bad_LnL, bad_topo_good_LnL, bad_asdsf))]



# For each analysis, get out a list of which chains passed and failed the diagnostic tests
pass_fail<-lapply(diagnostics, function(x) as.data.frame(x[[1]]))
pass_fail_combined<-do.call(rbind.fill, pass_fail) #Combine these into a single matrix where each row is a chain
rownames(pass_fail_combined)<-unlist(lapply(diagnostics, function(x) names(x$"ACT params")))

## Get the proportion of anlyses that pass each threshold for convergence
prop_pass<-apply(pass_fail_combined, 2, function(x) length(which(x==TRUE))/length(x))
prop_pass


# Split things out based on which set of reanalyses they were in - remembering that one analysis is in more than one category
setwd(orig_storage)

## Fail asdsf
# Pull out the analyses that match these
indices_bad_asdsf<-unlist(lapply(bad_asdsf, function(x) which(gsub("\\.nex.*$|\\.run.*$|_processed_mbReady", "",rownames(pass_fail_combined))==x)))
prop_pass_bad_asdsf<-apply(pass_fail_combined[indices_bad_asdsf,], 2, function(x) length(which(x==TRUE))/length(x))

## fail one or more ESS
# Pull out the analyses that match these
indices_bad_ESS<-unlist(lapply(bad_ESS,  function(x) which(gsub("\\.nex.*$|\\.run.*$|_processed_mbReady", "",rownames(pass_fail_combined))==x)))
prop_pass_bad_ESS<-apply(pass_fail_combined[indices_bad_ESS,], 2, function(x) length(which(x==TRUE))/length(x))

## fail topological ESS and LnL ESS
# Pull out the analyses that match these
indices_bad_topo_bad_LnL<-unlist(lapply(bad_topo_bad_LnL,  function(x) which(gsub("\\.nex.*$|\\.run.*$|_processed_mbReady", "",rownames(pass_fail_combined))==x)))
prop_pass_bad_topo_bad_LnL<-apply(pass_fail_combined[indices_bad_topo_bad_LnL,], 2, function(x) length(which(x==TRUE))/length(x))

## fail topological ESS but pass LnL ESS
# Pull out the analyses that match these
indices_bad_topo_good_LnL<-unlist(lapply(bad_topo_good_LnL,  function(x) which(gsub("\\.nex.*$|\\.run.*$|_processed_mbReady", "",rownames(pass_fail_combined))==x)))
prop_pass_bad_topo_good_LnL<-apply(pass_fail_combined[indices_bad_topo_good_LnL,], 2, function(x) length(which(x==TRUE))/length(x))

####### Identify how many parameters fail ESS for each chain
# Identify parameters that failed ESS in each chain & number
failed_ESS_pre<-unlist(lapply(diagnostics, function(x) x$Params_ESS_less_200), recursive=FALSE)
failed_ESS<-lapply(failed_ESS_pre, names)   # get the actual names of these parameters as the value in the object
failed_ESS_num<-sapply(failed_ESS, length)
# Identify how many diags are failed for each chain
pass_fail_reduced<-pass_fail_combined[,!colnames(pass_fail_combined) %in% c("Pass topo pseudo ESS", "Pass ASDSF 2 chains", "Pass SF corr 2 chains", "Pass PSRF 2 chains")]  #cut out redundant diags
num_failed_diags<-apply(pass_fail_reduced, 1, function(x) length(which(x==FALSE)))


# Pass/fail diags
pass_fail_orig<-orig_diags_LARGE[[1]] #pass/fail table for original diagnostics
pass_fail_reduced_orig<-pass_fail_orig[,!colnames(pass_fail_orig) %in% c("Pass topo pseudo ESS", "Pass ASDSF 2 chains", "Pass SF corr 2 chains", "Pass PSRF 2 chains")]  #cut out redundant diags
num_failed_diags_orig<-apply(pass_fail_reduced_orig, 1, function(x) length(which(x==FALSE)))
## Pars fail ESS
failed_ESS_orig<-orig_diags_LARGE[[2]] #pars fail ESS for original diags
failed_ESS_num_orig<-sapply(failed_ESS_orig, length)# get numnber of pars that originally failed ESS
## Raw diags
raw_diags_orig<-orig_diags_LARGE[[3]] #all raw original diags

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

      ## see if we can figure out which are better this time around, specifically
      
      props_pass_newHeat<-rbind((apply(pass_fail_reduced_orig, 2, sum)/7004), (apply(pass_fail_reduced, 2, sum)/7004))
      rownames(props_pass_newHeat)<-c("original", "new heating")


## Let's compare how many parameters fail the ESS between original and reanalyses - compare the proportion of parameters failing, because this changes with some reanalyses
fail_ESS_per_ana<-lapply(analysis_base_names, function(x) sum(failed_ESS_num[grep(x, names(failed_ESS_num))]))
names(fail_ESS_per_ana)<-analysis_base_names
fail_ESS_per_ana_orig<-lapply(analysis_base_names, function(x) sum(failed_ESS_num_orig[grep(x, names(failed_ESS_num_orig))]))
names(fail_ESS_per_ana_orig)<-analysis_base_names
num_params<-lapply(diagnostics, function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))  ## Get the total number of parameters that ESS is calculated for by getting the length of ESS values for first chain in an analysis (subtract 1 because topo us not used in calculation of number of params that fail ESS) - multiply by length of $ESS_values - number of chains
ratio_ESS_failed<-unlist(fail_ESS_per_ana)/unlist(num_params) # ratio of failed ESS to total # pars
## Get the ratio for the original analyses
num_params_orig<-lapply(raw_diags_orig[names(fail_ESS_per_ana_orig)], function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))
ratio_ESS_failed_orig<-unlist(fail_ESS_per_ana_orig)/unlist(num_params_orig)
# Difference in the ratios
diff_num_failed_ESS<-unlist(ratio_ESS_failed_orig)-unlist(ratio_ESS_failed)
mean(diff_num_failed_ESS)
median(diff_num_failed_ESS)
sd(diff_num_failed_ESS)
hist(diff_num_failed_ESS)

## combine these all into a table
diffs_reanalysis_heat<-cbind(analysis_base_names, diff_num_failed_diags, diff_num_failed_ESS, rep("", length(diff_num_failed_ESS)))  ## Add empty placeholder column that will be replaced with which subset each analysis belongs to
## add in a column that specifies which of the subset analyses each of these was from
diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_asdsf), 4]<-sapply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_asdsf), 4], function(x) paste(x, "bad_asdsf", sep=""))
diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_ESS), 4]<-sapply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_ESS), 4], function(x) paste(x, "bad_ESS", sep=""))
diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_topo_bad_LnL), 4]<-sapply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_topo_bad_LnL), 4], function(x) paste(x, "bad_topo_bad_LnL", sep=""))
diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_topo_good_LnL), 4]<-sapply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_topo_good_LnL), 4], function(x) paste(x, "bad_topo_good_LnL", sep=""))


## Write out a csv with this table in it
setwd(write_to)
write.csv(diffs_reanalysis_heat, file="Diffs_heat_large.csv")

## Get the mean, median, and standard deviation for each of the subsets
# Means
summ_diffs_bad_asdsf<-apply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_asdsf), 2:3], 2, function(x) mean(as.numeric(x)))
summ_diffs_bad_ESS<-apply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_ESS), 2:3], 2, function(x) mean(as.numeric(x)))
summ_diffs_bad_topo_bad_LnL<-apply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_topo_bad_LnL), 2:3], 2, function(x) mean(as.numeric(x)))
summ_diffs_bad_topo_good_LnL<-apply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_topo_good_LnL), 2:3], 2, function(x) mean(as.numeric(x)))
summ_diffs<-rbind(summ_diffs_bad_asdsf, summ_diffs_bad_ESS, summ_diffs_bad_topo_bad_LnL, summ_diffs_bad_topo_good_LnL)
colnames(summ_diffs)<-c("Failed diags", "failed ESS")
rownames(summ_diffs)<-c("fail asdsf", "fail ESS", "fail topo ESS fail LnL ESS", "fail topo ESS pass Lnl ESS")
write.csv(summ_diffs, file="Mean_diffs_summ_heat.csv")
# Medians
summ_diffs_bad_asdsf_med<-apply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_asdsf), 2:3], 2, function(x) median(as.numeric(x)))
summ_diffs_bad_ESS_med<-apply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_ESS), 2:3], 2, function(x) median(as.numeric(x)))
summ_diffs_bad_topo_bad_LnL_med<-apply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_topo_bad_LnL), 2:3], 2, function(x) median(as.numeric(x)))
summ_diffs_bad_topo_good_LnL_med<-apply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_topo_good_LnL), 2:3], 2, function(x) median(as.numeric(x)))
summ_diffs_med<-rbind(summ_diffs_bad_asdsf_med, summ_diffs_bad_ESS_med, summ_diffs_bad_topo_bad_LnL_med, summ_diffs_bad_topo_good_LnL_med)
colnames(summ_diffs_med)<-c("Failed diags", "failed ESS")
rownames(summ_diffs_med)<-c("fail asdsf", "fail ESS", "fail topo ESS fail LnL ESS", "fail topo ESS pass Lnl ESS")
write.csv(summ_diffs_med, file="Median_diffs_summ_heat.csv")
# Std Deviation
summ_diffs_bad_asdsf_sd<-apply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_asdsf), 2:3], 2, function(x) sd(as.numeric(x)))
summ_diffs_bad_ESS_sd<-apply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_ESS), 2:3], 2, function(x) sd(as.numeric(x)))
summ_diffs_bad_topo_bad_LnL_sd<-apply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_topo_bad_LnL), 2:3], 2, function(x) sd(as.numeric(x)))
summ_diffs_bad_topo_good_LnL_sd<-apply(diffs_reanalysis_heat[which(diffs_reanalysis_heat[,"analysis_base_names"] %in% bad_topo_good_LnL), 2:3], 2, function(x) sd(as.numeric(x)))
summ_diffs_sd<-rbind(summ_diffs_bad_asdsf_sd, summ_diffs_bad_ESS_sd, summ_diffs_bad_topo_bad_LnL_sd, summ_diffs_bad_topo_good_LnL_sd)
colnames(summ_diffs_sd)<-c("Failed diags", "failed ESS")
rownames(summ_diffs_sd)<-c("fail asdsf", "fail ESS", "fail topo ESS fail LnL ESS", "fail topo ESS pass Lnl ESS")
write.csv(summ_diffs_sd, file="StdDev_diffs_summ_heat.csv")


########################################################################################################################
########################################################################################################################
###########  Compare analyses with new nst=mixed to originals
########################################################################################################################
########################################################################################################################
## Bring in summaries
## set the working directory to where summaries of the chains are stored
setwd(storage_NstMixed)
# Read in all output - need to ensure that don't overwrite each time, as each file has same object names for summ
to_read<-list.files(pattern="summary_")
out_all_nstMixed_1<-list()    #make a list to dump the batches into
for(i in 1:length(to_read)){   # Dump the summ object from each batch into an element of a new list
  load(to_read[[i]])
  out_all_nstMixed_1[[i]]<-summ
}
## Combine everything into a single list, where each element of the list is a single chain
out_all_nstMixed<-unlist(out_all_nstMixed_1, recursive=FALSE)

## Now that analyses are all loaded up, we can set the working directory to where we want to output figures, etc.
setwd(write_to)

# Pull out the analyses that failed
failed_summs_nstMixed<-out_all_nstMixed[which(lapply(out_all_nstMixed, function(x) x[[1]])=="FAILED")]  # Get the analyses that failed

# Then get the analyses that were successfully summarized
diagnostics<-out_all_nstMixed[which(!lapply(out_all_nstMixed, function(x) x[[1]])=="FAILED")]
## Name the diagnostics by analysis
analysis_names<-lapply(diagnostics, function(x) gsub("\\.nex|\\.run1.t.*|_processed_mbReady", "", names(x$"ACT params")[[1]]))
# Set these as the names for the elements of the diagnostics list
names(diagnostics)<-analysis_names
## Remove some diags that accidentally got included here
# diagnostics<-diagnostics[!gsub("_processed_mbReady", "", names(diagnostics)) %in% c("16190", "chr5_3418", "chr6_1350", "ensGene.ENST00000458043.1", "knownGene.uc001dig.3.1")]


## If any analyses are duplicated, remove them
dup_analysis<-analysis_names[which(duplicated(analysis_names))]
if(length(dup_analysis)>0){
  diagnostics<-diagnostics[!names(diagnostics)==dup_analysis]
}

setwd(orig_storage)
bad_ESS<-read.csv("fail_oneormore_ESS_LARGE.csv", stringsAsFactors = FALSE, row.names=1)[[1]]
bad_topo_bad_LnL<-read.csv("fail_topo_bad_LnL_ESS_LARGE.csv", stringsAsFactors = FALSE, row.names=1)[[1]]
bad_topo_good_LnL<-read.csv("fail_topo_good_LnL_ESS_LARGE.csv", stringsAsFactors = FALSE, row.names=1)[[1]]
bad_asdsf<-read.csv("fail_ASDSF_LARGE.csv", stringsAsFactors = FALSE, row.names=1)[[1]]

diagnostics<-diagnostics[which(names(diagnostics) %in% c(bad_asdsf, bad_ESS, bad_topo_bad_LnL, bad_topo_good_LnL))]

# For each analysis, get out a list of which chains passed and failed the diagnostic tests
pass_fail<-lapply(diagnostics, function(x) as.data.frame(x[[1]]))
pass_fail_combined<-do.call(rbind.fill, pass_fail) #Combine these into a single matrix where each row is a chain
rownames(pass_fail_combined)<-unlist(lapply(diagnostics, function(x) names(x$"ACT params")))


## Get the proportion of anlyses that pass each threshold for convergence
prop_pass<-apply(pass_fail_combined, 2, function(x) length(which(x==TRUE))/length(x))
prop_pass

# Split things out based on which set of reanalyses they were in - remembering that one analysis is in more than one category
setwd(orig_storage)

## Fail asdsf
# Pull out the analyses that match these
indices_bad_asdsf<-unlist(lapply(bad_asdsf, function(x) which(gsub("\\.nex.*$|\\.run.*$|_processed_mbReady", "",rownames(pass_fail_combined))==x)))
prop_pass_bad_asdsf<-apply(pass_fail_combined[indices_bad_asdsf,], 2, function(x) length(which(x==TRUE))/length(x))

## fail one or more ESS
# Pull out the analyses that match these
indices_bad_ESS<-unlist(lapply(bad_ESS,  function(x) which(gsub("\\.nex.*$|\\.run.*$|_processed_mbReady", "",rownames(pass_fail_combined))==x)))
prop_pass_bad_ESS<-apply(pass_fail_combined[indices_bad_ESS,], 2, function(x) length(which(x==TRUE))/length(x))

## fail topological ESS and LnL ESS
# Pull out the analyses that match these
indices_bad_topo_bad_LnL<-unlist(lapply(bad_topo_bad_LnL,  function(x) which(gsub("\\.nex.*$|\\.run.*$|_processed_mbReady", "",rownames(pass_fail_combined))==x)))
prop_pass_bad_topo_bad_LnL<-apply(pass_fail_combined[indices_bad_topo_bad_LnL,], 2, function(x) length(which(x==TRUE))/length(x))

## fail topological ESS but pass LnL ESS
# Pull out the analyses that match these
indices_bad_topo_good_LnL<-unlist(lapply(bad_topo_good_LnL,  function(x) which(gsub("\\.nex.*$|\\.run.*$|_processed_mbReady", "",rownames(pass_fail_combined))==x)))
prop_pass_bad_topo_good_LnL<-apply(pass_fail_combined[indices_bad_topo_good_LnL,], 2, function(x) length(which(x==TRUE))/length(x))

####### Identify how many parameters fail ESS for each chain
# Identify parameters that failed ESS in each chain & number
failed_ESS_pre<-unlist(lapply(diagnostics, function(x) x$Params_ESS_less_200), recursive=FALSE)
failed_ESS<-lapply(failed_ESS_pre, names)   # get the actual names of these parameters as the value in the object
failed_ESS_num<-sapply(failed_ESS, length)
# Identify how many diags are failed for each chain
pass_fail_reduced<-pass_fail_combined[,!colnames(pass_fail_combined) %in% c("Pass topo pseudo ESS", "Pass ASDSF 2 chains", "Pass SF corr 2 chains", "Pass PSRF 2 chains")]  #cut out redundant diags
num_failed_diags<-apply(pass_fail_reduced, 1, function(x) length(which(x==FALSE)))


# Pass/fail diags
pass_fail_orig<-orig_diags_LARGE[[1]] #pass/fail table for original diagnostics
pass_fail_reduced_orig<-pass_fail_orig[,!colnames(pass_fail_orig) %in% c("Pass topo pseudo ESS", "Pass ASDSF 2 chains", "Pass SF corr 2 chains", "Pass PSRF 2 chains")]  #cut out redundant diags
num_failed_diags_orig<-apply(pass_fail_reduced_orig, 1, function(x) length(which(x==FALSE)))
## Pars fail ESS
failed_ESS_orig<-orig_diags_LARGE[[2]] #pars fail ESS for original diags
failed_ESS_num_orig<-sapply(failed_ESS_orig, length)# get numnber of pars that originally failed ESS
## Raw diags
raw_diags_orig<-orig_diags_LARGE[[3]] #all raw original diags

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


    ## see if we can figure out which are better this time around, specifically
    props_pass_nst_mixed<-rbind((apply(pass_fail_reduced_orig, 2, sum)/7004), (apply(pass_fail_reduced, 2, sum)/7004))
    rownames(props_pass_nst_mixed)<-c("original", "Nst=mixed")
    
    props_pass_changes<-rbind(props_pass_nst_mixed, props_pass_newHeat[2,])
    rownames(props_pass_changes)[3]<-"New heat"
    setwd(write_to)
    write.csv(props_pass_changes, "proportion_diags_passNstNewHeat.csv")

## Let's compare how many parameters fail the ESS between original and reanalyses - compare the proportion of parameters failing, because this changes with some reanalyses
fail_ESS_per_ana<-lapply(analysis_base_names, function(x) sum(failed_ESS_num[grep(x, names(failed_ESS_num))]))
names(fail_ESS_per_ana)<-analysis_base_names
fail_ESS_per_ana_orig<-lapply(analysis_base_names, function(x) sum(failed_ESS_num_orig[grep(x, names(failed_ESS_num_orig))]))
names(fail_ESS_per_ana_orig)<-analysis_base_names
num_params<-lapply(diagnostics, function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))  ## Get the total number of parameters that ESS is calculated for by getting the length of ESS values for first chain in an analysis (subtract 1 because topo us not used in calculation of number of params that fail ESS) - multiply by length of $ESS_values - number of chains
ratio_ESS_failed<-unlist(fail_ESS_per_ana)/unlist(num_params) # ratio of failed ESS to total # pars
## Get the ratio for the original analyses
num_params_orig<-lapply(raw_diags_orig[names(fail_ESS_per_ana_orig)], function(x) length(x$ESS_values)*(length(x$ESS_values[[1]])-1))
ratio_ESS_failed_orig<-unlist(fail_ESS_per_ana_orig)/unlist(num_params_orig)
# Difference in the ratios
diff_num_failed_ESS<-unlist(ratio_ESS_failed_orig)-unlist(ratio_ESS_failed)
mean(diff_num_failed_ESS)
median(diff_num_failed_ESS)
sd(diff_num_failed_ESS)
hist(diff_num_failed_ESS)



## combine these all into a table
diffs_reanalysis_nstMixed<-cbind(analysis_base_names, diff_num_failed_diags, diff_num_failed_ESS, rep("", length(diff_num_failed_ESS)))  ## Add empty placeholder column that will be replaced with which subset each analysis belongs to
## add in a column that specifies which of the subset analyses each of these was from
diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_asdsf), 4]<-sapply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_asdsf), 4], function(x) paste(x, "bad_asdsf", sep=""))
diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_ESS), 4]<-sapply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_ESS), 4], function(x) paste(x, "bad_ESS", sep=""))
diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_topo_bad_LnL), 4]<-sapply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_topo_bad_LnL), 4], function(x) paste(x, "bad_topo_bad_LnL", sep=""))
diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_topo_good_LnL), 4]<-sapply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_topo_good_LnL), 4], function(x) paste(x, "bad_topo_good_LnL", sep=""))


## Write out a csv with this table in it
setwd(write_to)
write.csv(diffs_reanalysis_nstMixed, file="Diffs_nstMixed_large.csv")

## Get the mean, median, and standard deviation for each of the subsets
# Means
summ_diffs_bad_asdsf<-apply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_asdsf), 2:3], 2, function(x) mean(as.numeric(x)))
summ_diffs_bad_ESS<-apply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_ESS), 2:3], 2, function(x) mean(as.numeric(x)))
summ_diffs_bad_topo_bad_LnL<-apply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_topo_bad_LnL), 2:3], 2, function(x) mean(as.numeric(x)))
summ_diffs_bad_topo_good_LnL<-apply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_topo_good_LnL), 2:3], 2, function(x) mean(as.numeric(x)))
summ_diffs<-rbind(summ_diffs_bad_asdsf, summ_diffs_bad_ESS, summ_diffs_bad_topo_bad_LnL, summ_diffs_bad_topo_good_LnL)
colnames(summ_diffs)<-c("Failed diags", "failed ESS")
rownames(summ_diffs)<-c("fail asdsf", "fail ESS", "fail topo ESS fail LnL ESS", "fail topo ESS pass Lnl ESS")
write.csv(summ_diffs, file="Mean_diffs_summ_nstMixed.csv")
# Medians
summ_diffs_bad_asdsf_med<-apply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_asdsf), 2:3], 2, function(x) median(as.numeric(x)))
summ_diffs_bad_ESS_med<-apply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_ESS), 2:3], 2, function(x) median(as.numeric(x)))
summ_diffs_bad_topo_bad_LnL_med<-apply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_topo_bad_LnL), 2:3], 2, function(x) median(as.numeric(x)))
summ_diffs_bad_topo_good_LnL_med<-apply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_topo_good_LnL), 2:3], 2, function(x) median(as.numeric(x)))
summ_diffs_med<-rbind(summ_diffs_bad_asdsf_med, summ_diffs_bad_ESS_med, summ_diffs_bad_topo_bad_LnL_med, summ_diffs_bad_topo_good_LnL_med)
colnames(summ_diffs_med)<-c("Failed diags", "failed ESS")
rownames(summ_diffs_med)<-c("fail asdsf", "fail ESS", "fail topo ESS fail LnL ESS", "fail topo ESS pass Lnl ESS")
write.csv(summ_diffs_med, file="Median_diffs_summ_nstMixed.csv")
# Std Deviation
summ_diffs_bad_asdsf_sd<-apply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_asdsf), 2:3], 2, function(x) sd(as.numeric(x)))
summ_diffs_bad_ESS_sd<-apply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_ESS), 2:3], 2, function(x) sd(as.numeric(x)))
summ_diffs_bad_topo_bad_LnL_sd<-apply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_topo_bad_LnL), 2:3], 2, function(x) sd(as.numeric(x)))
summ_diffs_bad_topo_good_LnL_sd<-apply(diffs_reanalysis_nstMixed[which(diffs_reanalysis_nstMixed[,"analysis_base_names"] %in% bad_topo_good_LnL), 2:3], 2, function(x) sd(as.numeric(x)))
summ_diffs_sd<-rbind(summ_diffs_bad_asdsf_sd, summ_diffs_bad_ESS_sd, summ_diffs_bad_topo_bad_LnL_sd, summ_diffs_bad_topo_good_LnL_sd)
colnames(summ_diffs_sd)<-c("Failed diags", "failed ESS")
rownames(summ_diffs_sd)<-c("fail asdsf", "fail ESS", "fail topo ESS fail LnL ESS", "fail topo ESS pass Lnl ESS")
write.csv(summ_diffs_sd, file="StdDev_diffs_summ_nstMixed.csv")

