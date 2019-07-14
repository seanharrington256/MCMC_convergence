## Script to compare Geweke's diagnostics with an initial window of 0.1 and 0.3

# Set working directory
setwd("Gewekes_wide_output")

# read in the summaries
to_read<-list.files(pattern="summary")

out_all_pre<-list()    #make a list to dump the batches into
for(i in 1:length(to_read)){   # Dump the summ object from each batch into an element of a new list
  load(to_read[[i]])
  out_all_pre[[i]]<-summ
}

## Combine everything into a list where each element of the list is a single analysis
out_all<-unlist(out_all_pre, recursive=FALSE)

# Pull out the analyses that failed
failed_summs<-out_all[which(lapply(out_all, function(x) x[[1]])=="FAILED")]  # Get the analyses that failed

# Then get the analyses that were successfully summarized
wide_gew<-out_all[which(!lapply(out_all, function(x) x[[1]])=="FAILED")]
# Which pass and fail
pass_fail<-lapply(wide_gew, function(x) x[[1]])


# For each analyses, get out a list of which chains passed and failed the geweke's tests
pass_fail<-unlist(pass_fail)
# Get the proportion of chains that pass Geweke's test
prop_gew_wide_pass<-length(which(pass_fail==TRUE))/length(pass_fail)
prop_gew_wide_pass  # This is slightly higher than the proprtion of chains that fail with the default window, but not by much - ~61% vs. ~57% with the default - this isn't the cause of the failures

# Identify the parameters that failed the Geweke's test in each chain - for chains that passed the test, this will return 0 length vectors, rather than NULL or NA
failed_geweke<-unlist(lapply(wide_gew, function(x) x$Params_fail_Geweke), recursive=FALSE)
names(failed_geweke)<-unlist(lapply(wide_gew, function(x) names(x$Params_fail_Geweke)), recursive=FALSE)

# Get only the chains that actually failed the test for one or more parameters
failed_geweke<-failed_geweke[which(lapply(failed_geweke, length)!=0)]
# Names of parameters that failed only
raw_names_failed_geweke<-unlist(lapply(failed_geweke, names))
# Which ones are these?
names_failed_geweke<-unique(raw_names_failed_geweke)
# Get number of times each parameter failed Geweke
num_failed_geweke<-lapply(names_failed_geweke, function(x) length(which(raw_names_failed_geweke==x)))

# Get the number of times each parameter shows up across analyses - parameters like kappa are only in a subset of analyses
pars_per_analysis<-unlist(lapply(wide_gew, function(x) names(x$Geweke_z_scores[[1]]$z))) # Start by getting the full list of the parameters that are included in all analyses
num_params_total<-lapply(names_failed_geweke, function(x) length(which(pars_per_analysis==x)))
# Proportion_fail_geweke
prop_fail_geweke<-unlist(num_failed_geweke)/unlist(num_params_total)
# Make a dataframe to plot the proportion with which each parameter fails
failed_geweke_to_plot<-cbind(unlist(names_failed_geweke), prop_fail_geweke)
## Plot this out
pdf(width=5, height=4.5, file="Fig_S2.4_proportion_fail_Geweke_WIDE.pdf")
barplot(as.numeric(failed_geweke_to_plot[,2]), names.arg=failed_geweke_to_plot[,1], las=2, col="blue")
dev.off()
pdf(width=5, height=4.5, file="Fig_S2.3_num_fail_perChain_Geweke_WIDE.pdf")
hist(unlist(lapply(failed_geweke, length)), main=NULL, xlab="# failed parameters/chain", col="blue")
dev.off()

## Look at how many chains fail multiple parameters for Geweke's
failed_geweke_multi<-failed_geweke[which(lapply(failed_geweke, length)>1)]  # Identify the analyses that failed Geweke's diagnostic for more than 1 parameter
length(failed_geweke_multi)## How many analyses fail multiple Geweke's?
length(failed_geweke_multi)/(length(wide_gew)*2) # Proportion of chains that fail multiple Geweke's
# this is still a lot of chains failing

