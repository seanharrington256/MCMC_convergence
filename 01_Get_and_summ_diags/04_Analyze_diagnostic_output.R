### Script to summarize diagnotic output generated from the 03_Calc_diags.R script

## Set up some directories
storage_all<-"Diag_output"  # location of the summaries of diagnostic output output by 03_Calc_diags.R
write_to<-"Figs"  # directory to dump various figures into 
figs_ms<-"Figs/figs_ms/"   # diretory to dump figures specifically used for the manuscript into
manual_dir<-"../05A_Manual_Reanalyses"  # location of the folder containing the Manually_reanalyze.csv file

setwd(write_to)


# Load up the packages
library(ggplot2)
library(plyr)
library(plotly)
library(scatterplot3d)
library(rgl)
library(VennDiagram)
library(dplyr)
library(pcaMethods)
library(ggridges)





## Bring in summaries of chains
  ## set the working directory to where summaries of the chains are stored
  setwd(storage_all)
  
  # Read in all output - need to ensure that don't overwrite each time, as each file has same object names for summ
  to_read<-list.files(pattern="summary.*.RData")
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
diagnostics<-out_all[which(!lapply(out_all, function(x) x[[1]])=="FAILED")]
## Name the diagnostics by analysis
analysis_names<-lapply(diagnostics, function(x) gsub("\\.nex|\\.run1.t.*", "", names(x$Burnin[1])))
# Set these as the names for the elements of the diagnostics list
names(diagnostics)<-analysis_names

## If any analyses are duplicated, remove them
dup_analysis<-analysis_names[which(duplicated(analysis_names))]
if(length(dup_analysis)>0){
  diagnostics<-diagnostics[!names(diagnostics)==dup_analysis]
}
# Also remove morf_448818_11codonalign.trim_interval.ctx, becuase it had no Nexus file 
diagnostics<-diagnostics[!names(diagnostics) %in% "morf_448818_11codonalign.trim_interval.ctx"]


# For each analyses, get out a list of which chains passed and failed the diagnostic tests
pass_fail<-lapply(diagnostics, function(x) x[[1]])
pass_fail_combined<-do.call(rbind, pass_fail) #Combine these into a single matrix where each row is a chain
rownames(pass_fail_combined)<-unlist(lapply(diagnostics, function(x) names(x$Burnin)))

## Get the proportion of anlyses that pass each threshold for convergence
prop_pass<-apply(pass_fail_combined, 2, function(x) length(which(x==TRUE))/length(x))
prop_pass

# Make a csv of this
  prop_pass_table1<-rbind(gsub("Pass ", "", names(prop_pass)), prop_pass)
  prop_pass_table2<-gsub("ESS", "ESS > 200", prop_pass_table1)
  prop_pass_table3<-gsub("ASDSF", "ASDSF < 0.01", prop_pass_table2)
  prop_pass_table4<-gsub("SF corr", "SF corr > 0.9", prop_pass_table3)
  prop_pass_table5<-gsub("PSRF", "PSRF < 1.02", prop_pass_table4)
  prop_pass_table6<-gsub("Geweke", "Geweke p > 0.05", prop_pass_table5)
  prop_pass_table<-prop_pass_table6[,!colnames(prop_pass_table6) %in% "Pass H&W"]
  colnames(prop_pass_table)<-rep(NA, ncol(prop_pass_table))
  colnames(prop_pass_table)[[1]]<-"Proportion of chains passing each diagnostic"
  setwd(figs_ms)
  write.csv(prop_pass_table, file="table_pass_diags.csv")
  setwd(write_to)

## Get the average number of generations to ESS value of 200 ##
    gens_200_all<-unlist(lapply(diagnostics, function(x) x$Generations_to_ESS_200), recursive=FALSE)    # Get a list where every element contains the gens to 200 for a single chain
    param_names<-unique(unlist(lapply(gens_200_all, FUN=function(x) names(x))))    ## Different analyses included different parameters (e.g., GTR vs HKY), need to find the full list of parameters that were used
    gens_200_table<-do.call(rbind, lapply(gens_200_all, function(x) x[match(param_names, names(x))]))    # Use this to order all parameters the same and make a table of gens to ess 200 for all for each chain
    colnames(gens_200_table)<-param_names
    
    # Get mean number of gens and standard deviation for each param
    mean_gens_200<-apply(gens_200_table, 2, mean, na.rm=TRUE)
    median_gens_200<-apply(gens_200_table, 2, median, na.rm=TRUE)
    sd_gens_200<-apply(gens_200_table, 2, sd, na.rm=TRUE)
    # Rank these parameters and see what converges fastest
    sorted_means_gens_200<-sort(mean_gens_200)
    sorted_median_gens_200<-sort(median_gens_200)
    sorted_sd_gens_200<-sort(sd_gens_200)
    sorted_means_gens_200
    sorted_median_gens_200
    sorted_sd_gens_200
  
    ## Plot out gens to ESS 200 with ridgeline plot
    # First transform the dataframe into a format easily readable by ggplot
    ess_200_melted<-reshape2::melt(gens_200_table, id.vars = NULL)
    # Not every parameter is in every analysis, so need to prune out the NA's
    ess_200_melted<-ess_200_melted[which(is.na(ess_200_melted$value)==FALSE),]
    ## See how many times each paramter shows up - 
    numpars_ESS<-lapply(levels(ess_200_melted$Var2), function(x) length(which(ess_200_melted$Var2==x)))
    names(numpars_ESS)<-levels(ess_200_melted$Var2)
    ## add number of times par shows up into par names for plotting
    plot_ESS_200<-ess_200_melted
    for(i in 1:length(levels(plot_ESS_200$Var2))){
      levels(plot_ESS_200$Var2)[[i]]<-paste(levels(plot_ESS_200$Var2)[[i]], " (", numpars_ESS[levels(plot_ESS_200$Var2)[[i]]], ")", sep="") 
    }
    # Not logged
    setwd(figs_ms)
    pdf(file="Fig_1_Gens_ESS_not_logged.pdf", width=6, height=4)
    ggplot(plot_ESS_200, aes(y=Var2, x=value, fill = Var2)) +
      labs(title=NULL,y="Parameter (# analyses/parameter)", x = "Generations to ESS of 200")+
      scale_fill_cyclical(values = c("aquamarine2", "darkviolet"))+
      scale_x_continuous(limits = c(2800000, 5500000))+ # cut off the long right tail
      stat_density_ridges(quantile_lines = TRUE, quantiles = 2, scale=2.5)+  # add in a line for the median of each
      scale_y_discrete(expand = expand_scale(add = c(0.2, 2.8)))
    dev.off()
    
    setwd(write_to)
    

    
## Get the average ESS for each parameter ##
    # Unlist so that every individual chain is an element of the list
    ess_all<-unlist(lapply(diagnostics, function(x) x$ESS_values), recursive=FALSE)
    ## Different analyses included different parameters (e.g., GTR vs HKY), need to find the full list of parameters that were used
    param_names_ess<-unique(unlist(lapply(ess_all, FUN=function(x) names(x))))
    # Use this to order all parameters the same and make a table of ess for each chain
    ess_table<-do.call(rbind, lapply(ess_all, function(x) x[match(param_names_ess, names(x))]))
    rownames(ess_table)<-names(ess_all)
    colnames(ess_table)<-param_names_ess
    # Get mean ess and standard deviation for each param
    mean_ess<-apply(ess_table, 2, mean, na.rm=TRUE)
    median_ess<-apply(ess_table, 2, median, na.rm=TRUE)
    sd_ess<-apply(ess_table, 2, sd, na.rm=TRUE)
    # Rank these parameters and see what converges fastest
    sorted_means_ess<-sort(mean_ess)
    sorted_median_ess<-sort(median_ess)
    sorted_sd_ess<-sort(sd_ess)
    sorted_means_ess
    sorted_median_ess
    sorted_sd_ess
    
    ##plot out ridgeline density plots
    # First transform the dataframe into a format easily readable by ggplot
    ess_melted<-reshape2::melt(ess_table, id.vars = NULL)
    # Not every parameter is in every analysis, so need to prune out the NA's
    ess_melted<-ess_melted[which(is.na(ess_melted$value)==FALSE),]
    setwd(figs_ms)
    # ESS values
    pdf(file="Fig_2_Final_ESS.pdf", width=6, height=4)
    ggplot(ess_melted, aes(y=Var2, x=value, fill = Var2)) +
      labs(title=NULL, y ="Parameter", x = "Final ESS")+
      stat_density_ridges(quantile_lines = TRUE, quantiles = 2, scale=3.5)+
      scale_fill_cyclical(values = c("aquamarine2", "darkviolet"))+
      scale_y_discrete(expand = expand_scale(add = c(0.2, 2.8)))+
      scale_x_continuous(limits = c(500, 1010))
    dev.off()
    setwd(write_to)
    
# Get out the autocorrelation times for each parameter
    act_params_all_pre<-unlist(lapply(diagnostics, function(x) x$"ACT params"), recursive=FALSE)
    act_params_all<-lapply(act_params_all_pre, unlist)  # formatting to unlist within each element, now we have a list that contains named numeric vectors of act for parameters in each single chain
 # Get out the topological autocorrelation time that is used in calculation of approximate topo ESS
    topo_act<-unlist(lapply(diagnostics, function(x) lapply(x$"topo act", function(y) y[[1]])), recursive=FALSE)
    # For each chain, concatenate the parameter act with the topology act
    for(i in 1:length(topo_act)){
      names(topo_act[[i]])<-"topo"    # add on topo names
    } 
    params_topo_act<-list()
    for(i in 1:length(act_params_all)){
      params_topo_act[[i]]<-c(act_params_all[[i]], topo_act[[i]])
    }
    
    # Handle the act similarly to ESS above - first have to deal with the fact that different analyses have different parameters
    param_names_act<-unique(unlist(lapply(params_topo_act, FUN=function(x) names(x))))
    # Use this to order all parameters the same and make a table of act for each chain
    act_table<-do.call(rbind, lapply(params_topo_act, function(x) x[match(param_names_act, names(x))]))
    colnames(act_table)<-param_names_act
    rownames(act_table)<-gsub("^.*?\\.|nex\\.run|\\.t", "",names(topo_act))
    # Get mean ess and standard deviation for each param
    mean_act<-apply(act_table, 2, mean, na.rm=TRUE)
    median_act<-apply(act_table, 2, median, na.rm=TRUE)
    sd_act<-apply(act_table, 2, sd, na.rm=TRUE)
    # Rank these parameters and see what converges fastest
    sorted_means_act<-sort(mean_act)
    sorted_median_act<-sort(median_act)
    sorted_sd_act<-sort(sd_act)
    sorted_means_act
    sorted_median_act
    sorted_sd_act
    
    ## Make ridgeline plots of act values
    # First transform the dataframe into a format easily readable by ggplot
    act_melted<-reshape2::melt(act_table, id.vars = NULL)
    # Not every parameter is in every analysis, so need to prune out the NA's
    act_melted<-act_melted[which(is.na(act_melted$value)==FALSE),]
    ##Raw values
    ggplot(act_melted, aes(y=Var2, x=value, fill = Var2)) +
      labs(title="ACT", y ="Parameter", x = "ACT")+
      geom_density_ridges(scale=3)+ 
      scale_fill_cyclical(values = c("aquamarine2", "darkviolet"))
    ## Set the x limit to not include the long right tail
    setwd(figs_ms)
    pdf(file="Fig_3_ACT_plot.pdf", width=6, height=4)
    ggplot(act_melted, aes(y=Var2, x=value, fill = Var2)) +
      labs(title=NULL, y ="Parameter", x = "ACT")+
      scale_fill_cyclical(values = c("aquamarine2", "darkviolet"))+
      scale_x_continuous(limits = c(0.95, 2.05), breaks=c(1,1.5,2))+
      stat_density_ridges(quantile_lines = TRUE, quantiles = 2, scale=3)+
      scale_y_discrete(expand = expand_scale(add = c(0.2, 2.8)))+
      xlab("Autocorrelation time")
    dev.off()
    setwd(write_to)
    
    
    
# Compare topo pseudo ESS to approximate ESS
    approx_ess_pre<-unlist(lapply(diagnostics, function(x) x$"Approx topo ESS"), recursive=FALSE) # Get out just the approx ESS part of the diagnostics list
    approx_ess_operator<-lapply(approx_ess_pre, function(x) x$operator)  # Get out the operator for each element
    approx_ess_value<-lapply(approx_ess_pre, function(x) x$approx.ess)  # Get out the value for each element
    to_compare<-which(approx_ess_operator=="=") # make an object for indexing to compare only those analyses for which approx ESS was calculated with "="
    approx_for_compare<-approx_ess_value[to_compare]  # Pull out just the ESS values that have "="
    pseudo_ESS_conf<-unlist(lapply(diagnostics, function(x) x$"Topo pseudo ess confidence"), recursive=FALSE) # Get out just the approx ESS part of the diagnostics list
    pseudo_ESS_med_for_compare<-lapply(pseudo_ESS_conf, function(x) x$median.ess)[to_compare]
    diffs_approx_pseudo<-unlist(approx_for_compare)-unlist(pseudo_ESS_med_for_compare)
    mean(diffs_approx_pseudo)
  #Compare the approximate ESS to the maximum pseudo ESS  
    pseudo_ESS_upper_for_compare<-lapply(pseudo_ESS_conf, function(x) x$ci.upper)[to_compare]
    diffs_approx_pseudo_upper<-unlist(approx_for_compare)-unlist(pseudo_ESS_upper_for_compare)
    mean(diffs_approx_pseudo_upper)

    ## Plot out as overlapping histograms
    setwd(figs_ms)
    pdf(width=5, height=4.5, file="Fig_S1_Diffs_approx_pseudo_ESS.pdf")
    hist(diffs_approx_pseudo_upper, main=NULL, xlab="Approximate minus psuedo topological ESS", ylab="Frequency", breaks=50, col=rgb(1,0,0,0.5))
    hist(diffs_approx_pseudo, add=TRUE, breaks=50, col=rgb(0,0,1,0.5))
    legend("topright", c("Upper 95% pseudo", "median Pseudo"), col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch=15, cex=0.7)
    dev.off()
    setwd(write_to)    
    
    
## Identify which parameters failed Geweke's test
############
############
  # Identify the parameters that failed the Geweke's test in each chain - for chains that passed the test, this will return 0 length vectors, rather than NULL or NA
  failed_geweke<-unlist(lapply(diagnostics, function(x) x$Params_fail_Geweke), recursive=FALSE)
  names(failed_geweke)<-unlist(lapply(diagnostics, function(x) names(x$Params_fail_Geweke)), recursive=FALSE)
  failed_geweke<-failed_geweke[which(lapply(failed_geweke, length)!=0)]  # Get only the chains that actually failed the test for one or more parameters
  raw_names_failed_geweke<-unlist(lapply(failed_geweke, names))  # Names of parameters that failed only
  names_failed_geweke<-unique(raw_names_failed_geweke)  # Which ones are these?
  num_failed_geweke<-lapply(names_failed_geweke, function(x) length(which(raw_names_failed_geweke==x))) # Get number of times each parameter failed Geweke
  # Get the number of times each parameter shows up across analyses - parameters like kappa are only in a subset of analyses
  pars_per_analysis<-unlist(lapply(diagnostics, function(x) names(x$ESS_values[[1]]))) # Start by getting the full list of the parameters that are included in all analyses
  num_params_total<-lapply(names_failed_geweke, function(x) length(which(pars_per_analysis==x)))
  prop_fail_geweke<-unlist(num_failed_geweke)/unlist(num_params_total)# Proportion_fail_geweke
  failed_geweke_to_plot<-cbind(unlist(names_failed_geweke), prop_fail_geweke) # Make a dataframe to plot the proportion with which each parameter fails
  ## Plot this out
  setwd(figs_ms)
  pdf(width=5, height=4.5, file="Fig_S2.2_proportion_fail_Geweke.pdf")
  barplot(as.numeric(failed_geweke_to_plot[,2]), names.arg=failed_geweke_to_plot[,1], las=2, col="blue")
  dev.off()
  setwd(write_to)    
  # How many parameters did each chain fail Geweke's for?
  setwd(figs_ms)
  pdf(width=5, height=4.5, file="Fig_S2.1_num_fail_perChain_Geweke.pdf")
  hist(unlist(lapply(failed_geweke, length)), main=NULL, xlab="# failed parameters/chain", col="blue")
  dev.off()
  setwd(write_to)    
  
  
  ## Look at Geweke's failures in more detail
  failed_geweke_multi<-failed_geweke[which(lapply(failed_geweke, length)>1)]  # Identify the analyses that failed Geweke's diagnostic for more than 1 parameter
  length(failed_geweke_multi)## How many analyses fail multiple Geweke's?
  length(failed_geweke_multi)/(length(diagnostics)*2) # Proportion of chains that fail multiple Geweke's
  # this is still a lot of chains failing
  
  ## Look at the normality of parameters in each chain - maybe non-normal parameter distributions cause more common failures
  normality<-unlist(lapply(diagnostics, function(x) x$Shap_test), recursive=FALSE) # Get just the info from the Shapiro-Wilk tests of normality for each chain
  names(normality)<-unlist(lapply(diagnostics, function(x) names(x$Shap_test)), recursive=FALSE)
  norm_pvals<-lapply(normality, function(y)  lapply(y, function(x)  x$p.value)) # get just the p value for each parameter using a rather unpleasant but functional set of nested lapplys
  # Use a loop to determine if each parameter that has failed Geweke's test is normal or not
  is_norm_fail_Gew<-list()
  for(i in 1:length(failed_geweke)){
    is_norm_fail_Gew[[i]]<-norm_pvals[names(failed_geweke)[[i]]][[1]][names(failed_geweke[[i]])]>0.05
  }
  names(is_norm_fail_Gew)<-names(failed_geweke)
  is_norm_FG_all<-unlist(is_norm_fail_Gew) # unlist this and just look at parameters without caring what analysis or chain they come from
  length(which(is_norm_FG_all==FALSE)) # 30,390 of the parameters that failed Geweke's are non-normal
  length(which(is_norm_FG_all==TRUE))  # but still 16,595 of the parameters that failed Geweke's are normal (or at least not detected as non-normal by Shapiro-Wilk test)
  length(which(is_norm_FG_all==FALSE))/(length(which(is_norm_FG_all==FALSE))+length(which(is_norm_FG_all==TRUE))) # about 65% of parameters that fail Geweke's are normal
  ## What is the overall proportion of chains that fail normality?
  length(which(unlist(norm_pvals)<=0.5))/(length(which(unlist(norm_pvals)<=0.5))+length(which(unlist(norm_pvals)>0.5))) # 86% of parameters are non-normal...
  
  ## Look at normality of chains that fail Geweke's but pass other diagnostics
  fg_pass_fail<-pass_fail_combined[names(failed_geweke),] # get the table of passes and failures for all chains that fail Geweke's
  fg_pass_fail<-fg_pass_fail[,!colnames(fg_pass_fail) %in% c("Pass Geweke", "Pass H&W")]  # Exclude H&W and Geweke
  num_faildiags_fail_g<-apply(fg_pass_fail, 1, function(x) length(which(x==FALSE))) # find out how many diagnostics are failed for each chain that fails at least 1 Geweke's test
  failG_pass_others<-rownames(fg_pass_fail[names(which(num_faildiags_fail_g==0)),])  # get out only the analyses that fail Geweke's for at least 1 parameter but pass all other diagnostic tests (excludign H&W test)
  ## Let's look and see if parameters are normal that fail Geweke's for a parameter or more but otherwise pass other diagnostics
  is_norm_fail_Gew_passOthers<-unlist(is_norm_fail_Gew[failG_pass_others])
  length(which(is_norm_fail_Gew_passOthers==FALSE)) # 25,456 are non-normal
  length(which(is_norm_fail_Gew_passOthers==TRUE))  # still 14,569 are normal (or at least not detected as non-normal by Shapiro-Wilk test)
  
  ## Identify the analyses that failed LnL Geweke's
  lnl_fail_geweke<-failed_geweke[which(lapply(failed_geweke, function(x) "LnL" %in% names(x))==TRUE)]
  length(lnl_fail_geweke) # number of chains that fail LnL for Geweke
  length(lnl_fail_geweke)/(length(diagnostics)*2)# proportion of chains that fail LnL for Geweke
  
  # Get out the pass/fail table for only those analyses that fail more than one parameter for Geweke
  fg_pass_fail<-pass_fail_combined[names(failed_geweke_multi),]
  fg_prop_pass<-apply(fg_pass_fail, 2, function(x) length(which(x==TRUE))/length(x))
  fg_prop_pass
  # no real indication in any of the other diagnostics that anything problematic is going on
  
  # See if ones that fail LnL for Geweke fail any other diagnostics?
  fglnl_pass_fail<-pass_fail_combined[names(lnl_fail_geweke),]
  fglnl_prop_pass<-apply(fglnl_pass_fail, 2, function(x) length(which(x==TRUE))/length(x))
  fglnl_prop_pass
  # No real indication of any problems here, either
  
  ## Table for supp that has prop pass for failing multiple Geweke's or for LnL
  Gew_supp_tab_1<-rbind(fg_prop_pass, fglnl_prop_pass)
  rownames(Gew_supp_tab_1)<-c("Fail multiple parameters", "Fail LnL")
  setwd(figs_ms)
  write.csv(Gew_supp_tab_1, file="Table_S2.1_Gew_multi_LnL_prop_pass.csv")
  setwd(write_to)
  
  # Identify analyses that fail Geweke's for >0.05 of parameters included (including LnL and LnPr, although these aren't actually parameters)
  names(ess_all)<-unlist(lapply(diagnostics, function(x) names(x$Params_fail_Geweke))) # assign names to Ess values split up by each chain that will match Geweke's names
  num_params_chains<-lapply(ess_all, length)[names(failed_geweke)] # Get out the number of parameters in each chain - indexing gets this for just the analyses that failed at least one parameter for Geweke
  num_params_fail_G<-lapply(failed_geweke, length)
  ## Loop to get out the proportion of parameters in each analysis that fail Geweke's test
  prop_params_fail_G<-list()
  for(i in 1:length(num_params_fail_G)){
    prop_params_fail_G[[i]]<-num_params_fail_G[[i]]/num_params_chains[[i]]
  }
  names(prop_params_fail_G)<-names(num_params_chains)
  which(prop_params_fail_G>0.05) # identify which of these fail Geweke's for more than 5% of parameters
  length(which(prop_params_fail_G>0.05)) # this is every chain that failed Geweke's test for any parameters 

  ## Let's see if z scores are correlated with shapiro wilk test
  # First let's get the z scores all into a big table
    gew_zscores_pre<-unlist(lapply(diagnostics, function(x) x$Geweke_z_scores), recursive=FALSE)
    names(gew_zscores_pre)<-unlist(lapply(diagnostics, function(x) names(x$Geweke_z_scores)), recursive=FALSE)
    gew_zscores<-lapply(gew_zscores_pre, function(x) x$z)
    ## find the full list of parameters that were used
    param_names_gew<-unique(unlist(lapply(gew_zscores, FUN=function(x) names(x))))
    # Use this to order all parameters the same and make a table
    gewZ_table<-do.call(rbind, lapply(gew_zscores, function(x) x[match(param_names_gew, names(x))]))
    colnames(gewZ_table)<-param_names_gew
  # Do something similar for normality W scores
    norm_W<-lapply(normality, function(y) sapply(y, function(x) x$statistic)) # ugly sapply in an lapply to get a vector of W scores for each chain - each chain as an element of a list
    norm_W<-lapply(norm_W, function(x) setNames(x, nm=gsub("\\.W", "", names(x)))) # delete a .W off the end of all parameter names
    # order all parameters the same as for Geweke's just above and make a table
    norm_W_table<-do.call(rbind, lapply(norm_W, function(x) x[match(param_names_gew, names(x))]))
    colnames(norm_W_table)<-param_names_gew
  # Get the correlations
    all(rownames(gewZ_table)==rownames(norm_W_table))  ##Check that the rows are in the same order
    all(colnames(gewZ_table)==colnames(norm_W_table))
    Z_W_cors<-list() # list to dump correlations into
    for(i in 1:ncol(gewZ_table)){   ## loop to go column by column of both tables and run correlations
      Z_W_cors[[i]]<-cor(gewZ_table[,i], norm_W_table[,i], use="complete.obs")
    }
    names(Z_W_cors)<-colnames(gewZ_table)
    Z_W_cors   # The correlation between z scores and W scores is very, very low. Effectively non-existent
  ## Let's see if z scores are correlated with ESS values at all
    all(rownames(gewZ_table)==rownames(ess_table))  # check that rows are the same - they're not
    rownames(gewZ_table)<-gsub("\\.run|\\.t$", "", rownames(gewZ_table))  ## Edit the rownames for gewZ_table
    all(rownames(gewZ_table)==rownames(ess_table)) # Now these match up
    z_ESS_cors<-list()
    for(i in 1:ncol(gewZ_table)){
      z_ESS_cors[[i]]<-cor(gewZ_table[,i], ess_table[,colnames(gewZ_table)[i]], use="complete.obs")
    }
    names(z_ESS_cors)<-colnames(gewZ_table)
    z_ESS_cors   ## Also pretty much completely unrelated
    
    ## Table of z score correlations with W scores and ESS
    tavble_z_cors<-rbind(Z_W_cors, z_ESS_cors)
    setwd(figs_ms)
    write.csv(tavble_z_cors, file="Table_S2.2_Z_score_cors.csv")
    setwd(write_to)
    
##########################################################################################
##########################################################################################
##########################################################################################
### Compare convergence to dataset properties
##########################################################################################
##########################################################################################
##########################################################################################

## Read in the csv that states which study each analysis is from
study_names<-read.csv("Study_names.csv")
study_names[,1]<-as.character(study_names[,1])
levels(study_names[,2])<-c("Barcoding", "Amniotes", "PhyLoTA", "mtDNA") # replace the citations with description of the study
    
## Get all of the ESS values into a single dataframe for easier manipulation
    ess_for_cors_pre<-lapply(diagnostics, function(x) x$ESS_values) # get the ess values out of the main diagnostics object
    ess_setup_1<-lapply(ess_for_cors_pre, function(x) do.call(rbind, x)) # within each element of the list, rbind the 2 chains together
    # Use a for loop to assign the name of each analysis as a column for each element of the list
    ess_setup_2<-list()
    for(i in 1:length(ess_setup_1)){
      ess_setup_2[[i]]<-cbind.data.frame(names(ess_setup_1)[[i]], ess_setup_1[[i]]) # cbind the name of the analysis into the matrix as the first column
      colnames(ess_setup_2[[i]])[[1]]<-"analysis_name" # name the first column
    }
    ess_setup_3<-do.call(rbind.fill, ess_setup_2) # rbind this all together filling in empty columns with NA using rbind.fill
    # Then, combine this with dataset characteristics
    data_chars<-do.call(rbind, lapply(diagnostics, function(x) x$analysis_pars[1:2]))
    ess_data_chars_1<-merge(ess_setup_3, data_chars, by.x="analysis_name", by.y="row.names", all.x=TRUE) # finally, merge together the ESS and dataset data
    ess_data_chars<-merge(ess_data_chars_1, study_names, by.x="analysis_name", by.y="V1", all.x=TRUE) # finally finally, add in the study citation
    colnames(ess_data_chars)[which(colnames(ess_data_chars)=="V2")]<-"study"
    
    cor(ess_data_chars$nchar, ess_data_chars$ntax, use="complete.obs")  ## are number of taxa and characters correlated? -- nope
    
    
## look at the numbers of taxa and characters across analyses
    col.index<-as.factor(ess_data_chars$study)
    setwd(figs_ms)
    pdf(file="Dataset_props.pdf", height=4, width=4)
    plot(ess_data_chars$ntax, ess_data_chars$nchar, cex=0.3, col=c("darkviolet", "darkgoldenrod3", "deepskyblue", "coral2")[col.index], pch=20, xlab="# taxa", ylab="# characters")
    legend(x="topright", legend = levels(col.index), col=c("darkviolet", "darkgoldenrod3", "deepskyblue", "coral2"), pch=20, pt.cex=1, cex=0.7)
    dev.off()
    setwd(write_to)

# Correlations and plots of ess for different parameters against Ntax and Nchars
    # Make a list of the different parameters that have ESS values
    par_ess<-colnames(ess_data_chars)[!colnames(ess_data_chars) %in% c("analysis_name", "ntax", "nchar", "study")]
    
    # Get the correlations of each of these parameters with ntax
    ess_cors_Ntax<-sapply(par_ess, function(x) cor(ess_data_chars[,"ntax"], ess_data_chars[,x], use="complete.obs"))
    # Get the correlations of each of these parameters with nchars
    ess_cors_Nchars<-sapply(par_ess, function(x) cor(ess_data_chars[,"nchar"], ess_data_chars[,x], use="complete.obs"))
        
    ## Check out if these correlations vary by dataset
      # Amniotes - LnPr gets removed out during the correlation step because it wasn't sampled for any of these analyses
        amn_ess_data_chars<-ess_data_chars[which(ess_data_chars[,"study"]=="Amniotes"),]
        amn_params_remove<-names(which(sapply(amn_ess_data_chars, function(x)all(is.na(x)))))  # Find any columns that are completely NA to also remove these from the list of params to correlate
        par_amn_ess<-colnames(amn_ess_data_chars)[!colnames(amn_ess_data_chars) %in% c("analysis_name", "ntax", "nchar", "study", amn_params_remove)]
        # Get the correlations of each of these parameters with ntax
        amn_ess_cors_Ntax<-sapply(par_amn_ess, function(x) cor(amn_ess_data_chars[,"ntax"], amn_ess_data_chars[,x], use="complete.obs"))
        # Get the correlations of each of these parameters with nchars
        amn_ess_cors_Nchars<-sapply(par_amn_ess, function(x) cor(amn_ess_data_chars[,"nchar"], amn_ess_data_chars[,x], use="complete.obs"))
        ## Correlations here are very low, probably because of the low variation in the number of taxa
      # Barcoding - pinvar gets removed out during the correlation step because it wasn't sampled for any of these analyses
        bc_ess_data_chars<-ess_data_chars[which(ess_data_chars[,"study"]=="Barcoding"),]
        bc_params_remove<-names(which(sapply(bc_ess_data_chars, function(x)all(is.na(x)))))  # Find any columns that are completely NA to also remove these from the list of params to correlate
        par_bc_ess<-colnames(bc_ess_data_chars)[!colnames(bc_ess_data_chars) %in% c("analysis_name", "ntax", "nchar", "study", bc_params_remove)]
        # Get the correlations of each of these parameters with ntax
        bc_ess_cors_Ntax<-sapply(par_bc_ess, function(x) cor(bc_ess_data_chars[,"ntax"], bc_ess_data_chars[,x], use="complete.obs"))
        # Get the correlations of each of these parameters with nchars
        bc_ess_cors_Nchars<-sapply(par_bc_ess, function(x) cor(bc_ess_data_chars[,"nchar"], bc_ess_data_chars[,x], use="complete.obs"))
      # mtDNA
        mt_ess_data_chars<-ess_data_chars[which(ess_data_chars[,"study"]=="mtDNA"),]
        mt_params_remove<-names(which(sapply(mt_ess_data_chars, function(x)all(is.na(x)))))  # Find any columns that are completely NA to also remove these from the list of params to correlate
        par_mt_ess<-colnames(mt_ess_data_chars)[!colnames(mt_ess_data_chars) %in% c("analysis_name", "ntax", "nchar", "study", mt_params_remove)]
        # Get the correlations of each of these parameters with ntax
        mt_ess_cors_Ntax<-sapply(par_mt_ess, function(x) cor(mt_ess_data_chars[,"ntax"], mt_ess_data_chars[,x], use="complete.obs"))
        # Get the correlations of each of these parameters with nchars
        mt_ess_cors_Nchars<-sapply(par_mt_ess, function(x) cor(mt_ess_data_chars[,"nchar"], mt_ess_data_chars[,x], use="complete.obs"))
      # Phylota
        ph_ess_data_chars<-ess_data_chars[which(ess_data_chars[,"study"]=="PhyLoTA"),]
        ph_params_remove<-names(which(sapply(ph_ess_data_chars, function(x)all(is.na(x)))))  # Find any columns that are completely NA to also remove these from the list of params to correlate
        par_ph_ess<-colnames(ph_ess_data_chars)[!colnames(ph_ess_data_chars) %in% c("analysis_name", "ntax", "nchar", "study", ph_params_remove)]
        # Get the correlations of each of these parameters with ntax
        ph_ess_cors_Ntax<-sapply(par_ph_ess, function(x) cor(ph_ess_data_chars[,"ntax"], ph_ess_data_chars[,x], use="complete.obs"))
        # Get the correlations of each of these parameters with nchars
        ph_ess_cors_Nchars<-sapply(par_ph_ess, function(x) cor(ph_ess_data_chars[,"nchar"], ph_ess_data_chars[,x], use="complete.obs"))
        
    ## Make a big table with the correlations of each parameter's ESS with number of taxa and characters - include the correlation when considering all datasets together and when looking at each individually
        ESS_correlation_table<-bind_rows(ess_cors_Ntax, amn_ess_cors_Ntax, bc_ess_cors_Ntax, mt_ess_cors_Ntax, ph_ess_cors_Ntax, ess_cors_Nchars, amn_ess_cors_Nchars, bc_ess_cors_Nchars, mt_ess_cors_Nchars, ph_ess_cors_Nchars)
        rownames(ESS_correlation_table)<-c("Full Ntax", "Amniotes Ntax", "Barcoding Ntax", "mtDNA Ntax", "PhyLoTA Ntax", "Full Nchars", "Amniotes Nchars", "Barcoding Nchars", "mtDNA Nchars", "PhyLoTA Nchars")
        
    # Plot these out to pdf
    pdf(width=4, height=7, file="ESS correlations.pdf")
    par(mfrow=c(2,1))
    ## Loop to get the correlations of the paramters with the num tax and num chars
    for(i in 1:length(par_ess)){
      plot(ess_data_chars[,"ntax"], ess_data_chars[,par_ess[[i]]], cex=0.01, xlab="# taxa", ylab=paste(par_ess[[i]], "ESS", sep=" "), main=paste("r =", round(ess_cors_Ntax[[i]], 3), sep=" ")) # for each parameter, plot it against the number of taxa
      plot(ess_data_chars[,"nchar"], ess_data_chars[,par_ess[[i]]], cex=0.01, xlab="# characters", ylab=paste(par_ess[[i]], "ESS", sep=" "), main=paste("r =", round(ess_cors_Nchars[[i]], 3), sep=" ")) # for each parameter, plot it against the number of characters
    }
    dev.off()
    

    # Look at the correlations of ESS for LnL with all other params' ESS
    # First make list of non-lnl ESS parameters
    par_ess_not_lnl<-colnames(ess_data_chars)[!colnames(ess_data_chars) %in% c("analysis_name", "ntax", "nchar", "LnL", "study")]
    # Get the correlations of each of these parameters with LnL ESS
    ess_cors_ess_lnl<-sapply(par_ess_not_lnl, function(x) cor(ess_data_chars[,x], ess_data_chars[,"LnL"], use="complete.obs"))
    # Plot these out to pdf
    pdf(width=4, height=9.5, file="ESS LnL correlations.pdf")
    par(mfrow=c(2,1))
    ## Loop to get the correlations of the paramters with the num tax and num chars
    for(i in 1:length(par_ess_not_lnl)){
      plot(ess_data_chars[,par_ess_not_lnl[[i]]], ess_data_chars[,"LnL"], cex=0.01, xlab=paste(par_ess_not_lnl[[i]], "ESS", sep=" "), ylab="LnL ESS", main=paste("corr =", round(ess_cors_ess_lnl[[i]], 3), sep=" ")) # for each parameter, plot it against the number of taxa
    }
    dev.off()

    setwd(figs_ms)
    write.csv(ess_cors_ess_lnl, file="ESS_cors_LnL_ESS.csv")
    setwd(write_to)

## extract out some multi-chain diagnostics
  # Pull out ASDSF from all analyses
    asdsf_only<-sapply(diagnostics, function(x) x$ASDSF_and_SF_corr$ASDSF)
    # merge the asdsf with ESS and dataset properties
    asdsf_ess_chars<-merge(ess_data_chars, as.data.frame(asdsf_only), by.x="analysis_name", by.y="row.names")

    # Correlations of asdsf with nchars and ntax
    cor_asdsf_ntax<-cor(asdsf_ess_chars[,"ntax"], asdsf_ess_chars[,"asdsf_only"], use="complete.obs")
    cor_asdsf_nchars<-cor(asdsf_ess_chars[,"nchar"], asdsf_ess_chars[,"asdsf_only"], use="complete.obs")

    ## Check out and see if these correlations vary by dataset
      # Amniotes - LnPr gets removed out during the correlation step because it wasn't sampled for any of these analyses
        amn_asdsf_data_chars<-asdsf_ess_chars[which(asdsf_ess_chars[,"study"]=="Amniotes"),]
        amn_cor_asdsf_ntax<-cor(amn_asdsf_data_chars[,"ntax"], amn_asdsf_data_chars[,"asdsf_only"], use="complete.obs")
        amn_cor_asdsf_nchars<-cor(amn_asdsf_data_chars[,"nchar"], amn_asdsf_data_chars[,"asdsf_only"], use="complete.obs")
        plot(amn_asdsf_data_chars[,"ntax"], amn_asdsf_data_chars[,"asdsf_only"], cex=0.2, xlab="N taxa", ylab="ASDSF")
        plot(amn_asdsf_data_chars[,"nchar"], amn_asdsf_data_chars[,"asdsf_only"], cex=0.2, xlab="N chars", ylab="ASDSF")
        # Correlations are low - low variation
      # Barcoding - pinvar gets removed out during the correlation step because it wasn't sampled for any of these analyses
        bc_asdsf_data_chars<-asdsf_ess_chars[which(asdsf_ess_chars[,"study"]=="Barcoding"),]
        bc_cor_asdsf_ntax<-cor(bc_asdsf_data_chars[,"ntax"], bc_asdsf_data_chars[,"asdsf_only"], use="complete.obs")
        bc_cor_asdsf_nchars<-cor(bc_asdsf_data_chars[,"nchar"], bc_asdsf_data_chars[,"asdsf_only"], use="complete.obs")
        plot(bc_asdsf_data_chars[,"ntax"], bc_asdsf_data_chars[,"asdsf_only"], cex=0.2, xlab="N taxa", ylab="ASDSF")
        plot(bc_asdsf_data_chars[,"nchar"], bc_asdsf_data_chars[,"asdsf_only"], cex=0.2, xlab="N chars", ylab="ASDSF")
      # mtDNA
        mt_asdsf_data_chars<-asdsf_ess_chars[which(asdsf_ess_chars[,"study"]=="mtDNA"),]
        mt_cor_asdsf_ntax<-cor(mt_asdsf_data_chars[,"ntax"], mt_asdsf_data_chars[,"asdsf_only"], use="complete.obs")
        mt_cor_asdsf_nchars<-cor(mt_asdsf_data_chars[,"nchar"], mt_asdsf_data_chars[,"asdsf_only"], use="complete.obs")
        plot(mt_asdsf_data_chars[,"ntax"], mt_asdsf_data_chars[,"asdsf_only"], cex=0.2, xlab="N taxa", ylab="ASDSF")
        plot(mt_asdsf_data_chars[,"nchar"], mt_asdsf_data_chars[,"asdsf_only"], cex=0.2, xlab="N chars", ylab="ASDSF")
      # Phylota
        ph_asdsf_data_chars<-asdsf_ess_chars[which(asdsf_ess_chars[,"study"]=="PhyLoTA"),]
        ph_cor_asdsf_ntax<-cor(ph_asdsf_data_chars[,"ntax"], ph_asdsf_data_chars[,"asdsf_only"], use="complete.obs")
        ph_cor_asdsf_nchars<-cor(ph_asdsf_data_chars[,"nchar"], ph_asdsf_data_chars[,"asdsf_only"], use="complete.obs")
        plot(ph_asdsf_data_chars[,"ntax"], ph_asdsf_data_chars[,"asdsf_only"], cex=0.2, xlab="N taxa", ylab="ASDSF")
        plot(ph_asdsf_data_chars[,"nchar"], ph_asdsf_data_chars[,"asdsf_only"], cex=0.2, xlab="N chars", ylab="ASDSF")
        ## For phylota, let's try removing any asdsf > 0.02, these look like outliers
        ph_asdsf_data_chars_pruned<-ph_asdsf_data_chars[ph_asdsf_data_chars[,"asdsf_only"]<=0.02,]
        # look at the correlations without those outliers:
          ph_cor_asdsf_ntax_pruned<-cor(ph_asdsf_data_chars_pruned[,"ntax"], ph_asdsf_data_chars_pruned[,"asdsf_only"], use="complete.obs")
          ph_cor_asdsf_nchars_pruned<-cor(ph_asdsf_data_chars_pruned[,"nchar"], ph_asdsf_data_chars_pruned[,"asdsf_only"], use="complete.obs")
          plot(ph_asdsf_data_chars_pruned[,"ntax"], ph_asdsf_data_chars_pruned[,"asdsf_only"], cex=0.2, xlab="N taxa", ylab="ASDSF")
          plot(ph_asdsf_data_chars_pruned[,"nchar"], ph_asdsf_data_chars_pruned[,"asdsf_only"], cex=0.2, xlab="N chars", ylab="ASDSF")
          hist(ph_asdsf_data_chars_pruned[,"ntax"], breaks=50)
          hist(ph_asdsf_data_chars_pruned[,"nchar"], breaks=50)
          
          ## Let's see if we sample out a smaller number of chains with certain numbrs of taxa if the correlation with asdsf goes away
          ## sample out 200 from each 50 taxa window
          less_50<-sample(which(ph_asdsf_data_chars_pruned[,"ntax"]<50), 200)
          bt_50_100<-sample(which(ph_asdsf_data_chars_pruned[,"ntax"]<100 & ph_asdsf_data_chars_pruned[,"ntax"]>=50), 200)
          bt_100_150<-sample(which(ph_asdsf_data_chars_pruned[,"ntax"]<150 & ph_asdsf_data_chars_pruned[,"ntax"]>=100), 200)
          bt_150_200<-sample(which(ph_asdsf_data_chars_pruned[,"ntax"]<200 & ph_asdsf_data_chars_pruned[,"ntax"]>=150), 200)
          above_200<-sample(which(200<=ph_asdsf_data_chars_pruned[,"ntax"]), 200)
          asdsf_subsample_ph<-ph_asdsf_data_chars_pruned[c(less_50, bt_50_100, bt_100_150, bt_150_200, above_200),]
          ## See what the ASDSF/Ntax correlation is now:
          ph_cor_asdsf_ntax_pruned_subsample<-cor(asdsf_subsample_ph[,"ntax"], asdsf_subsample_ph[,"asdsf_only"], use="complete.obs")
          plot(asdsf_subsample_ph[,"ntax"], asdsf_subsample_ph[,"asdsf_only"], cex=0.2, xlab="N taxa", ylab="ASDSF")
          
          
          
        ## For Phylota, let's see what the distribution of ntax and N chars looks like for asdsf above and below 0.01
        high_asdsf_ntax<-ph_asdsf_data_chars[which(ph_asdsf_data_chars$asdsf_only>0.01),]
        median(high_asdsf_ntax$ntax)
        low_asdsf_ntax<-ph_asdsf_data_chars[which(ph_asdsf_data_chars$asdsf_only<=0.01),]
        median(low_asdsf_ntax$ntax)
        hist(high_asdsf_ntax$ntax, xlim=c(4,255))
        hist(high_asdsf_ntax$ntax, xlim=c(4,255), ylim=c(0,10000))
        hist(low_asdsf_ntax$ntax, xlim=c(4,255), add=TRUE)
        
    # Plot out the relationships to pdf
    pdf(width=4, height=7, file="ASDSF_chars-tax_cors.pdf")
    par(mfrow=c(2,1))
    plot(asdsf_ess_chars[,"ntax"], asdsf_ess_chars[,"asdsf_only"], cex=0.01, xlab="N taxa", ylab="ASDSF", main=paste("corr =", round(cor_asdsf_ntax, 3), sep=" "))
    plot(asdsf_ess_chars[,"nchar"], asdsf_ess_chars[,"asdsf_only"], cex=0.01, xlab="N chars", ylab="ASDSF", main=paste("corr =", round(cor_asdsf_nchars, 3), sep=" "))
    
    # Zoom in the y axis - there seem to be a small number of analyses with very high ASDSF - zoom in to get a better picture of what most are doing
    plot(asdsf_ess_chars[,"ntax"], asdsf_ess_chars[,"asdsf_only"], cex=0.01, xlab="N taxa", ylab="ASDSF", ylim=c(0, 0.007), main=paste("Zoomed in - corr =", round(cor_asdsf_ntax, 3), sep=" "))
    plot(asdsf_ess_chars[,"nchar"], asdsf_ess_chars[,"asdsf_only"], cex=0.01, xlab="N chars", ylab="ASDSF", ylim=c(0, 0.007), main=paste("Zoomed in - corr =", round(cor_asdsf_nchars, 3), sep=" "))
    
    dev.off()
    
  ## Do the same for split frequency correlations
    sfcors_only<-sapply(diagnostics, function(x) x$ASDSF_and_SF_corr$sf_corrs$cor) # Pull out just the split frequency correlations for each analysis
    # merge the sfcors with the dataset properties and asdsf
    sfcors_asdsf_chars<-merge(asdsf_ess_chars, as.data.frame(sfcors_only), by.x="analysis_name", by.y="row.names")
    
    # Plot out the correlations between split frequency correlations and n taxa and n characters
    cor_sfcors_ntax<-cor(sfcors_asdsf_chars[,"ntax"], sfcors_asdsf_chars[,"sfcors_only"], use="complete.obs")
    cor_sfcors_nchars<-cor(sfcors_asdsf_chars[,"nchar"], sfcors_asdsf_chars[,"sfcors_only"], use="complete.obs")
    
    # Plot out the relationships to pdf
    pdf(width=4, height=7, file="SFcors_chars-tax_cors.pdf")
    par(mfrow=c(2,1))
    plot(sfcors_asdsf_chars[,"ntax"], sfcors_asdsf_chars[,"sfcors_only"], cex=0.01, xlab="N taxa", ylab="SF cors", main=paste("corr =", round(cor_sfcors_ntax, 3), sep=" "))
    plot(sfcors_asdsf_chars[,"nchar"], sfcors_asdsf_chars[,"sfcors_only"], cex=0.01, xlab="N chars", ylab="SF cors", main=paste("corr =", round(cor_sfcors_nchars, 3), sep=" "))
    # Zoom in the y axis - there seem to be a small number of analyses with very low SF correlations - zoom in to get a better picture of what most are doing
    plot(sfcors_asdsf_chars[,"ntax"], sfcors_asdsf_chars[,"sfcors_only"], cex=0.01, xlab="N taxa", ylab="SF cors", ylim=c(0.85, 1), main=paste("Zoomed in - corr =", round(cor_sfcors_ntax, 3), sep=" "))
    plot(sfcors_asdsf_chars[,"nchar"], sfcors_asdsf_chars[,"sfcors_only"], cex=0.01, xlab="N chars", ylab="SF cors", ylim=c(0.85, 1), main=paste("Zoomed in - corr =", round(cor_sfcors_nchars, 3), sep=" "))
    
    dev.off()
    
    ## Check out correlation between asdsf and sf cor
    cor_asdsf_sdcor<-cor(sfcors_asdsf_chars[,"asdsf_only"], sfcors_asdsf_chars[,"sfcors_only"], use="complete.obs")
    plot(sfcors_asdsf_chars[,"asdsf_only"], sfcors_asdsf_chars[,"sfcors_only"], cex=0.01, xlab="ASDSF", ylab="SF cors", ylim=c(0.85, 1), xlim=c(0, 0.007), main=paste("Zoomed in - corr =", round(cor_asdsf_sdcor, 3), sep=" "))
    
  ## See if asdsf and sf cors are correlated with the topo ess
    cor_sfcors_topoESS<-cor(sfcors_asdsf_chars[,"sfcors_only"], sfcors_asdsf_chars[,"topo"], use="complete.obs")
    cor_asdsf_topoESS<-cor(sfcors_asdsf_chars[,"asdsf_only"], sfcors_asdsf_chars[,"topo"], use="complete.obs")
    # almost no correlation for either
    plot(sfcors_asdsf_chars[,"sfcors_only"], sfcors_asdsf_chars[,"topo"], cex=0.01, xlab="sf_cors", ylab="Topo ESS", main=paste("Corr =", round(cor_sfcors_topoESS, 3), sep=" "))
    plot(sfcors_asdsf_chars[,"asdsf_only"], sfcors_asdsf_chars[,"topo"], cex=0.01, xlab="ASDSF", ylab="Topo ESS", main=paste("Corr =", round(cor_asdsf_topoESS, 3), sep=" "))
    plot(sfcors_asdsf_chars[,"sfcors_only"], sfcors_asdsf_chars[,"topo"], cex=0.01, xlab="sf_cors", ylab="Topo ESS", xlim=c(0.9,1), main=paste("Zoomed in - corr =", round(cor_sfcors_topoESS, 3), sep=" "))
    plot(sfcors_asdsf_chars[,"asdsf_only"], sfcors_asdsf_chars[,"topo"], cex=0.01, xlab="ASDSF", ylab="Topo ESS", xlim=c(0,0.01), main=paste("Zoomed in - corr =", round(cor_asdsf_topoESS, 3), sep=" "))
    
## Find out which analyses have an asdsf >0.01, low split frequency correlation, 
    # or low topological ess values 
    high_asdsf<-sfcors_asdsf_chars[which(sfcors_asdsf_chars$asdsf_only>0.01),]
    low_sf_cors<-sfcors_asdsf_chars[which(sfcors_asdsf_chars$sfcors_only<0.9),] # really not sure what's an appropriate cutoff, 0.9 is probably too low, many more fail at higher cutoffs
    low_topo_ess<-sfcors_asdsf_chars[which(sfcors_asdsf_chars$topo<200),]

  # It seems that many analyses are failing to get a topological ESS value above 200, but still passing most or all other diagnostic checks
    ## What are the dataset properties of the analyses that fail topo ESS compared to those that don't:
    high_topo_ess<-sfcors_asdsf_chars[which(sfcors_asdsf_chars$topo>200),] # Make a second object that has all of the analyses that have topo ESS above 200
    
    ## Medians are probably much more useful than means
    median(low_topo_ess$ntax)
    median(high_topo_ess$ntax)
    median(sfcors_asdsf_chars$ntax) # also look at the overall median N taxa
    # N taxa median for chains that have topo ESS >200 seems to basically match overall median of all analyses, N taxa for analyses with topo ESS <200 is much higher
    
    median(low_topo_ess$nchar)
    median(high_topo_ess$nchar)
    median(sfcors_asdsf_chars$nchar) # also look at the overall median N chars
    # Looks like the median N chars for all analyses and analyses that pass are similar, n chars for failing is lower
    
    # Plot out histograms of N taxa for when topo ESS is <200 and >200
    pdf(width=6, height=9, file="Topo_ESS_failed_dataset properties.pdf")
    par(mfrow=c(2,1))
    hist(low_topo_ess$ntax, col="red", xlab="# taxa", main="Topo ESS < 200", xlim=c(0,300), ylab=NULL)
    hist(high_topo_ess$ntax, col="blue", xlab="# taxa", main="Topo ESS > 200", xlim=c(0,300), ylab=NULL)
    
    hist(low_topo_ess$nchar, col="red", xlab="# chars", main="Topo ESS < 200", xlim=c(0,5000), ylab=NULL)
    hist(high_topo_ess$nchar, col="blue", xlab="# chars", main="Topo ESS > 200", xlim=c(0,5000), ylab=NULL)
    dev.off()
    
    # How many of the chains that fail topo ESS pass LnL ESS?
    n_low_topoESS<-nrow(low_topo_ess) # total number of analyses that have low topo ESS
    n_low_topo_and_LnL_ESS<-length(low_topo_ess[which(low_topo_ess$LnL<200),][,1]) # fail topo ESS, fail LnL ESS
    n_low_LnL_ESS<-nrow(sfcors_asdsf_chars[which(sfcors_asdsf_chars$LnL<200),])
    ## Seems that the vast majority of chains that fail topo ESS have a LnL ESS > 200
    
    # Make a Venn diagram to show this
    setwd(figs_ms)
    pdf(file="Fig_5A_LnL_topo_ESS_venn.pdf", height=5, width=5)
    draw.pairwise.venn(n_low_LnL_ESS, n_low_topoESS, n_low_topo_and_LnL_ESS, category = c("LnL ESS < 200", "Topological ESS < 200"), lty = rep("blank", 
          2), fill = c("aquamarine2", "darkviolet"), alpha = rep(0.5, 2), cat.pos = c(0, 
              0), cat.dist = rep(0.025, 2), scaled = TRUE, cex=1.5, cat.cex=1.5)                     
    dev.off()
    setwd(write_to)

    # How many of the chains that fail topo ess pass asdsf
    n_low_topoESS
    n_low_topo_bad_asdsf<-nrow(low_topo_ess[which(low_topo_ess$asdsf_only>0.01),]) # fail topo ESS, fail asdsf
    n_bad_asdsf<-nrow(sfcors_asdsf_chars[which(sfcors_asdsf_chars$asdsf_only>0.01),])
    ## Make the Venn diagram for this
    setwd(figs_ms)
    pdf(file="Fig_5B_Topo_ESS_ASDSF_venn.pdf", height=5, width=5)
    draw.pairwise.venn(n_low_topoESS, n_bad_asdsf, n_low_topo_bad_asdsf, category = c("Topological ESS < 200", "ASDSF > 0.01"), lty = rep("blank", 
          2), fill = c("darkviolet", "aquamarine2"), alpha = rep(0.5, 2), cat.pos = c(0, 
             0), cat.dist = rep(0.025, 2), scaled = TRUE, cex=1.5, cat.cex=1.5)                     
    dev.off()
    setwd(write_to)
    
    # compare the split frequency correlations for the analyses that pass vs. fail topo ESS
    median(low_topo_ess$sfcors_only)
    median(high_topo_ess$sfcors_only)
    median(sfcors_asdsf_chars$sfcors_only) # also look at the overall median N chars
    # These histograms below are both zoomed in because of low frequencies of really low sf correlations
    hist(low_topo_ess$sfcors_only, col="red", xlab="SF corrs", main="Topo ESS < 200", breaks=200, xlim=c(0.95, 1))
    hist(high_topo_ess$sfcors_only, col="blue", xlab="SF corrs", main="Topo ESS > 200", breaks=200, xlim=c(0.95, 1))

    ## Better way of looking at the analyses that fail topo ESS
    pass_fail_fail_asdsf<-pass_fail_combined[which(pass_fail_combined[,"Pass ASDSF"] == FALSE),]
    # Get the proportions of these analyses that fail other diagnostics
    prop_pass_asdsf_fail<-apply(pass_fail_fail_asdsf, 2, function(x) length(which(x==TRUE))/length(x))
    prop_pass_asdsf_fail # Elevated chances of failing other diagnostics, but pass a lot of them

    
    ## Write this out as a csv table for the ms 
    prop_pass_asdsf_fail_table1<-rbind(gsub("Pass ", "", names(prop_pass_asdsf_fail)), prop_pass_asdsf_fail)
    prop_pass_asdsf_fail_table2<-gsub("ESS", "ESS > 200", prop_pass_asdsf_fail_table1)
    prop_pass_asdsf_fail_table3<-gsub("ASDSF", "ASDSF < 0.01", prop_pass_asdsf_fail_table2)
    prop_pass_asdsf_fail_table4<-gsub("SF corr", "SF corr > 0.9", prop_pass_asdsf_fail_table3)
    prop_pass_asdsf_fail_table5<-gsub("PSRF", "PSRF < 1.02", prop_pass_asdsf_fail_table4)
    prop_pass_asdsf_fail_table<-gsub("Geweke", "Geweke p > 0.05", prop_pass_asdsf_fail_table5)
    colnames(prop_pass_asdsf_fail_table)<-rep(NA, ncol(prop_pass_asdsf_fail_table))
    colnames(prop_pass_asdsf_fail_table)[[1]]<-"Proportion of chains passing each diagnostic for analyses that have ASDSF > 0.01"
    setwd(figs_ms)
    write.csv(prop_pass_asdsf_fail_table, file="Prop_Pass_FAILED_asdsf.csv")
    setwd(write_to)    
  
    ## Look at the analyses that fail topo ESS
    pass_fail_fail_topo<-pass_fail_combined[which(pass_fail_combined[,"Pass Approx topo ESS"] == FALSE),]
    # Get the proportions of these analyses that fail other diagnostics
    prop_pass_topoESS_fail<-apply(pass_fail_fail_topo, 2, function(x) length(which(x==TRUE))/length(x))
    prop_pass_topoESS_fail # Elevated chances of failing other diagnostics, but pass a lot of them
    

## Look at correlations among parameters through the chains    
    ## Pull out the parameter correlations for each individual chain
    param_cors<-unlist(lapply(diagnostics, function(x) lapply(x$Param_corrs, reshape2::melt)), recursive=FALSE)
    # Process the parameter correlations to get a dataframe with one column as the name of the two parameters
    # and the another column of the correlation coefficient - remove the duplicates that are induced by having both above/below diagonal
    param_cors_processed<-list()
    for(j in 1:length(param_cors)){
      param_cors[[j]][,1]<-as.character(param_cors[[j]][,1])
      param_cors[[j]][,2]<-as.character(param_cors[[j]][,2])
      param_cors[[j]]<-param_cors[[j]][-which(param_cors[[j]][,1]==param_cors[[j]][,2]),]  # remove the correlations of a parameter with itself
      
      names_pars_cors<-vector()
      for(i in 1:length(param_cors[[j]][,1])){
        names_pars_cors[[i]]<-paste(sort(param_cors[[j]][i,1:2])[[1]], sort(param_cors[[j]][i,1:2])[[2]], sep="_")
      }
      pars_cors_names<-cbind.data.frame(param_cors[[j]], names_pars_cors, stringsAsFactors=FALSE)
      pars_cors_names<-pars_cors_names[order(pars_cors_names[,4]),]
      pars_cors_names<-pars_cors_names[-which(duplicated(pars_cors_names[,4])),]
      param_cors_processed[[j]]<-pars_cors_names[,c(4,3)]
    }
    
    # Find the full set of parameters that were compared and correlated across all analyses
    params_cor_names<-unique(unlist(lapply(param_cors_processed, FUN=function(x) x[,1])))
    param_cors_table<-do.call(rbind, lapply(param_cors_processed, function(x) x[match(params_cor_names, x[,1]),2])) # make a new table where rows are the correlations from each chain
    colnames(param_cors_table)<-params_cor_names # assign the proper column names
    
    # Get mean and median correlations among parameters
    mean_par_cors<-apply(param_cors_table, 2, mean, na.rm=TRUE)
    median_par_cors<-apply(param_cors_table, 2, median, na.rm=TRUE)
    sd_par_cors<-apply(param_cors_table, 2, sd, na.rm=TRUE)
    # Rank these and see which correlations are strongest
    sorted_means_par_cors<-sort(mean_par_cors)
    sorted_median_par_cors<-sort(median_par_cors)
    sorted_sd_par_cors<-sort(sd_par_cors)
    sorted_means_par_cors
    sorted_median_par_cors
    sorted_sd_par_cors
    ## Check out the analyses with particularly correlated pinvar and alpha (I + G; invariants + gamma)
    alpha_pinvar_cors<-param_cors_table[,"alpha_pinvar"]
    hist(alpha_pinvar_cors)  #take a quick look at how correlated these are
    ess_pinvar_alpha<-do.call(rbind, lapply(ess_all, function(x) x[c("alpha", "pinvar")]))   # Get out ESS values for only pinvar and alpha
    colnames(ess_pinvar_alpha)<-c("alpha", "pinvar")
    ess_high_ap_cor<-ess_pinvar_alpha[which(alpha_pinvar_cors>0.6),]  # Pull out the ESS values for analyses with a pinvar/alpha correlation of >0.6
    ess_low_ap_cor<-ess_pinvar_alpha[which(alpha_pinvar_cors<0.6),]  # Pull out the ESS values for analyses with a pinvar/alpha correlation of <0.6

    # Get the means and medians
    mean(ess_low_ap_cor[,"alpha"]) 
    median(ess_low_ap_cor[,"alpha"]) 
    mean(ess_low_ap_cor[,"pinvar"])
    median(ess_low_ap_cor[,"pinvar"]) 
    
    mean(ess_high_ap_cor[,"alpha"])
    median(ess_high_ap_cor[,"alpha"])
    mean(ess_high_ap_cor[,"pinvar"])
    median(ess_high_ap_cor[,"pinvar"])
      # Basically the same

    # put the median correlations into a big matrix of values
    ## rearrange some names to make things prettier down lower
    names(median_par_cors)<-gsub("kappa_LnL","LnL_kappa", names(median_par_cors))
    names(median_par_cors)<-gsub("kappa_pi.A.","pi.A._kappa", names(median_par_cors))
    names(median_par_cors)<-gsub("kappa_pi.C.","pi.C._kappa", names(median_par_cors))
    names(median_par_cors)<-gsub("kappa_pi.G.","pi.G._kappa", names(median_par_cors))
    names(median_par_cors)<-gsub("kappa_pi.T.","pi.T._kappa", names(median_par_cors))
    names(median_par_cors)<-gsub("kappa_pinvar","pinvar_kappa", names(median_par_cors))
    names(median_par_cors)<-gsub("kappa_TL","TL_kappa", names(median_par_cors))

    pars_cors_mat_names<-unique(unlist(strsplit(names(median_par_cors), "_")))  # get the names of all the parameters in this 
    pars_cors_mat<-matrix(NA, 16, 16) # make a blank matrix
    colnames(pars_cors_mat)<-pars_cors_mat_names
    rownames(pars_cors_mat)<-pars_cors_mat_names

    for(i in rownames(pars_cors_mat)){ # nested loops to replace each element of the matrix with the corresponding entry from median_par_cors
      for(j in colnames(pars_cors_mat)){
        if(i=="LnPr"){ # this if is to get things looking nice and triangular - switching around i and j to pull out values
          pars_cors_mat[j,i]<-median_par_cors[paste(i,j,sep="_")]
        }else{
          pars_cors_mat[i,j]<-median_par_cors[paste(i,j,sep="_")]
        }
      }
    }
    ## Not sure why this is necessary, but things aren't all on top of the diagonal otherwise
    pars_cors_mat["kappa", "LnPr"]<-median_par_cors["kappa_LnPr"]
    pars_cors_mat["alpha", "LnPr"]<-median_par_cors["alpha_LnPr"] 
    pars_cors_mat["LnL", "LnPr"]<-median_par_cors["LnL_LnPr"] 
    pars_cors_mat
    
    ## Make this into a heat map - this uses code from here: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
    melted_pars_cormat <- reshape2::melt(round(pars_cors_mat, 2))
    # Create a ggheatmap
    ggheatmap <- ggplot(melted_pars_cormat, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value = "gray90",
                           midpoint = 0, limit = c(-1,1), space = "Lab", 
                           name="Pearson\nCorrelation") +
      theme_minimal()+ # minimal theme
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 12, hjust = 1))+
      coord_fixed()
    # Print the heatmap
    print(ggheatmap)
    # Write it to pdf with the actual correlation values
    setwd(figs_ms)
    pdf(file="pars_corrs.pdf", height=8, width=8)
    ggheatmap + 
      geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.6, 0.7),
        legend.direction = "horizontal",
        legend.title=element_text(size=18),
        legend.text=element_text(size=15))+
      guides(fill = guide_colorbar(barwidth = 12, barheight = 2, title.position = "top", title.hjust = 0.6))
    
    dev.off()
    setwd(write_to)
    
## Look at chains that fail either ASDSF or PSRF and see how they do for single-chian diags
    fail_ASDSForPSRF<-apply(pass_fail_combined, 1, function(x) x["Pass ASDSF"]==FALSE | x["Pass PSRF"]==FALSE) # see which analyses failed either ASDSF or PSRF
    pass_fail_fail_multi<-pass_fail_combined[fail_ASDSForPSRF,]  # pull out the pass/fail table for just the analyses that failed ASDSF or PSRF
    prop_pass_fail_multi<-apply(pass_fail_fail_multi, 2, function(x) length(which(x==TRUE))/length(x)) # get out the proportion of analyses that pass other diags that fail either of the multichain diags
    prop_pass_fail_multi    
    
    ## Make a Venn diagram of how many analyses fail ASDSF or PSRF
    n_fail_asdsf<-length(which(pass_fail_combined[,"Pass ASDSF"]==FALSE))
    n_fail_psrf<-length(which(pass_fail_combined[,"Pass PSRF"]==FALSE))
    n_fail_psrf_and_asdsf<-length(which(apply(pass_fail_combined, 1, function(x) x["Pass ASDSF"]==FALSE & x["Pass PSRF"]==FALSE)==TRUE))
    setwd(figs_ms)
    pdf(file="Fig_5C_ASDSF_PSRF_venn.pdf", height=5, width=5)
    draw.pairwise.venn(n_fail_asdsf, n_fail_psrf, n_fail_psrf_and_asdsf, category = c("ASDSF > 0.01", "PRSF > 1.02"), lty = rep("blank", 
        2), fill = c("aquamarine2", "darkviolet"), alpha = rep(0.5, 2), cat.pos = c(0, 
            0), cat.dist = rep(0.025, 2), scaled = TRUE, cex=1.5, cat.cex=1.5)                     
    dev.off()
    setwd(write_to)
    
    ## make a Venn diagram for how many analyses fail either ASDSF or PSRF compared to that those that fail ESS for any parameters
    n_fail_asdsfORpsrf<-length(which(fail_ASDSForPSRF==TRUE))
    n_fail_any_ESS<-length(which(pass_fail_combined[,"Pass ESS"]==FALSE))
    n_fail_adsfORpsrfANDanyESS<-length(which(pass_fail_fail_multi[,"Pass ESS"]==FALSE))
    setwd(figs_ms)
    pdf(file="Fig_5D_multi_chain_vs_ESS_venn.pdf", height=5, width=5)
    draw.pairwise.venn(n_fail_asdsfORpsrf, n_fail_any_ESS, n_fail_adsfORpsrfANDanyESS, category = c("ASDSF > 0.01 OR PRSF > 1.02", "Any non-topo ESS < 200"), lty = rep("blank", 
        2), fill = c("darkviolet", "aquamarine2"), alpha = rep(0.5, 2), cat.pos = c(0, 
          0), cat.dist = rep(0.025, 2), scaled = TRUE, cex=1.5, cat.cex=1.5)                     
    dev.off()
    setwd(write_to)
    

## Get the acceptance rates for each move:
acceptances<-lapply(diagnostics, function(x) x$acc_rates)
# get only the analyses that successfully had acceptance rates pulled out - not all analyses had the file necessary to pull out acceptance rates
acceptances<-acceptances[which(!lapply(acceptances, function(x) x[[1]])=="NA")]


## Different analyses included different moves, need to find the full list of moves that were used
move_names<-unique(unlist(lapply(acceptances, FUN=function(x) colnames(x))))
# Use this to order all parameters the same and make a table of gens to ess 200 for all for each chain
acceptance_table<-do.call(rbind, lapply(acceptances, function(x) x[,match(move_names, colnames(x))]))
colnames(acceptance_table)<-move_names

# Get mean, median, and sd acceptance for each move
mean_acc<-apply(acceptance_table, 2, mean, na.rm=TRUE)
median_acc<-apply(acceptance_table, 2, median, na.rm=TRUE)
sd_gens_acc<-apply(acceptance_table, 2, sd, na.rm=TRUE)
range_acc<-apply(acceptance_table, 2, range, na.rm=TRUE)
# Rank these parameters and see what converges fastest
sorted_means_acc<-sort(mean_acc)
sorted_median_acc<-sort(median_acc)
sorted_sd_acc<-sort(sd_gens_acc)
sorted_means_acc
sorted_median_acc
sorted_sd_acc
range_acc

# Histograms of rates for topology moves
hist(acceptance_table[,"ExtSPR(Tau,V)"])
hist(acceptance_table[,"ExtTBR(Tau,V)"])
hist(acceptance_table[,"NNI(Tau,V)"])
hist(acceptance_table[,"ParsSPR(Tau,V)"])


## make a ridgeline plot of the acceptance rates for parameters
# First transform the dataframe into a format easily readable by ggplot
accept_melted<-reshape2::melt(acceptance_table, id.vars = NULL)
# Not every parameter is in every analysis, so need to prune out the NA's
accept_melted<-accept_melted[which(is.na(accept_melted$value)==FALSE),]
setwd(figs_ms)
pdf(file="Fig_4_Move_acceptance_ridges.pdf", width=6, height=4)
ggplot(accept_melted, aes(y=Var2, x=value, fill = Var2)) +
  labs(title=NULL,y="Move", x = "Acceptance %")+
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, scale=4)+ 
  scale_fill_cyclical(values = c("aquamarine2", "darkviolet"))+
  scale_y_discrete(expand = expand_scale(add = c(0.2, 3.3)))
dev.off()
setwd(write_to)

  
############################################################################### 
############################################################################### 
############################################################################### 
######    Look at tree length, chain swap rates, and credible size info
############################################################################### 
  
  ## Not all analyses had TL, chain swaps, or cred size that could be extracted -- prune down lists of diagnostics that have these only
  # Tree length first
  failed_TL<-which(sapply(diagnostics, function(x) "FAILED" %in% x$TL[[1]]))  #which analyses didn't pull out tree length?
  if(length(failed_TL)>0){  # if any don't have TL, remove them and get list of diags that have TL only
    diag_TL<-diagnostics[-failed_TL]  
  }else{
    diag_TL<-diagnostics
  }
  # Same for credible set size
  failed_cred_size<-which(sapply(diagnostics, function(x) is.na(x$cred_size[[1]])))
  if(length(failed_cred_size)>0){ 
    diag_cred<-diagnostics[-failed_cred_size]  
  }else{
    diag_cred<-diagnostics
  }
  # Same for chains swaps size
  failed_cswaps<-which(sapply(diagnostics, function(x) "No log file found" %in% x$chain_swaps[[2]]))
  if(length(failed_cswaps)>0){ 
    diag_cswaps<-diagnostics[-failed_cswaps]  
  }else{
    diag_cswaps<-diagnostics
  }
  
  ############## Correlations of TL
  ## Need to add in chain information to ess_data_chars
  ess_chars_by_chain_pre<-ess_data_chars # make new object
  ess_chars_by_chain_pre$analysis_name<-paste(ess_data_chars$analysis_name,  rep(1:2, length.out=length(ess_data_chars$analysis_name)), sep=".") # paste .1 and .2 onto each analysis name to indicate the chain
  colnames(ess_chars_by_chain_pre)[which(colnames(ess_chars_by_chain_pre)=="analysis_name")]<-"chain_name"
  ess_chars_by_chain<-cbind(ess_data_chars[,"analysis_name"], ess_chars_by_chain_pre)
  colnames(ess_chars_by_chain)[[1]]<-"analysis_name"
  
  ## Set up TL to be formatted for the correlations and merge it with ESS, etc. and then merge them together
  TL<-lapply(diag_TL, function(x) x$TL[which(names(x$TL)=="median")]) # get out just the median TL
  TL_cor_form<-cbind.data.frame(gsub(".median", "", names(unlist(TL))),  unlist(TL), stringsAsFactors=FALSE)  # make this into a dataframe with the analysis name as a column
  colnames(TL_cor_form)<-c("analysis_name", "med_TL") # add in column names
  TL_cor_form$analysis_name<-paste(gsub(".median", "", names(unlist(TL))), rep(1:2, length.out=length(TL_cor_form$analysis_name)), sep=".")
  ess_data_TL<-merge(ess_chars_by_chain, TL_cor_form, by.x="chain_name", by.y="analysis_name", all.x=TRUE)
  # Add asdsf into this
  asdsf_for_TL<-sapply(diagnostics, function(x) x$ASDSF_and_SF_corr$ASDSF) # Pull out ASDSF from all alanyses
  # merge the asdsf with ESS & other data characteristics
  ESS_asdsf_TL<-merge(ess_data_TL, as.data.frame(asdsf_for_TL), by.x="analysis_name", by.y="row.names", all.x=TRUE)
  
  # Make a list of the different parameters that have ESS values
  par_ess_TL<-colnames(ESS_asdsf_TL)[!colnames(ESS_asdsf_TL) %in% c("analysis_name", "med_TL", "chain_name", "ntax", "nchar", "study")]
  # Get the correlations of each of these parameters with TL
  ess_cors_TL<-sapply(par_ess_TL, function(x) cor(ESS_asdsf_TL[,"med_TL"], ESS_asdsf_TL[,x], use="complete.obs"))
  # Plot these out to pdf
  pdf(width=4, height=4, file="ESS TL correlations.pdf")
  ## Loop to get the correlations of the paramters with the num tax and num chars
  for(i in 1:length(par_ess_TL)){
    plot(ESS_asdsf_TL[,"med_TL"], ESS_asdsf_TL[,par_ess_TL[[i]]], cex=0.01, xlab="median TL", ylab=paste(par_ess_TL[[i]], "ESS", sep=" "), main=paste("r =", round(ess_cors_TL[[i]], 3), sep=" ")) # for each parameter, plot it against the TL
  }
  dev.off()
  
  
  ############## Correlations of credible set size
  ## Set up cred formatted in a similar way as TL for merge with ESS, etc
  cred<-lapply(diag_cred, function(x) x$cred_size)
  cred_cor_form<-cbind.data.frame(names(cred), as.numeric(cred), stringsAsFactors = FALSE)
  colnames(cred_cor_form)<-c("analysis_name", "cred_size")
  ess_data_cred<-merge(ESS_asdsf_TL, cred_cor_form, by.x="analysis_name", by.y="analysis_name", all.x=TRUE) # merge together the ESS and credible set size
  # Make a list of the different parameters that have ESS values, as well as the asdsf
  par_ess_cred<-colnames(ess_data_cred)[!colnames(ess_data_cred) %in% c("chain_name", "analysis_name", "cred_size", "ntax", "nchar", "study", "med_TL")]
  # Get the correlations of each of these parameters with credible set size
  ess_asdsf_cors_cred<-sapply(par_ess_cred, function(x) cor(as.numeric(ess_data_cred[,"cred_size"]), as.numeric(ess_data_cred[,x]), use="complete.obs"))
  
  # Plot these out to pdf
  pdf(width=4, height=4, file="ESS Credible Set correlations.pdf")
  ## Loop to get the correlations of the paramters with the num tax and num chars
  for(i in 1:length(par_ess_cred)){
    plot(ess_data_cred[,"cred_size"], ess_data_cred[,par_ess_cred[[i]]], cex=0.01, xlab="Credible set size", ylab=paste(par_ess_cred[[i]], "ESS", sep=" "), main=paste("r =", round(ess_cors_cred[[i]], 3), sep=" ")) # for each parameter, plot it against the credible set size
  }
  dev.off()
  
  ## Compare cred set size to the number of taxa - add ntax into the asdsf and cred size object
  cor(ess_data_cred[,"ntax"], ess_data_cred[,"cred_size"], use="complete.obs")
  plot(ess_data_cred[,"ntax"], ess_data_cred[,"cred_size"], cex=0.01)
        ## The above looks horrendous - split it up by dataset and see if it looks any better
        # Amniotes
        amn_asdsf_cred_ntax<-ess_data_cred[which(ess_data_cred[,"study"]=="Amniotes"),]
        amn_cor_ntax_cred<-cor(amn_asdsf_cred_ntax[,"ntax"], amn_asdsf_cred_ntax[,"cred_size"])
        plot(amn_asdsf_cred_ntax[,"ntax"], amn_asdsf_cred_ntax[,"cred_size"], cex=0.2, xlab="N tax", ylab="Cred set")
        
        amn_cor_cred_asdsf<-cor(amn_asdsf_cred_ntax[,"cred_size"], amn_asdsf_cred_ntax[,"asdsf_for_TL"])
        plot(amn_asdsf_cred_ntax[,"cred_size"], amn_asdsf_cred_ntax[,"asdsf_for_TL"], cex=0.2, xlab="Cred set", ylab="asdsf")
        
        ## DO NOT HAVE CRED SIZE FOR mtDNA or Barcoding
        
        ## Phylota
        ph_asdsf_cred_ntax<-ess_data_cred[which(ess_data_cred[,"study"]=="PhyLoTA"),]
        ph_cor_ntax_cred<-cor(ph_asdsf_cred_ntax[,"ntax"], ph_asdsf_cred_ntax[,"cred_size"])
        plot(ph_asdsf_cred_ntax[,"ntax"], ph_asdsf_cred_ntax[,"cred_size"], cex=0.01, xlab="N tax", ylab="Cred set")
        # compare asdsf and cred size for Phylota only
        ph_cor_ntax_cred<-cor(ph_asdsf_cred_ntax[,"cred_size"], ph_asdsf_cred_ntax[,"asdsf_for_TL"])
        plot(ph_asdsf_cred_ntax[,"cred_size"], ph_asdsf_cred_ntax[,"asdsf_for_TL"], cex=0.01, xlab="Cred size", ylab="ASDSF")
            # Take a look with any asdsf > 0.02 removed
            ph_asdsf_cred_ntax_pruned<-ph_asdsf_cred_ntax[which(ph_asdsf_cred_ntax$asdsf_for_TL<=0.02),]
            ph_cor_cred_asdsf_pruned<-cor(ph_asdsf_cred_ntax_pruned[,"cred_size"], ph_asdsf_cred_ntax_pruned[,"asdsf_for_TL"])
            plot(ph_asdsf_cred_ntax_pruned[,"cred_size"], ph_asdsf_cred_ntax_pruned[,"asdsf_for_TL"], cex=0.01, xlab="Cred size", ylab="asdsf")
            # With asdsf pruned this way, let's try looking at asdsf, cred size, and number of taxa all together in a 3D plot
              # scatterplot3d(x=ph_asdsf_cred_ntax_pruned[,"ntax"], y=ph_asdsf_cred_ntax_pruned[,"cred_size"], z=ph_asdsf_cred_ntax_pruned[,"asdsf_for_TL"], cex.symbols=0.01)
              # plot3d(x=ph_asdsf_cred_ntax_pruned[,"ntax"], y=ph_asdsf_cred_ntax_pruned[,"cred_size"], z=ph_asdsf_cred_ntax_pruned[,"asdsf_for_TL"], size=0.01, xlab="ntax", ylab="cred_size", zlab="ASDSF")
            # Plot just cred size against number of taxa
            ph_cor_ntax_cred_pruned<-cor(ph_asdsf_cred_ntax_pruned[,"ntax"], ph_asdsf_cred_ntax_pruned[,"cred_size"])
            plot(ph_asdsf_cred_ntax_pruned[,"ntax"], ph_asdsf_cred_ntax_pruned[,"cred_size"], cex=0.01, xlab="N Tax", ylab="Cred size")
              ## This quickly hits the maximum cred size for the chain length/sampling interval
  
  
  ############## Correlations of chain swaps
  ## Set up cswaps formatted in a similar way as cred set, TL for merge with ESS, etc
  cswaps<-lapply(diag_cswaps, function(x) x$chain_swaps) # get out chain swaps
  cswaps12<-unlist(lapply(cswaps, function(x) lapply(x, function(y) y[1,2])))  ## Get just the swap acceptance between chains 1 and 2
  cswap_cor_form<-cbind.data.frame(names(cswaps12), as.numeric(cswaps12), stringsAsFactors = FALSE)
  colnames(cswap_cor_form)<-c("analysis_name", "cswap_12")
  cswap_cor_form$analysis_name<-paste(gsub(".$", "", names(cswaps12)), rep(1:2, length.out=length(cswap_cor_form$analysis_name)), sep=".")
  ess_data_cswap<-merge(ess_data_cred, cswap_cor_form, by.x="chain_name", by.y="analysis_name", all.x=TRUE) # finally, merge together the ESS and cswap data and everything else
  # Make a list of the different parameters that have ESS values
  par_ess_cswap<-colnames(ess_data_cswap)[!colnames(ess_data_cswap) %in% c("analysis_name", "chain_name", "cswap_12", "cred_size", "ntax", "nchar", "study", "med_TL")]
  # Get the correlations of each of these parameters with chain swaps
  ess_cors_cswap<-sapply(par_ess_cswap, function(x) cor(as.numeric(ess_data_cswap[,"cswap_12"]), as.numeric(ess_data_cswap[,x]), use="complete.obs"))
  # Plot these out to pdf
  pdf(width=4, height=4, file="ESS Chain 1&2 swap correlations.pdf")
  ## Loop to get the correlations of the paramters with the num tax and num chars
  for(i in 1:length(par_ess_cswap)){
    plot(ess_data_cswap[,"cswap_12"], ess_data_cswap[,par_ess_TL[[i]]], cex=0.01, xlab="Chain 1&2 swap acceptance", ylab=paste(par_ess_cswap[[i]], "ESS", sep=" "), main=paste("r =", round(ess_cors_cswap[[i]], 3), sep=" "))
  }
  dev.off()
  
  
  ### Look at TL/number of branches rather than raw TL - this is probably a more informative measure
  tl_branch<-ess_data_cswap[,"med_TL"]/((ess_data_cswap[,"ntax"]-1)*2)  # divide the tree length by (#taxa-1)*2  (this is the number of branches)
  ess_data_tl_branch<-cbind(ess_data_cswap, tl_branch)
  # Make a list of the different parameters to be correlated
  par_ess_TL_branch<-colnames(ess_data_tl_branch)[!colnames(ess_data_tl_branch) %in% c("analysis_name", "chain_name", "cswap_12", "cred_size", "ntax", "nchar", "study", "med_TL", "tl_branch")]
  # Get the correlations of each of these parameters with TL
  ess_cors_TL_branch<-sapply(par_ess_TL_branch, function(x) cor(ess_data_tl_branch[,"tl_branch"], ess_data_tl_branch[,x], use="complete.obs"))
  # Plot these out to pdf
  pdf(width=4, height=4, file="ESS TLperBranch correlations.pdf")
  ## Loop to get the correlations of the paramters with the num tax and num chars
  for(i in 1:length(par_ess_TL_branch)){
    plot(ess_data_tl_branch[,"tl_branch"], ess_data_tl_branch[,par_ess_TL_branch[[i]]], cex=0.01, xlab="TL/n branches", ylab=paste(par_ess_TL_branch[[i]], "ESS", sep=" "), main=paste("r =", round(ess_cors_TL_branch[[i]], 3), sep=" ")) # for each parameter, plot it against the TL
  }
  dev.off()
  
  ## Add these into the table of the other correlations with parameter value ESS
  # Combine the chain swap/parameter ESS correlations into the table with the correlations among parameter ESS and ntax and nchars
  ESS_corr_table_wswaps<-bind_rows(ESS_correlation_table, ess_cors_cswap)
  ESS_corr_table_full<-bind_rows(ESS_corr_table_wswaps, ess_cors_TL_branch)
  ESS_corr_table_full<-as.matrix(ESS_corr_table_full)
  rownames(ESS_corr_table_full)<-c("# taxa", "# taxa Amniotes", "# taxa Barcoding", "# taxa mtDNA", "# taxa PhyLoTA", "# chars", "# chars Amniotes", "# chars Barcoding", "# chars mtDNA", "# chars PhyLoTA", "Chain 1&2 swap rate", "TL/branch")

  
  ## This above table has ASDSF for some rows, need to add in for the rest 
  ESS_corr_table_full["# taxa", "asdsf_for_TL"]<-cor_asdsf_ntax
  ESS_corr_table_full["# taxa Amniotes", "asdsf_for_TL"]<-amn_cor_asdsf_ntax
  ESS_corr_table_full["# taxa Barcoding", "asdsf_for_TL"]<-bc_cor_asdsf_ntax
  ESS_corr_table_full["# taxa mtDNA", "asdsf_for_TL"]<-mt_cor_asdsf_ntax
  ESS_corr_table_full["# taxa PhyLoTA", "asdsf_for_TL"]<-ph_cor_asdsf_ntax
  ESS_corr_table_full["# chars", "asdsf_for_TL"]<-cor_asdsf_nchars
  ESS_corr_table_full["# chars Amniotes", "asdsf_for_TL"]<-amn_cor_asdsf_nchars
  ESS_corr_table_full["# chars Barcoding", "asdsf_for_TL"]<-bc_cor_asdsf_nchars
  ESS_corr_table_full["# chars mtDNA", "asdsf_for_TL"]<-mt_cor_asdsf_nchars
  ESS_corr_table_full["# chars PhyLoTA", "asdsf_for_TL"]<-ph_cor_asdsf_nchars
  
  colnames(ESS_corr_table_full)<-gsub("asdsf_for_TL", "ASDSF", colnames(ESS_corr_table_full))
  
  setwd(figs_ms)
  write.csv(ESS_corr_table_full, file="Table_3_ESS_ASDSF_correlations.csv")
  setwd(write_to)
  
  
  
  ## Make a table that has other correlations with ASDSF and one with the correlations with sfcors
  otherASDSF_cors<-c(ph_cor_asdsf_ntax_pruned, ph_cor_asdsf_nchars_pruned, cor_asdsf_topoESS)
  names(otherASDSF_cors)<-c( "# taxa PhyLoTA ASDSF =< 0.02", "# chars PhyLoTA ASDSF =< 0.02", "Topo ESS")
  sfcors_cors<-c(cor_sfcors_ntax, cor_sfcors_nchars, cor_sfcors_topoESS, cor_asdsf_sdcor)
  names(sfcors_cors)<-c("# taxa", "# chars", "Topo ESS", "ASDSF")
  # write these out
  setwd(figs_ms)
  write.csv(otherASDSF_cors, file="ASDSF_Other_corrs.csv")
  write.csv(sfcors_cors, file="SFcors_corrs.csv")
  setwd(write_to)  
  
#### Look at relationship between cred set size and chain swap acceptance rate
  cor_cred_cswaps<-cor(ess_data_cswap[,"cswap_12"], ess_data_cswap[,"cred_size"], use="complete.obs")
  # Plot out the relationships to pdf
  pdf(width=4, height=4, file="Cred_set_Chains_swap_cors.pdf")
  plot(ess_data_cswap[,"cswap_12"], ess_data_cswap[,"cred_size"], cex=0.01, xlab="Chain 1&2 acceptance", ylab="Cred set size", main=paste("corr =", round(cor_cred_cswaps, 3), sep=" "))
  dev.off()
  
#####################################################################################  
#####################################################################################  
###      Get out 1,000 of each type of analysis that we want to reestimate
#####################################################################################
#####################################################################################    
  #### Find analyses that fail Geweke's and don't fail other diags ####
  # start by making an object that contains relatively common diags for the analyses that fail Geweke
  pass_fail_common_fg<-pass_fail_combined[names(failed_geweke),c("Pass ESS", "Pass ASDSF", "Pass PSRF", "Pass Approx topo ESS")]
  # Then identify which of these pass all of the other common diags
  fail_g_pass_common<-rownames(pass_fail_common_fg[apply(pass_fail_common_fg, 1, function(x) !(FALSE %in% x)),])
  fail_g_pass_common  # these are the names of chains that fail Geweke's but pass ESS, ASDSF, PSRF, and Approx topo ESS
  
  # Which chains that fail topo ESS pass LnL ESS?
  fail_topo_good_lnlESS<-as.character(low_topo_ess[which(low_topo_ess$LnL>200),][,1]) # fail topo ESS, pass LnL ESS
  fail_topo_bad_lnlESS<-as.character(low_topo_ess[which(low_topo_ess$LnL<200),][,1]) # fail topo ESS, fail LnL ESS
  
  # Find chains that fail ESS and see what's up in some of these
  fail_oneormore_ESS<-names(pass_fail_combined[,"Pass ESS"][pass_fail_combined[,"Pass ESS"]==FALSE])
  
  
## List these out for ease of use:
  fail_g_pass_common # Fail Geweke's but pass common diagnostics
  fail_topo_good_lnlESS # Fail topo but have good LnL ESS
  fail_topo_bad_lnlESS # Fail topo and LnL ESS
  fail_oneormore_ESS # Fail one of more parameter ESS

  ## Right now these above all refer to single chains, and so an analysis could be in each of these lists twice if both chains failed, to make sure that we have 1000 analyses to reanalyses, get the base analysis name and only the unique ones
  fail_g_pass_common<-gsub("\\.nex.*$|\\.run.*$", "", fail_g_pass_common)  # strip off the run info
  fail_g_pass_common<-unique(fail_g_pass_common)  # get only the unique entries here
  
  fail_topo_good_lnlESS<-gsub("\\.nex.*$|\\.run.*$", "", fail_topo_good_lnlESS)  # strip off the run info
  fail_topo_good_lnlESS<-unique(fail_topo_good_lnlESS)  # get only the unique entries here
  
  fail_topo_bad_lnlESS<-gsub("\\.nex.*$|\\.run.*$", "", fail_topo_bad_lnlESS)  # strip off the run info
  fail_topo_bad_lnlESS<-unique(fail_topo_bad_lnlESS)  # get only the unique entries here

  fail_oneormore_ESS<-gsub("\\.nex.*$|\\.run.*$", "", fail_oneormore_ESS)  # strip off the run info
  fail_oneormore_ESS<-unique(fail_oneormore_ESS)  # get only the unique entries here
  
    
  # Sample 1,000 analyses out from each of these for each that have more than 1000 analyses - if not, then just use all
  if(length(fail_g_pass_common)>1000){
    set.seed(77)
    fail_g_pass_com_LARGE<-sample(fail_g_pass_common, 1000)
  }else{
    fail_g_pass_com_LARGE<-fail_g_pass_common
  }
  
  if(length(fail_topo_good_lnlESS)>1000){
    set.seed(77)
    fail_topo_good_lnlESS_LARGE<-sample(fail_topo_good_lnlESS, 1000)
  }else{
    fail_topo_good_lnlESS_LARGE<-fail_topo_good_lnlESS
  }
  
  if(length(fail_topo_bad_lnlESS)>1000){
    set.seed(77)
    fail_topo_bad_lnlESS_LARGE<-sample(fail_topo_bad_lnlESS, 1000)
  }else{
    fail_topo_bad_lnlESS_LARGE<-fail_topo_bad_lnlESS
  }

  if(length(fail_oneormore_ESS)>1000){
    set.seed(77)
    fail_oneormore_ESS_LARGE<-sample(fail_oneormore_ESS, 1000)    
  }else{
    fail_oneormore_ESS_LARGE<-fail_oneormore_ESS
  }

  
  # Write a csv file containing the lines for each of these
  ## This has already been done, so commented out
    # write.csv(fail_g_pass_com_LARGE, file="fail_Geweke_pass_common_LARGE.csv")
    # write.csv(fail_topo_good_lnlESS_LARGE, file="fail_topo_good_LnL_ESS_LARGE.csv")
    # write.csv(fail_topo_bad_lnlESS, file="fail_topo_bad_LnL_ESS_LARGE.csv")
    # write.csv(fail_oneormore_ESS_LARGE, file="fail_oneormore_ESS_LARGE.csv")
  
  ## Here, we'll read in the csv's
  set_1<-read.csv("fail_Geweke_pass_common_LARGE.csv", stringsAsFactors = FALSE)
  set_2<-read.csv("fail_oneormore_ESS_LARGE.csv", stringsAsFactors = FALSE)
  set_3<-read.csv("fail_topo_bad_LnL_ESS_LARGE.csv", stringsAsFactors = FALSE)
  set_4<-read.csv("fail_topo_good_LnL_ESS_LARGE.csv", stringsAsFactors = FALSE)
  
  
  
  
  ## Write out an R data file that has the diagnostic info on these subsets of analyses
  all_to_reestimate<-c(set_1[,2], set_2[,2], set_3[,2], set_4[,2])
  orig_diags_subsets_LARGE<-diagnostics[unique(all_to_reestimate)]  # the original diagnostics for the datasets that we reestimated trees for with different heating and with nst=mixed
  # get out the pass/fail table for these same sets of analyses
  to_pull_LARGE<-which(gsub("\\.nex.*$|\\.run.*$", "", rownames(pass_fail_combined)) %in% unique(all_to_reestimate))  # get out the indices of the rows to pull out of the combined pass/fail table
  pass_fail_orig_reestimated_LARGE<-pass_fail_combined[to_pull_LARGE,]
  # Identify parameters that fail Geweke's in each chain
  failed_geweke_orig_reest_pre_LARGE<-unlist(lapply(orig_diags_subsets_LARGE, function(x) x$Params_fail_Geweke), recursive=FALSE)
  failed_geweke_orig_reest_LARGE<-lapply(failed_geweke_orig_reest_pre_LARGE, names)   # get the actual names of these parameters as the value in the object
  # Identify parameters that failed ESS in each chain
  failed_ESS_orig_reest_pre_LARGE<-unlist(lapply(orig_diags_subsets_LARGE, function(x) x$Params_ESS_less_200), recursive=FALSE)
  failed_ESS_orig_reest_LARGE<-lapply(failed_ESS_orig_reest_pre_LARGE, names)   # get the actual names of these parameters as the value in the object
  
  # save these as a list
  orig_diags_LARGE<-list(pass_fail_orig_reestimated_LARGE, failed_geweke_orig_reest_LARGE, failed_ESS_orig_reest_LARGE, orig_diags_subsets_LARGE)
  names(orig_diags_LARGE)<-c("pass_fail", "pars_fail_Geweke", "Pars_failed_ESS", "raw_diags")
  save(orig_diags_LARGE, file="Orig_diags_for_reestimated_LARGE.RData")
  
    ## Get the chain swap acceptance for each of these so that we can adjust heating accordingly
    cswaps_for_reest<-cswap_cor_form[which(gsub("\\.1$|\\.2$", "",cswap_cor_form$analysis_name) %in% all_to_reestimate),]
    cswap_cor_form$analysis_name<-gsub("\\.1$|\\.2$", "",cswap_cor_form$analysis_name)
    # This so far includes a separate chain swap acceptance rate for each chain, just get the mean of these for each analysis
    cswap_cor_form[cswap_cor_form$analysis_name==test, 2]
    cswap_mean_reest<-sapply(unique(cswap_cor_form$analysis_name), function(x) mean(cswap_cor_form[cswap_cor_form$analysis_name==x, 2]))
    names(cswap_mean_reest)<-unique(cswap_cor_form$analysis_name)
    ## specify if we want to increase or decrease heating for each analysis - this is coarse, but let's decrease if less than or equal to 0.5 and increase otherwise
    heating_change<-cswap_mean_reest
    heating_change[cswap_mean_reest>0.5]<-"increase"
    heating_change[cswap_mean_reest<=0.5]<-"decrease"
    save(cswap_mean_reest, file="example_chains_cswaps_LARGE.RData")
    save(heating_change, file="example_chains_heatChange.RData")
  
    # Get the burnin for each chain
    example_chains_burn_LARGE<-lapply(all_to_reestimate, function(x) max(unlist(diagnostics[[x]]$Burnin)))
    names(example_chains_burn_LARGE)<-all_to_reestimate
    save(example_chains_burn_LARGE, file="example_chains_burnin_LARGE.RData")
    
## Same thing for anlayses with I+G and high correlations between I and G for reestimation with only G  
  ## Pull out studies with a high I+G correlation that either have bad ESS values or that have all good ESS - run these through without +I and see what changes
  # already have an object from up above of the analyses that have a correlation of >0.6 between I and G
  all_high_IG_cor_bad_I_ESS<-rownames(ess_high_ap_cor[which(ess_high_ap_cor[,"alpha"]<200),])  # get out the analyses that have high correlations between IG and also ESS values <200 for G
  all_high_IG_cor_bad_G_ESS<-rownames(ess_high_ap_cor[which(ess_high_ap_cor[,"pinvar"]<200),])  # get out the analyses that have high correlations between IG and also ESS values <200 for I
  all_high_IG_cor_bad_ESS<-c(all_high_IG_cor_bad_I_ESS, all_high_IG_cor_bad_G_ESS)
  
  pass_fail_common_diags<-pass_fail_combined[,c("Pass ESS", "Pass ASDSF", "Pass PSRF", "Pass Approx topo ESS")]  # make an object that contains the pass/fail stats for the common diagnostics (e.g., exclude Geweke's, pseudo ESS, SF Correlations)
  pass_common_diags<-names(which(apply(pass_fail_common_diags, 1, function(x) !(FALSE %in% x))))  # Get a list of the analyses that pass all common diags 
  high_IG_cor_good_conv_both_chains<-rownames(ess_high_ap_cor)[which(rownames(ess_high_ap_cor) %in% pass_common_diags)] # which of the analyses that converge well (as far as we can tell) have high IG correlations
  high_IG_cor_good_conv_both_chains<-gsub("\\.nex.*$|\\.run.*$", "", high_IG_cor_good_conv_both_chains) # strip off file extension stuff
  high_IG_cor_good_conv<-high_IG_cor_good_conv_both_chains[which(duplicated(high_IG_cor_good_conv_both_chains))] ## high_IG_cor_good_conv_both_chains looks at each chain independently, i.e., if one chain passes all ESS, but the other doesn't, then it can be included - instead remove any analyses that don't pass common diagnostics for both chains
  
  ## Strip off run/file extensions and get unique entries only for all_high_IG_cor_bad_ESS
  all_high_IG_cor_bad_ESS<-gsub("\\.nex.*$|\\.run.*$", "", all_high_IG_cor_bad_ESS)  # strip off the run info
  all_high_IG_cor_bad_ESS<-unique(all_high_IG_cor_bad_ESS)  # get only the unique entries here

## sample out 1000 analyses each that have high IG correlation and either good or bad ESS
  if(length(all_high_IG_cor_bad_ESS)>1000){
    set.seed(666)
    high_IG_cor_bad_ESS_LARGE<-sample(all_high_IG_cor_bad_ESS, 1000)
  }else{
    high_IG_cor_bad_ESS_LARGE<-all_high_IG_cor_bad_ESS 
  }
  
  if(length(high_IG_cor_good_conv)>1000){
    set.seed(138)
    high_IG_cor_good_conv_LARGE<-sample(high_IG_cor_good_conv, 1000)
  }else{
    high_IG_cor_good_conv_LARGE<-high_IG_cor_good_conv
  }
  
  
  # Write a csv file containing the lines for each of these
  ## I've done this already, so commented out for now
    # write.csv(high_IG_cor_bad_ESS_LARGE, file="high_IG_cor_bad_ESS_LARGE.csv")
    # write.csv(high_IG_cor_good_conv_LARGE, file="high_IG_cor_good_conv_LARGE.csv")
  
  ## Here, we'll read in the csv's
  IGset_1<-read.csv("high_IG_cor_bad_ESS_LARGE.csv", stringsAsFactors = FALSE)[,2]
  IGset_2<-read.csv("high_IG_cor_good_conv_LARGE.csv", stringsAsFactors = FALSE)[,2]
  
  ## Write out an R data file of the diagnostics for these IG subset analyses
  datasets_IG_reestimated_LARGE<-c(high_IG_cor_bad_ESS_LARGE, high_IG_cor_good_conv_LARGE)
  orig_diags_IG_subsets_LARGE<-diagnostics[unique(datasets_IG_reestimated_LARGE)]  # the original diagnostics for the datasets that we reestimated trees for with different heating and with nst=mixed
  # get out the pass/fail table for these same sets of analyses
  to_pull_IG_LARGE<-which(gsub("\\.nex.*$|\\.run.*$", "", rownames(pass_fail_combined)) %in% unique(datasets_IG_reestimated_LARGE))  # get out the indices of the rows to pull out of the combined pass/fail table
  pass_fail_IG_orig_reestimated_LARGE<-pass_fail_combined[to_pull_IG_LARGE,]
  # Identify parameters that fail Geweke's in each chain
  failed_geweke_orig_IG_reest_pre_LARGE<-unlist(lapply(orig_diags_IG_subsets_LARGE, function(x) x$Params_fail_Geweke), recursive=FALSE)
  failed_geweke_orig_IG_reest_LARGE<-lapply(failed_geweke_orig_IG_reest_pre_LARGE, names)   # get the actual names of these parameters as the value in the object
  # Identify parameters that failed ESS in each chain
  failed_ESS_orig_IG_reest_pre_LARGE<-unlist(lapply(orig_diags_IG_subsets_LARGE, function(x) x$Params_ESS_less_200), recursive=FALSE)
  failed_ESS_orig_IG_reest_LARGE<-lapply(failed_ESS_orig_IG_reest_pre_LARGE, names)   # get the actual names of these parameters as the value in the object
  
  # save these as a list
  orig_IG_diags_LARGE<-list(pass_fail_IG_orig_reestimated_LARGE, failed_geweke_orig_IG_reest_LARGE, failed_ESS_orig_IG_reest_LARGE, orig_diags_IG_subsets_LARGE)
  names(orig_IG_diags_LARGE)<-c("pass_fail", "pars_fail_Geweke", "Pars_failed_ESS", "raw_diags")
  save(orig_IG_diags_LARGE, file="Orig_diags_for_IG_reestimated.RData")
  
  ## Get the burnin for each of these analyses 
  high_IG_cor_bad_ESS_LARGE_burns<-lapply(high_IG_cor_bad_ESS_LARGE, function(x) max(unlist(diagnostics[[x]]$Burnin))) # get the max burnin of the two chains in an analysis
  names(high_IG_cor_bad_ESS_LARGE_burns)<-high_IG_cor_bad_ESS_LARGE
  # same for analyses below
  high_IG_cor_good_conv_LARGE_burns<-lapply(high_IG_cor_good_conv_LARGE, function(x) max(unlist(diagnostics[[x]]$Burnin))) # get the max burnin of the two chains in an analysis
  names(high_IG_cor_good_conv_LARGE_burns)<-high_IG_cor_good_conv_LARGE
  
  IG_ex_chains_burn_LARGE<-list(high_IG_cor_bad_ESS_LARGE_burns, high_IG_cor_good_conv_LARGE_burns)
  names(IG_ex_chains_burn_LARGE)<-c("high_IG_cor_bad_ESS_LARGE_burns", "high_IG_cor_good_conv_LARGE_burns")
  save(IG_ex_chains_burn_LARGE, file="example_IG_chains_burn_LARGE.RData")

  
##########################################################################################  
##########################################################################################
##########################################################################################
## Manually tweak some analyses that perform poorly to try to improve convergence

## Read in a set of analyses that were pulled out manually previously
setwd(manual_dir)
manual_reanalyze<-read.csv("Manually_reanalyze.csv", stringsAsFactors=FALSE)[,2]
## Find how many parameters fail ESS for each anaysis
  failed_ESS<-lapply(diagnostics, function(x) unlist(x$Params_ESS_less_200))
  num_failed_ESS<-lapply(failed_ESS, length) 
  num_failed_ESS[manual_reanalyze]  # among the highest number of ESS failures across the datasets

  # Get the original diagnostics
  orig_manual_raw_diags<-diagnostics[manual_reanalyze]
  
  ## Write out the diags in the same way as done for the large sets of reestimates
  # Get pass/fail table
  to_pull_manual<-which(gsub("\\.nex.*$|\\.run.*$", "", rownames(pass_fail_combined)) %in% names(orig_manual_raw_diags))  # get out the indices of the rows to pull out of the combined pass/fail table
  pass_fail_orig_manual<-pass_fail_combined[to_pull_manual,]
  # Identify parameters that fail Geweke's in each chain
  failed_geweke_orig_manual_pre<-unlist(lapply(orig_manual_raw_diags, function(x) x$Params_fail_Geweke), recursive=FALSE)
  failed_geweke_orig_manual<-lapply(failed_geweke_orig_manual_pre, names)   # get the actual names of these parameters as the value in the object
  # Identify parameters that failed ESS in each chain
  failed_ESS_orig_manual_pre<-unlist(lapply(orig_manual_raw_diags, function(x) x$Params_ESS_less_200), recursive=FALSE)
  failed_ESS_orig_manual<-lapply(failed_ESS_orig_manual_pre, names)   # get the actual names of these parameters as the value in the object
  
  # save these as a list
  setwd(write_to)
  orig_manual_diags<-list(pass_fail_orig_manual, failed_geweke_orig_manual, failed_ESS_orig_manual, orig_manual_raw_diags)
  names(orig_manual_diags)<-c("pass_fail", "pars_fail_Geweke", "Pars_failed_ESS", "raw_diags")
  save(orig_manual_diags, file="orig_manual_diags.RData")

  ## Get the burnin for each of these analyses 
  manual_burns<-lapply(manual_reanalyze, function(x) max(unlist(diagnostics[[x]]$Burnin))) # get the max burnin of the two chains in an analysis
  names(manual_burns)<-manual_reanalyze

  save(manual_burns, file="orig_manual_burns.RData")
  
####################################################################################################
####################################################################################################
############  PCA of convergence diagnostics to see what's up
####################################################################################################
####################################################################################################  
  
  
  ### Add in autocorrelation times to object with ESS values, etc.
  act_table_PCA<-act_table
  colnames(act_table_PCA)<-paste(colnames(act_table), "ACT", sep=".")
  data_for_pca_2<-merge(ess_data_tl_branch, act_table_PCA, by.x="chain_name", by.y="row.names", all.x=TRUE) # I know that data_for_pca_2 starts at 2 and it's confusing...there was a 1 before that I deleted and I'm afraid of messing up if I renumber everything--apologies to anyone who's bothered to look at this and is tremendously offended
  
  
  ## Get out the PRSF values and merge with the object with other relevant info
  psrf_1<-lapply(diagnostics, function(x) x$PSRF$psrf[,1])
  psrf_names<-unique(unlist(lapply(psrf_1, names)))    ## Different analyses included different parameters (e.g., GTR vs HKY), need to find the full list of parameters that were used
  psrf<-do.call(rbind, lapply(psrf_1, function(x) x[match(psrf_names, names(x))]))
  colnames(psrf)<-paste(psrf_names, "psrf", sep=".")
  data_for_pca_3<-merge(data_for_pca_2, psrf, by.x="analysis_name", by.y="row.names", all.x=TRUE)  # merge together
  colnames(data_for_pca_3)[colnames(data_for_pca_3) %in% c("LnL", "TL", "r.A...C.", "r.A...G.", "r.A...T.", "r.C...G.", "r.C...T.", "r.G...T.", "pi.A.", "pi.C.", "pi.G.", "pi.T.", "alpha", "pinvar", "topo", "kappa", "LnPr")]<-paste(colnames(data_for_pca_3)[colnames(data_for_pca_3) %in% c("LnL", "TL", "r.A...C.", "r.A...G.", "r.A...T.", "r.C...G.", "r.C...T.", "r.G...T.", "pi.A.", "pi.C.", "pi.G.", "pi.T.", "alpha", "pinvar", "topo", "kappa", "LnPr")], "ESS", sep=".")
  colnames(data_for_pca_3)[colnames(data_for_pca_3)=="asdsf_for_TL"]<-"asdsf"

  ## Add in the geweke's z scores
  gewZ_table_pca<-gewZ_table
  rownames(gewZ_table_pca)<-gsub("nex", "",rownames(gewZ_table))
  colnames(gewZ_table_pca)<-paste(colnames(gewZ_table), "Z", sep="__")
  data_for_pca_4<-merge(data_for_pca_3, gewZ_table_pca, by.x="chain_name", by.y="row.names", all.x=TRUE)  # merge together
  
  ## Trim out the chain and analysis name stuff from this before feeding it into the PCA
  data_for_pca<-data_for_pca_4[,which(!names(data_for_pca_4) %in% c("chain_name", "analysis_name"))]

  # This standard PCA doesn't work with missing data
  ## Use pcaMethods library to use methods that allow for imputation of missing data
  resPCA<-pca(data_for_pca, nPcs=13, method="nipals", scale="uv")  
  # these settings are temporary for exploration
  pdf(file="test_PCA.pdf", width=25, height=25)
  slplot(resPCA)
  plotPcs(resPCA, pcs=1:4)
  dev.off()
  
  ## Try removing columns that have missing data
  params_for_reduced<-c("LnL.ESS", "TL.ESS", "topo.ESS", "ntax", "nchar", "asdsf", "tl_branch", "LnL.psrf", "TL.psrf")
  reduced_for_PCA<-data_for_pca[,params_for_reduced]
  resPCAreduced<-pca(reduced_for_PCA, nPcs=6, method="svd", scale="uv")  
  pdf(file="test_PCA_Complete.pdf", width=20, height=20)
  slplot(resPCAreduced)
  plotPcs(resPCAreduced, pcs=1:4)
  dev.off()
  
  
  ## Try removing a few rows so that we can include some acutocorrelation times
  data_for_pca_reduced_rows<-data_for_pca[-which(is.na(data_for_pca$LnL.ACT)),]
  params_for_reduced_rows<-c("LnL.ESS", "TL.ESS", "topo.ESS", "ntax", "nchar", "asdsf", "tl_branch", "LnL.psrf", "TL.psrf", "LnL.ACT","TL.ACT", "LnL__Z", "TL__Z")
  reduced_rows_for_PCA<-data_for_pca_reduced_rows[,params_for_reduced_rows]
  resPCAreduced_rows<-pca(reduced_rows_for_PCA, nPcs=6, method="svd", scale="uv")  
  pdf(file="test_PCA_complete_reducedRows.pdf", width=20, height=20)
  slplot(resPCAreduced_rows)
  plotPcs(resPCAreduced_rows, pcs=1:4)
  dev.off()
  
# save.image("Diagnostic_summary_all.RData")
# load("Diagnostic_summary_all.RData")
  
  
