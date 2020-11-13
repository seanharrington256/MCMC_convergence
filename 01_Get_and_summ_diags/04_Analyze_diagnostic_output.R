### Script to summarize diagnotic output generated from the 03_Calc_diags.R script

# This part assumes that the parent directory of the directory containing this script contains a folder named
#    Diag_output, which contains the output from all batches from script 03_Calc_diags.R
#    in that same directory, a directory for various figures will be created, with a subdirectory inside it for 
#    to save figures that were used in the manuscript
#  Note that there is some exploratory stuff still in here, not everything made it into the manuscript

# Set the working directory to one level up - we will then be in the base folder containing all of the subfolders we need
# only do this if the basedir object doesn't exist - prevents accidental rerunning of this step from another directory and messing up all paths
if(!exists("basedir")){
  setwd("../")
  basedir<-getwd()  # get the path to the base directory
}


## Set up some directories
storage_all<-paste0(basedir, "/Diag_output")  # location of the summaries of diagnostic output output by 03_Calc_diags.R
write_to<-paste0(basedir, "/Figs")  # directory to dump various figures into 
figs_ms<-paste0(basedir, "/Figs/figs_ms/")   # diretory to dump figures specifically used for the manuscript into

if(!dir.exists(write_to)){ # create the figs directory if it doesn't exist yet 
  dir.create(write_to)
}
if(!dir.exists(figs_ms)){ # same for the figs_ms directory
  dir.create(figs_ms)
}

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
library(ggbiplot)


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
# Also remove morf_448818_11codonalign.trim_interval.ctx, becuase it had no Nexus file 
diagnosticsAll<-diagnosticsAll[!names(diagnosticsAll) %in% "morf_448818_11codonalign.trim_interval.ctx"]


# For each analysis, get out a list of which chains passed and failed the diagnostic tests
pass_fail<-lapply(diagnosticsAll, function(x) as.data.frame(x[[1]]))
pass_fail_combined<-do.call(rbind.fill, pass_fail) #Combine these into a single matrix where each row is a chain
rownames(pass_fail_combined)<-unlist(lapply(diagnosticsAll, function(x) names(x$"ACT params")))

## Get the proportion of anlyses that pass each threshold for convergence
prop_pass<-apply(pass_fail_combined, 2, function(x) length(which(x==TRUE))/length(x))
prop_pass

# Make a csv of this
  prop_pass_table1<-rbind(gsub("Pass ", "", names(prop_pass)), prop_pass)
  prop_pass_table2<-gsub("ESS", "ESS > 200", prop_pass_table1)
  prop_pass_table3<-gsub("ASDSF", "ASDSF < 0.01", prop_pass_table2)
  prop_pass_table4<-gsub("SF corr", "SF corr > 0.9", prop_pass_table3)
  prop_pass_table<-gsub("PSRF", "PSRF < 1.02", prop_pass_table4)
  colnames(prop_pass_table)<-rep(NA, ncol(prop_pass_table))
  colnames(prop_pass_table)[[1]]<-"Proportion of chains passing each diagnostic"
  setwd(figs_ms)
  write.csv(prop_pass_table, file="table_pass_diags.csv")
  setwd(write_to)

### From here down, start mostly working only with analyses that have passed all diagnostic thresholds
  ## i.e., that according to our best evidence, have converged - don't include columns for passing with only 2 chains- if analyses passed with all 4 chains, we want to use them
  pass_all<-lapply(pass_fail, function(x) !FALSE %in% as.matrix(x)[,!colnames(as.matrix(x)) %in% c("Pass ASDSF 2 chains", "Pass SF corr 2 chains", "Pass PSRF 2 chains")]) 
  pass_all<-pass_all[unlist(pass_all)]
  
  # get out the diognostics for analyses that passed all diagnostic tests
  diagnostics<-diagnosticsAll[names(pass_all)]

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
      scale_x_continuous(limits = c(4800000, 6250000))+ # cut off the long right tail
      stat_density_ridges(quantile_lines = TRUE, quantiles = 2, scale=2.5)+  # add in a line for the median of each
      scale_y_discrete(expand = expansion(add = c(0.2, 2.8)))
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
    # Rank these parameters and see what has highest ESS
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
      scale_y_discrete(expand = expansion(add = c(0.2, 2.8)))+
      scale_x_continuous(limits = c(500, 800))
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
      scale_y_discrete(expand = expansion(add = c(0.2, 2.8)))+
      xlab("Autocorrelation time")
    dev.off()
    setwd(write_to)
    
    
#     
# # Compare topo pseudo ESS to approximate ESS - this was exploratory and is no longer in the manuscript, but it's retained here if anyone digs this far and wants to see it, and also, because 
#     approx_ess_pre<-unlist(lapply(diagnostics, function(x) x$"Approx topo ESS"), recursive=FALSE) # Get out just the approx ESS part of the diagnostics list
#     approx_ess_operator<-lapply(approx_ess_pre, function(x) x$operator)  # Get out the operator for each element
#     approx_ess_value<-lapply(approx_ess_pre, function(x) x$approx.ess)  # Get out the value for each element
#     to_compare<-which(approx_ess_operator=="=") # make an object for indexing to compare only those analyses for which approx ESS was calculated with "="
#     approx_for_compare<-approx_ess_value[to_compare]  # Pull out just the ESS values that have "="
#     pseudo_ESS_conf<-unlist(lapply(diagnostics, function(x) x$"Topo pseudo ess confidence"), recursive=FALSE) # Get out just the approx ESS part of the diagnostics list
#     pseudo_ESS_med_for_compare<-lapply(pseudo_ESS_conf, function(x) x$median.ess)[to_compare]
#     diffs_approx_pseudo<-unlist(approx_for_compare)-unlist(pseudo_ESS_med_for_compare)
#     mean(diffs_approx_pseudo)
#   #Compare the approximate ESS to the maximum pseudo ESS  
#     pseudo_ESS_upper_for_compare<-lapply(pseudo_ESS_conf, function(x) x$ci.upper)[to_compare]
#     diffs_approx_pseudo_upper<-unlist(approx_for_compare)-unlist(pseudo_ESS_upper_for_compare)
#     mean(diffs_approx_pseudo_upper)
# 
#     ## Plot out as overlapping histograms - again, this didn't end up getting used in the manuscript
#     setwd(figs_ms)
#     pdf(width=5, height=4.5, file="Fig_S1_Diffs_approx_pseudo_ESS.pdf")
#     hist(diffs_approx_pseudo_upper, main=NULL, xlab="Approximate minus psuedo topological ESS", ylab="Frequency", breaks=50, col=rgb(1,0,0,0.5))
#     hist(diffs_approx_pseudo, add=TRUE, breaks=50, col=rgb(0,0,1,0.5))
#     legend("topright", c("Upper 95% pseudo", "median Pseudo"), col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch=15, cex=0.7)
#     dev.off()
#     setwd(write_to)    
#     
    
    names(ess_all)<-unlist(lapply(diagnostics, function(x) names(x$"ACT params")))
    
    
##########################################################################################
##########################################################################################
##########################################################################################
### Compare convergence to dataset properties
##########################################################################################
##########################################################################################
##########################################################################################

## Read in the csv that states which study each analysis is from
study_names<-read.csv(paste0(basedir, "/Study_names.csv"))
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
    ess_data_chars<-merge(ess_data_chars_1, study_names, by.x="analysis_name", by.y="Dataset", all.x=TRUE) # finally finally, add in the study citation

    cor(as.numeric(ess_data_chars$nchar), as.numeric(ess_data_chars$ntax), use="complete.obs")  ## are number of taxa and characters correlated? -- not really
    
    
## look at the numbers of taxa and characters across analyses
    col.index<-as.factor(ess_data_chars$Study)
    setwd(figs_ms)
    pdf(file="Dataset_props.pdf", height=4, width=4)
    plot(ess_data_chars$ntax, ess_data_chars$nchar, cex=0.3, col=c("darkviolet", "darkgoldenrod3", "deepskyblue", "coral2")[col.index], pch=20, xlab="# taxa", ylab="# characters")
    legend(x="topright", legend = levels(col.index), col=c("darkviolet", "darkgoldenrod3", "deepskyblue", "coral2"), pch=20, pt.cex=1, cex=0.7)
    dev.off()
    setwd(write_to)

# Correlations and plots of ess for different parameters against Ntax and Nchars
    # Make a list of the different parameters that have ESS values
    par_ess<-colnames(ess_data_chars)[!colnames(ess_data_chars) %in% c("analysis_name", "ntax", "nchar", "Study")]
    # For each of these parameters, do multiple regressions with number of taxa and number of characters
    ess_cors_tax_chars<-sapply(par_ess, function(x) lm(as.numeric(ess_data_chars[,x]) ~ as.numeric(ess_data_chars[,"ntax"]) + as.numeric(ess_data_chars[,"nchar"])))
    ess_cors_tax_charsCoeff<-lapply(ess_cors_tax_chars, function(x) summary(x)$coefficients)  # get the coefficients and p values for each of these 

    ## Check out if these correlations vary by dataset
    # Amniotes - LnPr gets removed out during the correlation step because it wasn't sampled for any of these analyses
    amn_ess_data_chars<-ess_data_chars[which(ess_data_chars[,"Study"]=="Amniotes"),]
    amn_params_remove<-names(which(sapply(amn_ess_data_chars, function(x)all(is.na(x)))))  # Find any columns that are completely NA to also remove these from the list of params to correlate
    par_amn_ess<-colnames(amn_ess_data_chars)[!colnames(amn_ess_data_chars) %in% c("analysis_name", "ntax", "nchar", "Study", amn_params_remove)]
    amn_ess_cors_tax_chars<-sapply(par_amn_ess, function(x) lm(as.numeric(amn_ess_data_chars[,x]) ~ as.numeric(amn_ess_data_chars[,"ntax"]) + as.numeric(amn_ess_data_chars[,"nchar"])))
    amn_ess_cors_tax_charsCoeff<-lapply(amn_ess_cors_tax_chars, function(x) summary(x)$coefficients)  # get the coefficients and p values for each of these 

    # Barcoding - pinvar gets removed out during the correlation step because it wasn't sampled for any of these analyses
    bc_ess_data_chars<-ess_data_chars[which(ess_data_chars[,"Study"]=="Barcoding"),]
    bc_params_remove<-names(which(sapply(bc_ess_data_chars, function(x)all(is.na(x)))))  # Find any columns that are completely NA to also remove these from the list of params to correlate
    par_bc_ess<-colnames(bc_ess_data_chars)[!colnames(bc_ess_data_chars) %in% c("analysis_name", "ntax", "nchar", "Study", bc_params_remove)]
    bc_ess_cors_tax_chars<-sapply(par_bc_ess, function(x) lm(as.numeric(bc_ess_data_chars[,x]) ~ as.numeric(bc_ess_data_chars[,"ntax"]) + as.numeric(bc_ess_data_chars[,"nchar"])))
    bc_ess_cors_tax_charsCoeff<-lapply(bc_ess_cors_tax_chars, function(x) summary(x)$coefficients)  # get the coefficients and p values for each of these 
    
    # mtDNA
    mt_ess_data_chars<-ess_data_chars[which(ess_data_chars[,"Study"]=="mtDNA"),]
    mt_params_remove<-names(which(sapply(mt_ess_data_chars, function(x)all(is.na(x)))))  # Find any columns that are completely NA to also remove these from the list of params to correlate
    par_mt_ess<-colnames(mt_ess_data_chars)[!colnames(mt_ess_data_chars) %in% c("analysis_name", "ntax", "nchar", "Study", mt_params_remove, "kappa")]  ### Also pull kappa out of this because only 2 analyses have kappa
    mt_ess_cors_tax_chars<-sapply(par_mt_ess, function(x) lm(as.numeric(mt_ess_data_chars[,x]) ~ as.numeric(mt_ess_data_chars[,"ntax"]) + as.numeric(mt_ess_data_chars[,"nchar"])))
    mt_ess_cors_tax_charsCoeff<-lapply(mt_ess_cors_tax_chars, function(x) summary(x)$coefficients)  # get the coefficients and p values for each of these 
    
    # Phylota
    ph_ess_data_chars<-ess_data_chars[which(ess_data_chars[,"Study"]=="PhyLoTA"),]
    ph_params_remove<-names(which(sapply(ph_ess_data_chars, function(x)all(is.na(x)))))  # Find any columns that are completely NA to also remove these from the list of params to correlate
    par_ph_ess<-colnames(ph_ess_data_chars)[!colnames(ph_ess_data_chars) %in% c("analysis_name", "ntax", "nchar", "Study", ph_params_remove)]
    ph_ess_cors_tax_chars<-sapply(par_ph_ess, function(x) lm(as.numeric(ph_ess_data_chars[,x]) ~ as.numeric(ph_ess_data_chars[,"ntax"]) + as.numeric(ph_ess_data_chars[,"nchar"])))
    ph_ess_cors_tax_charsCoeff<-lapply(ph_ess_cors_tax_chars, function(x) summary(x)$coefficients)  # get the coefficients and p values for each of these 
    
  
    ## Make a big table with the correlations of each parameter's ESS with number of taxa and characters - include the correlation when considering all datasets together and when looking at each individually
  
    # First need to rename a bunch of rows
    for(i in 1:length(ess_cors_tax_charsCoeff)){
      par<-names(ess_cors_tax_charsCoeff)[[i]]
      rownames(ess_cors_tax_charsCoeff[[i]])<-c(paste0("All ", par, " Intercept"), paste0("All ", par, " N taxa"), paste0("All ", par, " N chars"))
    }
    ess_taxChars_cors<-do.call(rbind, ess_cors_tax_charsCoeff)
    
    for(i in 1:length(amn_ess_cors_tax_charsCoeff)){
      par<-names(amn_ess_cors_tax_charsCoeff)[[i]]
      rownames(amn_ess_cors_tax_charsCoeff[[i]])<-c(paste0("Amniotes ", par, " Intercept"), paste0("Amniotes ", par, " N taxa"), paste0("Amniotes ", par, " N chars"))
    }
    amn_ess_taxChars_cors<-do.call(rbind, amn_ess_cors_tax_charsCoeff)

    for(i in 1:length(bc_ess_cors_tax_charsCoeff)){
      par<-names(bc_ess_cors_tax_charsCoeff)[[i]]
      rownames(bc_ess_cors_tax_charsCoeff[[i]])<-c(paste0("Barcoding ", par, " Intercept"), paste0("Barcoding ", par, " N taxa"), paste0("Barcoding ", par, " N chars"))
    }
    bc_ess_taxChars_cors<-do.call(rbind, bc_ess_cors_tax_charsCoeff)

    for(i in 1:length(mt_ess_cors_tax_charsCoeff)){
      par<-names(mt_ess_cors_tax_charsCoeff)[[i]]
      rownames(mt_ess_cors_tax_charsCoeff[[i]])<-c(paste0("mtDNA ", par, " Intercept"), paste0("mtDNA ", par, " N taxa"), paste0("mtDNA ", par, " N chars"))
    }
    mt_ess_taxChars_cors<-do.call(rbind, mt_ess_cors_tax_charsCoeff)
    
    for(i in 1:length(ph_ess_cors_tax_charsCoeff)){
      par<-names(ph_ess_cors_tax_charsCoeff)[[i]]
      rownames(ph_ess_cors_tax_charsCoeff[[i]])<-c(paste0("PhyLoTA ", par, " Intercept"), paste0("PhyLoTA ", par, " N taxa"), paste0("PhyLoTA ", par, " N chars"))
    }
    ph_ess_taxChars_cors<-do.call(rbind, ph_ess_cors_tax_charsCoeff)
    
    ESS_cor_table<-rbind(ess_taxChars_cors, amn_ess_taxChars_cors, bc_ess_taxChars_cors, mt_ess_taxChars_cors, ph_ess_taxChars_cors)
    ESS_cor_table<-ESS_cor_table[-grep("Intercept", rownames(ESS_cor_table)),]   # Drop out the intercept rows
    write.csv(ESS_cor_table, file="tax_chars_ESS_reg.csv")
    
    
    # Plot these out to pdf
    pdf(width=4, height=7, file="ESS correlations.pdf")
    par(mfrow=c(2,1))
    ## Loop to get the correlations of the paramters with the num tax and num chars
    for(i in par_ess){
      plot(ess_data_chars[,"ntax"], ess_data_chars[,i], cex=0.01, xlab="# taxa", ylab=paste(i, "ESS", sep=" "), main=paste("b =", round(ess_cors_tax_charsCoeff[[i]][2,1], 3), "p =", round(ess_cors_tax_charsCoeff[[i]][2,4], 3),  sep=" ")) # for each parameter, plot it against the number of taxa
      plot(ess_data_chars[,"nchar"], ess_data_chars[,i], cex=0.01, xlab="# characters", ylab=paste(i, "ESS", sep=" "), main=paste("b =", round(ess_cors_tax_charsCoeff[[i]][3,1], 3), "p =", round(ess_cors_tax_charsCoeff[[i]][3,4], 3),  sep=" ")) # for each parameter, plot it against the number of characters
    }
    dev.off()
    

    # # Look at the correlations of ESS for LnL with all other params' ESS
    # # First make list of non-lnl ESS parameters
    # par_ess_not_lnl<-colnames(ess_data_chars)[!colnames(ess_data_chars) %in% c("analysis_name", "ntax", "nchar", "LnL", "Study")]
    # # Run a multiple regression of all non-LnL ESS values with LnL ESS
    # ess_cors_ess_lnl_FULL<-summary(lm(ess_data_chars[,"LnL"] ~  ess_data_chars[,"TL"] + ess_data_chars[,"r.A...C."] + ess_data_chars[,"r.A...G."] + ess_data_chars[,"r.A...T."] + 
    #      ess_data_chars[,"r.C...G."] + ess_data_chars[,"r.C...T."] + ess_data_chars[,"r.G...T."] + ess_data_chars[,"pi.A."] +
    #      ess_data_chars[,"pi.C."] + ess_data_chars[,"pi.G."] + ess_data_chars[,"pi.T."] + ess_data_chars[,"alpha"] +
    #      ess_data_chars[,"pinvar"] + ess_data_chars[,"topo"]  + ess_data_chars[,"LnPr"]))
    # 
    # ### !!!!!! Important note  !!!!
    # ### Removed Kappa from this, because it has no overlap with the GTR parameters 
    # ### caused this error: Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
    # ###       0 (non-NA) cases
    # par_ess_not_lnl<-par_ess_not_lnl[!par_ess_not_lnl %in% "kappa"]
    # 
    # # rename some rows
    # rownames(ess_cors_ess_lnl_FULL$coefficients)<-c("Intercept", par_ess_not_lnl)
    # ess_cors_ess_lnl_FULL<-ess_cors_ess_lnl_FULL$coefficients
    # 
    # # Plot these out to pdf
    # pdf(width=4, height=9.5, file="ESS LnL correlations.pdf")
    # par(mfrow=c(2,1))
    # ## Loop to get the correlations of the paramters with the num tax and num chars
    # for(i in par_ess_not_lnl){
    #   plot(ess_data_chars[,i], ess_data_chars[,"LnL"], cex=0.01, xlab=paste(i, "ESS", sep=" "), ylab="LnL ESS", main=paste("b =", round(ess_cors_ess_lnl_FULL[i,1], 3),"p =", round(ess_cors_ess_lnl_FULL[i,4], 3), sep=" ")) # for each parameter, plot it against the number of taxa
    # }
    # dev.off()
    # 
    # setwd(figs_ms)
    # write.csv(ess_cors_ess_lnl_FULL, file="ESS_cors_LnL_ESS.csv")
    # setwd(write_to)

## extract out some multi-chain diagnostics
  # Pull out ASDSF from analyses
    ### analyses with 4 chains have $asdsf (for 4 chains) $asdsf2chains (for only two chains)
    # get out 2 chain analyses
    num_chains<-lapply(diagnostics, function(x) nrow(x$Convergence_Summary)) # get out an object stating the number of chains
    
    asdsf_only_2_chains<-lapply(diagnostics[num_chains==2], function(x) x$ASDSF_and_SF_corr$ASDSF)  # get the asdsf from analyses that only had 2 chains
    asdsf_only_2_chains_4<-lapply(diagnostics[num_chains==4], function(x) x$ASDSF_and_SF_corr$ASDSF$asdsf2chains) # get the asdsf from analyses that only had 4 chains but calculated based only on 2 chains
    asdsf_only<-c(asdsf_only_2_chains, asdsf_only_2_chains_4)  # combine all of the 2 chain asdsf estimates
    
    asdsf_only_4_chains<-lapply(diagnostics[num_chains==4], function(x) x$ASDSF_and_SF_corr$ASDSF$asdsf) # get those with 4 chains
    
    
    # merge the asdsf with ESS and dataset properties
    asdsf_ess_chars_2_chains<-merge(ess_data_chars, as.data.frame(unlist(asdsf_only)), by.x="analysis_name", by.y="row.names")
    colnames(asdsf_ess_chars_2_chains)[colnames(asdsf_ess_chars_2_chains)=="unlist(asdsf_only)"]<-"asdsf_only"

    asdsf_ess_chars_4_chains<-merge(ess_data_chars, as.data.frame(unlist(asdsf_only_4_chains)), by.x="analysis_name", by.y="row.names")
    colnames(asdsf_ess_chars_4_chains)[colnames(asdsf_ess_chars_4_chains)=="unlist(asdsf_only_4_chains)"]<-"asdsf_only"
    
    
    # Multiple regression of asdsf with nchars and ntax
    reg_asdsf_2_chains<-summary(lm(as.numeric(asdsf_ess_chars_2_chains[,"asdsf_only"]) ~ as.numeric(asdsf_ess_chars_2_chains[,"ntax"]) + as.numeric(asdsf_ess_chars_2_chains[,"nchar"])))
    reg_asdsf_2_chains_coeff<-reg_asdsf_2_chains$coefficients
    
    ## Check out and see if these correlations vary by dataset
      # Amniotes - LnPr gets removed out during the correlation step because it wasn't sampled for any of these analyses
        amn_asdsf_data_chars<-asdsf_ess_chars_2_chains[which(asdsf_ess_chars_2_chains[,"Study"]=="Amniotes"),]
        amn_reg_asdsf<-summary(lm(as.numeric(amn_asdsf_data_chars[,"asdsf_only"]) ~ as.numeric(amn_asdsf_data_chars[,"ntax"]) + as.numeric(amn_asdsf_data_chars[,"nchar"])))
        amn_reg_asdsf_coeff<-amn_reg_asdsf$coefficients
        plot(amn_asdsf_data_chars[,"ntax"], amn_asdsf_data_chars[,"asdsf_only"], cex=0.2, xlab="N taxa", ylab="ASDSF")
        plot(amn_asdsf_data_chars[,"nchar"], amn_asdsf_data_chars[,"asdsf_only"], cex=0.2, xlab="N chars", ylab="ASDSF")

        # Correlations are low - low variation
      # Barcoding - pinvar gets removed out during the correlation step because it wasn't sampled for any of these analyses
        bc_asdsf_data_chars<-asdsf_ess_chars_2_chains[which(asdsf_ess_chars_2_chains[,"Study"]=="Barcoding"),]
        bc_reg_asdsf<-summary(lm(as.numeric(bc_asdsf_data_chars[,"asdsf_only"]) ~ as.numeric(bc_asdsf_data_chars[,"ntax"]) + as.numeric(bc_asdsf_data_chars[,"nchar"])))
        bc_reg_asdsf_coeff<-bc_reg_asdsf$coefficients
        plot(bc_asdsf_data_chars[,"ntax"], bc_asdsf_data_chars[,"asdsf_only"], cex=0.2, xlab="N taxa", ylab="ASDSF")
        plot(bc_asdsf_data_chars[,"nchar"], bc_asdsf_data_chars[,"asdsf_only"], cex=0.2, xlab="N chars", ylab="ASDSF")
      # mtDNA
        mt_asdsf_data_chars<-asdsf_ess_chars_2_chains[which(asdsf_ess_chars_2_chains[,"Study"]=="mtDNA"),]
        mt_reg_asdsf<-summary(lm(as.numeric(mt_asdsf_data_chars[,"asdsf_only"]) ~ as.numeric(mt_asdsf_data_chars[,"ntax"]) + as.numeric(mt_asdsf_data_chars[,"nchar"])))
        mt_reg_asdsf_coeff<-mt_reg_asdsf$coefficients
        plot(mt_asdsf_data_chars[,"ntax"], mt_asdsf_data_chars[,"asdsf_only"], cex=0.2, xlab="N taxa", ylab="ASDSF")
        plot(mt_asdsf_data_chars[,"nchar"], mt_asdsf_data_chars[,"asdsf_only"], cex=0.2, xlab="N chars", ylab="ASDSF")
      # Phylota - all 4 chain analyese are Phylota - do 2 chain and 4 chain here
        ph_asdsf_4_data_chars<-asdsf_ess_chars_4_chains[which(asdsf_ess_chars_4_chains[,"Study"]=="PhyLoTA"),]
        ph_reg_asdsf_4<-summary(lm(as.numeric(ph_asdsf_4_data_chars[,"asdsf_only"]) ~ as.numeric(ph_asdsf_4_data_chars[,"ntax"]) + as.numeric(ph_asdsf_4_data_chars[,"nchar"])))
        ph_reg_asdsf_4_coeff<-ph_reg_asdsf_4$coefficients
        plot(ph_asdsf_4_data_chars[,"ntax"], ph_asdsf_4_data_chars[,"asdsf_only"], cex=0.2, xlab="N taxa", ylab="ASDSF")
        plot(ph_asdsf_4_data_chars[,"nchar"], ph_asdsf_4_data_chars[,"asdsf_only"], cex=0.2, xlab="N chars", ylab="ASDSF")
        
        # Now for 2 chains
        ph_asdsf_2_data_chars<-asdsf_ess_chars_2_chains[which(asdsf_ess_chars_2_chains[,"Study"]=="PhyLoTA"),]
        ph_reg_asdsf_2<-summary(lm(as.numeric(ph_asdsf_2_data_chars[,"asdsf_only"]) ~ as.numeric(ph_asdsf_2_data_chars[,"ntax"]) + as.numeric(ph_asdsf_2_data_chars[,"nchar"])))
        ph_reg_asdsf_2_coeff<-ph_reg_asdsf_2$coefficients
        plot(ph_asdsf_2_data_chars[,"ntax"], ph_asdsf_2_data_chars[,"asdsf_only"], cex=0.2, xlab="N taxa", ylab="ASDSF")
        plot(ph_asdsf_2_data_chars[,"nchar"], ph_asdsf_2_data_chars[,"asdsf_only"], cex=0.2, xlab="N chars", ylab="ASDSF")
        
  
    ## Get an object of SF corrs from just 2 chains for all analyses
    sfcors_only_2_chains<-lapply(diagnostics[num_chains==2], function(x) x$ASDSF_and_SF_corr$sf_corrs)  # get the sfcorrs from analyses that only had 2 chains
    sfcors_only_2_chains_4<-lapply(diagnostics[num_chains==4], function(x) x$ASDSF_and_SF_corr$sf_corrs$cor_freq_2chains) # get the sfcorrs from analyses that only had 4 chains but calculated based only on 2 chains
    sfcors_only<-c(sfcors_only_2_chains, sfcors_only_2_chains_4)  # combine all of the 2 chain asdsf estimates
    
    ## get an object of SF corrs from all 4 chains for analyses with 4 chains
    sfcors_only_4_chains<-lapply(diagnostics[num_chains==4], function(x) x$ASDSF_and_SF_corr$sf_corrs$avg_cor_freq)
    
    ## Merge sf correlations with the ess and dataset characteristics
    sfcors_asdsf_chars_2_chains<-merge(asdsf_ess_chars_2_chains, as.data.frame(unlist(sfcors_only)), by.x="analysis_name", by.y="row.names")
    colnames(sfcors_asdsf_chars_2_chains)[colnames(sfcors_asdsf_chars_2_chains)=="unlist(sfcors_only)"]<-"sfcors_only"
    
    sfcors_asdsf_chars_4_chains<-merge(asdsf_ess_chars_4_chains, as.data.frame(unlist(sfcors_only_4_chains)), by.x="analysis_name", by.y="row.names")
    colnames(sfcors_asdsf_chars_4_chains)[colnames(sfcors_asdsf_chars_4_chains)=="unlist(sfcors_only_4_chains)"]<-"sfcors_only"

    # multiple regression between split frequency correlations and n taxa and n characters
    reg_sfcors<-summary(lm(as.numeric(sfcors_asdsf_chars_2_chains[,"sfcors_only"]) ~ as.numeric(sfcors_asdsf_chars_2_chains[,"ntax"]) + as.numeric(sfcors_asdsf_chars_2_chains[,"nchar"])))
    reg_sfcors_coeff<-reg_sfcors$coefficients

    # Plot out the relationships to pdf
    pdf(width=4, height=7, file="SFcors_chars-tax_cors.pdf")
    par(mfrow=c(2,1))
    plot(sfcors_asdsf_chars_2_chains[,"ntax"], sfcors_asdsf_chars_2_chains[,"sfcors_only"], cex=0.01, xlab="N taxa", ylab="SF cors", main=paste("b =", round(reg_sfcors_coeff[2,1], 3), "p =",  round(reg_sfcors_coeff[2,4], 3),sep=" "))
    plot(sfcors_asdsf_chars_2_chains[,"nchar"], sfcors_asdsf_chars_2_chains[,"sfcors_only"], cex=0.01, xlab="N chars", ylab="SF cors", main=paste("b =", round(reg_sfcors_coeff[3,1], 3), "p =",  round(reg_sfcors_coeff[3,4], 3),sep=" "))
    # Zoom in the y axis - there seem to be a small number of analyses with very low SF correlations - zoom in to get a better picture of what most are doing
    plot(sfcors_asdsf_chars_2_chains[,"ntax"], sfcors_asdsf_chars_2_chains[,"sfcors_only"], cex=0.01, xlab="N taxa", ylab="SF cors", ylim=c(0.85, 1),main=paste("Zoomed in - b =", round(reg_sfcors_coeff[2,1], 3), "p =",  round(reg_sfcors_coeff[2,4], 3),sep=" "))
    plot(sfcors_asdsf_chars_2_chains[,"nchar"], sfcors_asdsf_chars_2_chains[,"sfcors_only"], cex=0.01, xlab="N chars", ylab="SF cors", ylim=c(0.85, 1), main=paste("Zoomed in - b =", round(reg_sfcors_coeff[3,1], 3), "p =",  round(reg_sfcors_coeff[3,4], 3),sep=" "))
    
    dev.off()
    
    ## Check out correlation between asdsf and sf cor - basically nonexistent: split frequency correlation is pretty much universally high
    cor_asdsf_sdcor<-cor(sfcors_asdsf_chars_2_chains[,"asdsf_only"], sfcors_asdsf_chars_2_chains[,"sfcors_only"], use="complete.obs")
    plot(sfcors_asdsf_chars_2_chains[,"asdsf_only"], sfcors_asdsf_chars_2_chains[,"sfcors_only"], cex=0.01, xlab="ASDSF", ylab="SF cors", ylim=c(0.85, 1), xlim=c(0, 0.007), main=paste("Zoomed in - corr =", round(cor_asdsf_sdcor, 3), sep=" "))
  
            
            
## Compare failure of some ASDSF, etc. 
## Find out which analyses have an asdsf >0.01, low split frequency correlation, or low topological ess values 
    ## Need to get out necessary diagnostics from full set, including failed analyses
            comp_fail1<-lapply(diagnosticsAll, function(x) x$ESS_values) # get the ess values out of the main diagnostics object
            comp_fail2<-lapply(comp_fail1, function(x) do.call(rbind, x)) # within each element of the list, rbind the 2 chains together
            # Use a for loop to assign the name of each analysis as a column for each element of the list
            comp_fail3<-list()
            for(i in 1:length(comp_fail2)){
              comp_fail3[[i]]<-cbind.data.frame(names(comp_fail2)[[i]], comp_fail2[[i]]) # cbind the name of the analysis into the matrix as the first column
              colnames(comp_fail3[[i]])[[1]]<-"analysis_name" # name the first column
            }
            comp_fail4<-do.call(rbind.fill, comp_fail3) # rbind this all together filling in empty columns with NA using rbind.fill
            # Then, combine this with dataset characteristics
            data_charsAll<-do.call(rbind, lapply(diagnosticsAll, function(x) x$analysis_pars[1:2]))
            comp_fail5<-merge(comp_fail4, data_charsAll, by.x="analysis_name", by.y="row.names", all.x=TRUE) # finally, merge together the ESS and dataset data
            comp_fail6<-merge(comp_fail5, study_names, by.x="analysis_name", by.y="Dataset", all.x=TRUE) # finally finally, add in the study citation
            
            ####### do comparisons completely separately for 2 & 4 chain analyses
              #####    don't want to compare just 2 chains of ASDSF because it's compared to passing ESS
              #####    for all 4 chains, etc. 
            num_chains_all<-lapply(diagnosticsAll, function(x) nrow(x$Convergence_Summary)) # get out an object stating the number of chains

            # split out the general data by 2 or 4 chain analyses
            comp_fail6_2_chains<-comp_fail6[which(comp_fail6$analysis_name %in% names(num_chains_all)[num_chains_all==2]),]
            comp_fail6_4_chains<-comp_fail6[which(comp_fail6$analysis_name %in% names(num_chains_all)[num_chains_all==4]),]
            
            ## Get asdsf based from 2 chain analyses
            asdsf_all_2_chains<-lapply(diagnosticsAll[num_chains_all==2], function(x) x$ASDSF_and_SF_corr$ASDSF)
            ## Get asdsf based on 4 chains for analyses that have it
            asdsf_all_4_chains<-lapply(diagnosticsAll[num_chains_all==4], function(x) x$ASDSF_and_SF_corr$ASDSF$asdsf)
              # also get asdsf from 2 chains for the 4 chain analyses
            asdsf_all_2from4<-lapply(diagnosticsAll[num_chains_all==4], function(x) x$ASDSF_and_SF_corr$ASDSF$asdsf2chains)
            
            # merge these in with the other dataset characteristics
            comp_fail7_2_chains<-merge(comp_fail6_2_chains, as.data.frame(unlist(asdsf_all_2_chains)), by.x="analysis_name", by.y="row.names")
              colnames(comp_fail7_2_chains)[colnames(comp_fail7_2_chains)=="unlist(asdsf_all_2_chains)"]<-"asdsf_All"
            comp_fail7_4_chainsA<-merge(comp_fail6_4_chains, as.data.frame(unlist(asdsf_all_4_chains)), by.x="analysis_name", by.y="row.names")
            comp_fail7_4_chains<-merge(comp_fail7_4_chainsA, as.data.frame(unlist(asdsf_all_2from4)), by.x="analysis_name", by.y="row.names")
              colnames(comp_fail7_4_chains)[colnames(comp_fail7_4_chains)=="unlist(asdsf_all_4_chains)"]<-"asdsf_All"
              colnames(comp_fail7_4_chains)[colnames(comp_fail7_4_chains)=="unlist(asdsf_all_2from4)"]<-"asdsf_All_2from4"
            
            
            ## Get sfcors based from 2 & 4 chain analyses
            sfcors_all_2_chains<-lapply(diagnosticsAll[num_chains_all==2], function(x) x$ASDSF_and_SF_corr$sf_corrs)  # get the sfcorrs from analyses that only had 2 chains
            sfcors_all_4_chains<-lapply(diagnosticsAll[num_chains_all==4], function(x) x$ASDSF_and_SF_corr$sf_corrs$avg_cor_freq)
            sfcors_all_2from4<-lapply(diagnosticsAll[num_chains_all==4], function(x) x$ASDSF_and_SF_corr$sf_corrs$cor_freq_2chains)
            
          
            # merge the sfcors with the dataset properties and asdsf
            comp_fail_2_chains<-merge(comp_fail7_2_chains, as.data.frame(unlist(sfcors_all_2_chains)), by.x="analysis_name", by.y="row.names")
              colnames(comp_fail_2_chains)[colnames(comp_fail_2_chains)=="unlist(sfcors_all_2_chains)"]<-"sfcors_All"
            
            
            comp_fail_4_chainsA<-merge(comp_fail7_4_chains, as.data.frame(unlist(sfcors_all_4_chains)), by.x="analysis_name", by.y="row.names")
            comp_fail_4_chains<-merge(comp_fail_4_chainsA, as.data.frame(unlist(sfcors_all_2from4)), by.x="analysis_name", by.y="row.names")
              colnames(comp_fail_4_chains)[colnames(comp_fail_4_chains)=="unlist(sfcors_all_4_chains)"]<-"sfcors_All"
              colnames(comp_fail_4_chains)[colnames(comp_fail_4_chains)=="unlist(sfcors_all_2from4)"]<-"sfcors_All_2from4"
              
            head(comp_fail_2_chains)
            head(comp_fail_4_chains)
            

    # Which fail topo diags???
    high_asdsf_2_chains<-comp_fail_2_chains[which(comp_fail_2_chains$asdsf_All>0.01),]
    high_asdsf_4_chains<-comp_fail_4_chains[which(comp_fail_4_chains$asdsf_All>0.01),]
    
    low_sf_cors_2_chains<-comp_fail_2_chains[which(comp_fail_2_chains$sfcors_All<0.9),] # really not sure what's an appropriate cutoff, 0.9 is probably too low, many more fail at higher cutoffs
    low_sf_cors_4_chains<-comp_fail_4_chains[which(comp_fail_4_chains$sfcors_All<0.9),]
    
    low_topo_ess_2_chains<-comp_fail_2_chains[which(comp_fail_2_chains$topo<200),]
    low_topo_ess_4_chains<-comp_fail_4_chains[which(comp_fail_4_chains$topo<200),]
    

  # It seems that many analyses are failing to get a topological ESS value above 200, but still passing most or all other diagnostic checks
    ## What are the dataset properties of the analyses that fail topo ESS compared to those that don't:
    high_topo_ess_2_chains<-comp_fail_2_chains[which(comp_fail_2_chains$topo>200),] # Make a second object that has all of the analyses that have topo ESS above 200
    high_topo_ess_4_chains<-comp_fail_4_chains[which(comp_fail_4_chains$topo>200),]
    
    ## Medians are probably much more useful than means
    median(unlist(high_topo_ess_2_chains$ntax))
    median(unlist(low_topo_ess_2_chains$ntax))
    median(unlist(comp_fail_2_chains$ntax)) # also look at the overall median N taxa
    # N taxa median for chains that have topo ESS >200 seems to basically match overall median of all analyses, N taxa for analyses with topo ESS <200 is much higher
    
    median(unlist(high_topo_ess_2_chains$nchar))
    median(unlist(low_topo_ess_2_chains$nchar))
    median(unlist(comp_fail_2_chains$nchar)) # also look at the overall median N chars
    # Looks like the median N chars for all analyses and analyses that pass are similar, n chars for failing is lower
    
    # Plot out histograms of N taxa for when topo ESS is <200 and >200
    pdf(width=6, height=9, file="Topo_ESS_failed_dataset properties.pdf")
    par(mfrow=c(2,1))
    hist(unlist(low_topo_ess_2_chains$ntax), col="red", xlab="# taxa", main="Topo ESS < 200", xlim=c(0,300), ylab=NULL)
    hist(unlist(high_topo_ess_2_chains$ntax), col="blue", xlab="# taxa", main="Topo ESS > 200", xlim=c(0,300), ylab=NULL)
    
    hist(unlist(low_topo_ess_2_chains$nchar), col="red", xlab="# chars", main="Topo ESS < 200", xlim=c(0,5000), ylab=NULL)
    hist(unlist(high_topo_ess_2_chains$nchar), col="blue", xlab="# chars", main="Topo ESS > 200", xlim=c(0,5000), ylab=NULL)
    dev.off()
    
        ## Same thing for 4 chains
              ## Medians are probably much more useful than means
              median(unlist(high_topo_ess_4_chains$ntax))
              median(unlist(low_topo_ess_4_chains$ntax))
              median(unlist(comp_fail_4_chains$ntax)) # also look at the overall median N taxa
              # N taxa median for chains that have topo ESS >200 seems to basically match overall median of all analyses, N taxa for analyses with topo ESS <200 is much higher
              
              median(unlist(high_topo_ess_4_chains$nchar))
              median(unlist(low_topo_ess_4_chains$nchar))
              median(unlist(comp_fail_4_chains$nchar)) # also look at the overall median N chars
              # Looks like the median N chars for all analyses and analyses that pass are similar, n chars for failing is lower
              
              # Plot out histograms of N taxa for when topo ESS is <200 and >200
              pdf(width=6, height=9, file="Topo_ESS_failed_dataset properties_4chains.pdf")
              par(mfrow=c(2,1))
              hist(unlist(low_topo_ess_4_chains$ntax), col="red", xlab="# taxa", main="Topo ESS < 200", xlim=c(0,300), ylab=NULL)
              hist(unlist(high_topo_ess_4_chains$ntax), col="blue", xlab="# taxa", main="Topo ESS > 200", xlim=c(0,300), ylab=NULL)
              
              hist(unlist(low_topo_ess_4_chains$nchar), col="red", xlab="# chars", main="Topo ESS < 200", xlim=c(0,5000), ylab=NULL)
              hist(unlist(high_topo_ess_4_chains$nchar), col="blue", xlab="# chars", main="Topo ESS > 200", xlim=c(0,5000), ylab=NULL)
              dev.off()
    
    

    
    # How many of the chains that fail topo ESS pass LnL ESS?
    n_low_topoESS_2_chains<-nrow(low_topo_ess_2_chains) # total number of analyses that have low topo ESS
    n_low_topo_and_LnL_ESS_2_chains<-length(low_topo_ess_2_chains[which(low_topo_ess_2_chains$LnL<200),][,1]) # fail topo ESS, fail LnL ESS
    n_low_LnL_ESS_2_chains<-nrow(comp_fail_2_chains[which(comp_fail_2_chains$LnL<200),])
    ## For 2 chains, all 174 chains that have low topo ESS have LnL ESS>200 and the 23 chains with LnL<200 all have topo>200
    
    # Make a Venn diagram to show this
    # setwd(figs_ms)
    pdf(file="Fig_5A_LnL_topo_ESS_venn_2_chains.pdf", height=5, width=5)
    draw.pairwise.venn(n_low_LnL_ESS_2_chains, n_low_topoESS_2_chains, n_low_topo_and_LnL_ESS_2_chains, category = c("LnL ESS < 200", "Topological ESS < 200"), lty = rep("blank", 
          2), fill = c("aquamarine2", "darkviolet"), alpha = rep(0.5, 2), cat.pos = c(0, 
              0), cat.dist = rep(0.025, 2), scaled = TRUE, cex=1.5, cat.cex=1.5)                     
    dev.off()
    # setwd(write_to)

    
    ############# Same for 4 chains with topo ESS and LnL ESS
          # How many of the chains that fail topo ESS pass LnL ESS?
          n_low_topoESS_4_chains<-nrow(low_topo_ess_4_chains) # total number of analyses that have low topo ESS
          n_low_topo_and_LnL_ESS_4_chains<-length(low_topo_ess_4_chains[which(low_topo_ess_4_chains$LnL<200),][,1]) # fail topo ESS, fail LnL ESS
          n_low_LnL_ESS_4_chains<-nrow(comp_fail_4_chains[which(comp_fail_4_chains$LnL<200),])
          ## For 2 chains, all 174 chains that have low topo ESS have LnL ESS>200 and the 23 chains with LnL<200 all have topo>200
          
          # Make a Venn diagram to show this
          setwd(figs_ms)
          pdf(file="Fig_5A_LnL_topo_ESS_venn_4_chains.pdf", height=5, width=5)
          draw.pairwise.venn(n_low_LnL_ESS_4_chains, n_low_topoESS_4_chains, n_low_topo_and_LnL_ESS_4_chains, category = c("LnL ESS < 200", "Topological ESS < 200"), lty = rep("blank", 
              2), fill = c("aquamarine2", "darkviolet"), alpha = rep(0.5, 2), cat.pos = c(0, 
                  0), cat.dist = rep(0.025, 2), scaled = TRUE, cex=1.5, cat.cex=1.5)                     
          dev.off()
          setwd(write_to)
    

    # How many of the chains that fail topo ess pass asdsf - start with 2 chains
    n_low_topoESS_2_chains
    n_low_topo_bad_asdsf_2_chains<-nrow(low_topo_ess_2_chains[which(low_topo_ess_2_chains$asdsf_All>0.01),]) # fail topo ESS, fail asdsf
    n_bad_asdsf_2_chains<-nrow(comp_fail_2_chains[which(comp_fail_2_chains$asdsf_All>0.01),])
    ## Make the Venn diagram for this
    # setwd(figs_ms)
    pdf(file="Fig_5B_Topo_ESS_ASDSF_venn_2_CHAINS.pdf", height=5, width=5)
    draw.pairwise.venn(n_low_topoESS_2_chains, n_bad_asdsf_2_chains, n_low_topo_bad_asdsf_2_chains, category = c("Topological ESS < 200", "ASDSF > 0.01"), lty = rep("blank",
          2), fill = c("darkviolet", "aquamarine2"), alpha = rep(0.5, 2), cat.pos = c(0,
             0), cat.dist = rep(0.025, 2), scaled = TRUE, cex=1.5, cat.cex=1.5)
    dev.off()
    # setwd(write_to)

    
            ## Do same thing for 4 chains
                n_low_topoESS_4_chains
                n_low_topo_bad_asdsf_4_chains<-nrow(low_topo_ess_4_chains[which(low_topo_ess_4_chains$asdsf_All>0.01),]) # fail topo ESS, fail asdsf
                n_bad_asdsf_4_chains<-nrow(comp_fail_4_chains[which(comp_fail_4_chains$asdsf_All>0.01),])
                ## Make the Venn diagram for this
                setwd(figs_ms)
                pdf(file="Fig_5B_Topo_ESS_ASDSF_venn_4_CHAINS.pdf", height=5, width=5)
                draw.pairwise.venn(n_low_topoESS_4_chains, n_bad_asdsf_4_chains, n_low_topo_bad_asdsf_4_chains, category = c("Topological ESS < 200", "ASDSF > 0.01"), lty = rep("blank",
                  2), fill = c("darkviolet", "aquamarine2"), alpha = rep(0.5, 2), cat.pos = c(0,
                    0), cat.dist = rep(0.025, 2), scaled = TRUE, cex=1.5, cat.cex=1.5)
                dev.off()
                setwd(write_to)
                
                
    # compare the split frequency correlations for the analyses that pass vs. fail topo ESS
    median(low_topo_ess_2_chains$sfcors_All)
    median(high_topo_ess_2_chains$sfcors_All)
    median(sfcors_asdsf_chars_2_chains$sfcors_only) # also look at the overall median N chars
    # These histograms below are both zoomed in because of low frequencies of really low sf correlations
    hist(low_topo_ess_2_chains$sfcors_All, col="red", xlab="SF corrs", main="Topo ESS < 200", breaks=200, xlim=c(0.95, 1))
    hist(high_topo_ess_2_chains$sfcors_All, col="blue", xlab="SF corrs", main="Topo ESS > 200", breaks=200, xlim=c(0.95, 1))

    
    
    ## Better way of looking at the analyses that fail topo ESS
    
    #### Important to note that "Pass ASDSF" in pass_fail_combined is based on 4 chains for analyses that have 4
      ## doing these next few parts so that the Pass XXX 2 chains percents are based only on entries that aren't NA 
      ## (i.e., they'll only apply to the 4 chain analyses)
    
    pass_fail_fail_asdsf<-pass_fail_combined[which(pass_fail_combined[,"Pass ASDSF"] == FALSE),]
    # Get the proportions of these analyses that fail other diagnostics
    prop_pass_asdsf_fail<-apply(pass_fail_fail_asdsf, 2, function(x) length(which(x==TRUE))/length(which(!is.na(x))))
    prop_pass_asdsf_fail # Elevated chances of failing other diagnostics, but pass a lot of them

    
    ## Write this out as a csv table for the ms 
    prop_pass_asdsf_fail_table1<-rbind(gsub("Pass ", "", names(prop_pass_asdsf_fail)), prop_pass_asdsf_fail)
    prop_pass_asdsf_fail_table2<-gsub("ESS", "ESS > 200", prop_pass_asdsf_fail_table1)
    prop_pass_asdsf_fail_table3<-gsub("ASDSF", "ASDSF < 0.01", prop_pass_asdsf_fail_table2)
    prop_pass_asdsf_fail_table4<-gsub("SF corr", "SF corr > 0.9", prop_pass_asdsf_fail_table3)
    prop_pass_asdsf_fail_table<-gsub("PSRF", "PSRF < 1.02", prop_pass_asdsf_fail_table4)

    colnames(prop_pass_asdsf_fail_table)<-rep(NA, ncol(prop_pass_asdsf_fail_table))
    colnames(prop_pass_asdsf_fail_table)[[1]]<-"Proportion of chains passing each diagnostic for analyses that have ASDSF > 0.01"
    setwd(figs_ms)
    write.csv(prop_pass_asdsf_fail_table, file="Prop_Pass_FAILED_asdsf.csv")
    setwd(write_to)    
  
    
    
    ## Look at the analyses that fail topo ESS
    pass_fail_fail_topo<-pass_fail_combined[which(pass_fail_combined[,"Pass Approx topo ESS"] == FALSE),]
    # Get the proportions of these analyses that fail other diagnostics
    prop_pass_topoESS_fail<-apply(pass_fail_fail_topo, 2, function(x) length(which(x==TRUE))/length(which(!is.na(x))))
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
    
## Look at chains that fail either ASDSF or PSRF and see how they do for single-chain diags
    fail_ASDSForPSRF<-apply(pass_fail_combined, 1, function(x) x["Pass ASDSF"]==FALSE | x["Pass PSRF"]==FALSE) # see which analyses failed either ASDSF or PSRF
    pass_fail_fail_multi<-pass_fail_combined[fail_ASDSForPSRF,]  # pull out the pass/fail table for just the analyses that failed ASDSF or PSRF
    prop_pass_fail_multi<-apply(pass_fail_fail_multi, 2, function(x) length(which(x==TRUE))/length(which(!is.na(x)))) # get out the proportion of analyses that pass other diags that fail either of the multichain diags
    prop_pass_fail_multi
    
    ## Make a Venn diagram of how many analyses fail ASDSF or PSRF
    n_fail_asdsf<-length(which(pass_fail_combined[,"Pass ASDSF"]==FALSE))
    n_fail_psrf<-length(which(pass_fail_combined[,"Pass PSRF"]==FALSE))
    n_fail_psrf_and_asdsf<-length(which(apply(pass_fail_combined, 1, function(x) x["Pass ASDSF"]==FALSE & x["Pass PSRF"]==FALSE)==TRUE))
    setwd(figs_ms)
    pdf(file="Fig_5C_ASDSF_PSRF_venn.pdf", height=5, width=5)
    draw.pairwise.venn(n_fail_asdsf, n_fail_psrf, n_fail_psrf_and_asdsf, category = c("ASDSF > 0.01", "PSRF > 1.02"), lty = rep("blank", 
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
    draw.pairwise.venn(n_fail_asdsfORpsrf, n_fail_any_ESS, n_fail_adsfORpsrfANDanyESS, category = c("ASDSF > 0.01 OR PSRF > 1.02", "Any non-topo ESS < 200"), lty = rep("blank", 
        2), fill = c("darkviolet", "aquamarine2"), alpha = rep(0.5, 2), cat.pos = c(0, 
          0), cat.dist = rep(0.025, 2), scaled = TRUE, cex=1.5, cat.cex=1.5)                     
    dev.off()
    setwd(write_to)
    

## Get the acceptance rates for each move:
    
acceptance_only_2_chains<-lapply(diagnostics[num_chains==2], function(x) x$acc_rates)  # get the acceptance from analyses that only had 2 chains
acceptance_only_2_chains<-acceptance_only_2_chains[which(!lapply(acceptance_only_2_chains, function(x) x[[1]])=="NA")] # remove any analyses that didn't have acceptance info
acceptance_4_chains<-lapply(diagnostics[num_chains==4], function(x) x$acc_rates) # get the acceptance from analyses that had 4 chains
acceptance_2_chains<-c(acceptance_only_2_chains, lapply(acceptance_4_chains, function(x) x[1:2,]))  # combine all of the 2 chain info - getting the first 2 runs out of the 4 chain analyses just requires taking the first 2 rows out of acceptance rate object
    

## Different analyses included different moves, need to find the full list of moves that were used
move_names_2_chains<-unique(unlist(lapply(acceptance_2_chains, FUN=function(x) colnames(x))))
move_names_4_chains<-unique(unlist(lapply(acceptance_4_chains, FUN=function(x) colnames(x))))

# Use this to order all parameters the same and make a table of gens to ess 200 for all for each chain
acceptance_table_2_chains<-do.call(rbind, lapply(acceptance_2_chains, function(x) x[,match(move_names_2_chains, colnames(x))]))
colnames(acceptance_table_2_chains)<-move_names_2_chains
acceptance_table_4_chains<-do.call(rbind, lapply(acceptance_4_chains, function(x) x[,match(move_names_4_chains, colnames(x))]))
colnames(acceptance_table_4_chains)<-move_names_4_chains


# Get mean, median, and sd acceptance for each move - 2 CHAINS
    mean_acc_2_chains<-apply(acceptance_table_2_chains, 2, mean, na.rm=TRUE)
    median_acc_2_chains<-apply(acceptance_table_2_chains, 2, median, na.rm=TRUE)
    sd_gens_acc_2_chains<-apply(acceptance_table_2_chains, 2, sd, na.rm=TRUE)
    range_acc_2_chains<-apply(acceptance_table_2_chains, 2, range, na.rm=TRUE)
    # Rank these parameters and see what converges fastest
    sorted_means_acc_2_chains<-sort(mean_acc_2_chains)
    sorted_median_acc_2_chains<-sort(median_acc_2_chains)
    sorted_sd_acc_2_chains<-sort(sd_gens_acc_2_chains)
    sorted_means_acc_2_chains
    sorted_median_acc_2_chains
    sorted_sd_acc_2_chains
    range_acc_2_chains

# Get mean, median, and sd acceptance for each move - 4 CHAINS
    mean_acc_4_chains<-apply(acceptance_table_4_chains, 2, mean, na.rm=TRUE)
    median_acc_4_chains<-apply(acceptance_table_4_chains, 2, median, na.rm=TRUE)
    sd_gens_acc_4_chains<-apply(acceptance_table_4_chains, 2, sd, na.rm=TRUE)
    range_acc_4_chains<-apply(acceptance_table_4_chains, 2, range, na.rm=TRUE)
    # Rank these parameters and see what converges fastest
    sorted_means_acc_4_chains<-sort(mean_acc_4_chains)
    sorted_median_acc_4_chains<-sort(median_acc_4_chains)
    sorted_sd_acc_4_chains<-sort(sd_gens_acc_4_chains)
    sorted_means_acc_4_chains
    sorted_median_acc_4_chains
    sorted_sd_acc_4_chains
    range_acc_4_chains
    
    
        
# Histograms of rates for topology moves
hist(acceptance_table_4_chains[,"ExtSPR(Tau,V)"])
hist(acceptance_table_4_chains[,"ExtTBR(Tau,V)"])
hist(acceptance_table_4_chains[,"NNI(Tau,V)"])
hist(acceptance_table_4_chains[,"ParsSPR(Tau,V)"])


## make a ridgeline plot of the acceptance rates for parameters
#  Note this is for 4 chains only for the sake of uniformity
# First transform the dataframe into a format easily readable by ggplot
accept_melted<-reshape2::melt(acceptance_table_4_chains, id.vars = NULL)
# Not every parameter is in every analysis, so need to prune out the NA's
accept_melted<-accept_melted[which(is.na(accept_melted$value)==FALSE),]
setwd(figs_ms)
pdf(file="Fig_4_Move_acceptance_ridges_4Chains.pdf", width=6, height=4)
ggplot(accept_melted, aes(y=Var2, x=value, fill = Var2)) +
  labs(title=NULL,y="Move", x = "Acceptance %")+
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, scale=4)+ 
  scale_fill_cyclical(values = c("aquamarine2", "darkviolet"))+
  scale_y_discrete(expand = expansion(add = c(0.2, 3.3)))
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
  ess_chars_by_chain_pre$analysis_name<-as.character(ess_chars_by_chain_pre$analysis_name)
  # this is somewhat tricky because different analyses have different numbers of chains
  for(i in unique(ess_data_chars$analysis_name)){ # for each unique analysis name:
    index<-which(ess_data_chars$analysis_name==i) # get an index of the analyses that match the i'th name
    ess_chars_by_chain_pre$analysis_name[index]<-paste(ess_data_chars$analysis_name[index], 1:length(index), sep=".") # for each chain in each analysis
  }
  colnames(ess_chars_by_chain_pre)[which(colnames(ess_chars_by_chain_pre)=="analysis_name")]<-"chain_name"
  ess_chars_by_chain<-cbind(ess_data_chars[,"analysis_name"], ess_chars_by_chain_pre)
  colnames(ess_chars_by_chain)[[1]]<-"analysis_name"
  
  ## Set up TL to be formatted for the correlations and merge it with ESS, etc. and then merge them together
  TL<-lapply(diag_TL, function(x) x$TL[which(names(x$TL)=="median")]) # get out just the median TL
  TL_cor_form<-cbind.data.frame(gsub(".median", "", names(unlist(TL))),  unlist(TL), stringsAsFactors=FALSE)  # make this into a dataframe with the analysis name as a column
  colnames(TL_cor_form)<-c("analysis_name", "med_TL") # add in column names
  ## Add in chain names
  TL_names_pre<-gsub(".median", "", names(unlist(TL))) # get the names of each analysis from TL
  for(i in unique(TL_names_pre)){ # for the unique TL names:
    index<-which(TL_names_pre==i) # get an index of chains matching that name
    TL_cor_form$analysis_name[index]<-paste(TL_names_pre[index], 1:length(index), sep=".") # append a chain index to each chain within an analysis
  }
  
  ess_data_TL<-merge(ess_chars_by_chain, TL_cor_form, by.x="chain_name", by.y="analysis_name", all.x=TRUE)
  # Add asdsf into this - use the ASDSF from 4 chains when there are 4
  asdsf_2and4_forTL<-c(asdsf_only_2_chains, asdsf_only_4_chains)
  
  # merge the asdsf with ESS & other data characteristics
  ESS_asdsf_TL<-merge(ess_data_TL, as.data.frame(unlist(asdsf_2and4_forTL)), by.x="analysis_name", by.y="row.names", all.x=TRUE)
  colnames(ESS_asdsf_TL)<-gsub("unlist\\(asdsf_2and4_forTL\\)", "asdsf", colnames(ESS_asdsf_TL))
  
  # Make a list of the different parameters that have ESS values include asdsf here
  par_ess_TL<-colnames(ESS_asdsf_TL)[!colnames(ESS_asdsf_TL) %in% c("analysis_name", "med_TL", "chain_name", "ntax", "nchar", "Study")]
  # Get the correlations of ESS values for each of these parameters and of asdsf with TL
  ### This is no longer used, just exploratory
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
  par_ess_cred<-colnames(ess_data_cred)[!colnames(ess_data_cred) %in% c("chain_name", "analysis_name", "cred_size", "ntax", "nchar", "Study", "med_TL")]
  # Get the correlations of each of these parameters with credible set size
  ess_cors_cred<-sapply(par_ess_cred, function(x) cor(as.numeric(ess_data_cred[,"cred_size"]), as.numeric(ess_data_cred[,x]), use="complete.obs"))
  
  # Plot these out to pdf
  pdf(width=4, height=4, file="ESS Credible Set correlations.pdf")
  ## Loop to get the correlations of the paramters with the num tax and num chars
  for(i in 1:length(par_ess_cred)){
    plot(ess_data_cred[,"cred_size"], ess_data_cred[,par_ess_cred[[i]]], cex=0.01, xlab="Credible set size", ylab=paste(par_ess_cred[[i]], "ESS", sep=" "), main=paste("r =", round(ess_cors_cred[[i]], 3), sep=" ")) # for each parameter, plot it against the credible set size
  }
  dev.off()
  
  ## Compare cred set size to the number of taxa - add ntax into the asdsf and cred size object
  cor(unlist(ess_data_cred[,"ntax"]), ess_data_cred[,"cred_size"], use="complete.obs")
  plot(ess_data_cred[,"ntax"], ess_data_cred[,"cred_size"], cex=0.01)
        ## The above looks horrendous - split it up by dataset and see if it looks any better
        # Amniotes
        amn_asdsf_cred_ntax<-ess_data_cred[which(ess_data_cred[,"Study"]=="Amniotes"),]
        amn_cor_ntax_cred<-cor(unlist(amn_asdsf_cred_ntax[,"ntax"]), amn_asdsf_cred_ntax[,"cred_size"])
        plot(amn_asdsf_cred_ntax[,"ntax"], amn_asdsf_cred_ntax[,"cred_size"], cex=0.2, xlab="N tax", ylab="Cred set")
        
        amn_cor_cred_asdsf<-cor(amn_asdsf_cred_ntax[,"cred_size"], amn_asdsf_cred_ntax[,"asdsf"])
        plot(amn_asdsf_cred_ntax[,"cred_size"], amn_asdsf_cred_ntax[,"asdsf"], cex=0.2, xlab="Cred set", ylab="asdsf")
        
        ## DO NOT HAVE CRED SIZE FOR mtDNA or Barcoding
        
        ## Phylota
        ph_asdsf_cred_ntax<-ess_data_cred[which(ess_data_cred[,"Study"]=="PhyLoTA"),]
        ph_cor_ntax_cred<-cor(unlist(ph_asdsf_cred_ntax[,"ntax"]), ph_asdsf_cred_ntax[,"cred_size"])
        plot(ph_asdsf_cred_ntax[,"ntax"], ph_asdsf_cred_ntax[,"cred_size"], cex=0.01, xlab="N tax", ylab="Cred set")
        # compare asdsf and cred size for Phylota only
        ph_cor_ntax_cred<-cor(ph_asdsf_cred_ntax[,"cred_size"], ph_asdsf_cred_ntax[,"asdsf"])
        plot(ph_asdsf_cred_ntax[,"cred_size"], ph_asdsf_cred_ntax[,"asdsf"], cex=0.01, xlab="Cred size", ylab="ASDSF")
        # This hits the maximum cred size pretty quickly
            
  ############# Process chain swaps for some correlation
  # Set up cswaps formatted in a similar way as cred set, TL for merge with ESS, etc
  cswaps<-lapply(diag_cswaps, function(x) x$chain_swaps) # get out chain swaps
  cswaps12<-unlist(lapply(cswaps, function(x) lapply(x, function(y) y[1,2])))  ## Get just the swap acceptance between chains 1 and 2
  cswap_cor_form<-cbind.data.frame(names(cswaps12), as.numeric(cswaps12), stringsAsFactors = FALSE)
  colnames(cswap_cor_form)<-c("analysis_name", "cswap_12")
  # loop like above to tack on a chain index for each chain within an analysis
  cswap_names_pre<-gsub(".$", "", names(cswaps12)) # names for each
  for(i in unique(cswap_names_pre)){ # for the unique cswap names names:
    index<-which(cswap_names_pre==i) # get an index of chains matching that name
    cswap_cor_form$analysis_name[index]<-paste(cswap_names_pre[index], 1:length(index), sep=".") # append a chain index to each chain within an analysis
  }
  ess_data_cswap<-merge(ess_data_cred, cswap_cor_form, by.x="chain_name", by.y="analysis_name", all.x=TRUE) # finally, merge together the ESS and cswap data and everything else

  
  ### Look at TL/number of branches rather than raw TL - this is probably a more informative measure
  tl_branch<-ess_data_cswap[,"med_TL"]/((unlist(ess_data_cswap[,"ntax"])-1)*2)  # divide the tree length by (#taxa-1)*2  (this is the number of branches)
  ess_data_tl_branch<-cbind(ess_data_cswap, tl_branch)
  # Make a list of the different parameter to be correlated
  par_ess_TL_branch<-colnames(ess_data_tl_branch)[!colnames(ess_data_tl_branch) %in% c("analysis_name", "chain_name", "cswap_12", "cred_size", "ntax", "nchar", "Study", "med_TL", "tl_branch")]
  # Get the correlations of each of these parameters with TL
  ess_cors_TL_branch<-sapply(par_ess_TL_branch, function(x) cor(ess_data_tl_branch[,"tl_branch"], ess_data_tl_branch[,x], use="complete.obs"))
  
  
  # Plot these out to pdf
  pdf(width=4, height=4, file="ESS TLperBranch correlations.pdf")
  ## Loop to get the correlations of the paramters with the num tax and num chars
  for(i in 1:length(par_ess_TL_branch)){
    plot(ess_data_tl_branch[,"tl_branch"], ess_data_tl_branch[,par_ess_TL_branch[[i]]], cex=0.01, xlab="TL/n branches", ylab=paste(par_ess_TL_branch[[i]], "ESS", sep=" "), main=paste("r =", round(ess_cors_TL_branch[[i]], 3), sep=" ")) # for each parameter, plot it against the TL
  }
  dev.off()

  
#### Look at relationship between cred set size and chain swap acceptance rate
  cor_cred_cswaps<-cor(ess_data_cswap[,"cswap_12"], ess_data_cswap[,"cred_size"], use="complete.obs")
  # Plot out the relationship to pdf
  pdf(width=4, height=4, file="Cred_set_Chains_swap_cors.pdf")
  plot(ess_data_cswap[,"cswap_12"], ess_data_cswap[,"cred_size"], cex=0.01, xlab="Chain 1&2 acceptance", ylab="Cred set size", main=paste("corr =", round(cor_cred_cswaps, 3), sep=" "))
  dev.off()
  
#####################################################################################  
#####################################################################################  
###      Get out 1,000 of each type of analysis that we want to reestimate
#####################################################################################
#####################################################################################    

  # Which chains that fail topo ESS pass LnL ESS?
  # Combine low topo objects with 2 & 4 chains
  all_low_topo_ess<-rbind.fill(low_topo_ess_2_chains, low_topo_ess_4_chains)
  
  fail_topo_good_lnlESS<-as.character(all_low_topo_ess[which(all_low_topo_ess$LnL>200),][,1]) # fail topo ESS, pass LnL ESS
  fail_topo_bad_lnlESS<-as.character(all_low_topo_ess[which(all_low_topo_ess$LnL<200),][,1]) # fail topo ESS, fail LnL ESS
  
  # Find chains that fail ESS and see what's up in some of these
  fail_oneormore_ESS<-rownames(pass_fail_combined[pass_fail_combined[,"Pass ESS"]==FALSE,])
  
  ## Find some chains that fail asdsf
  fail_asdsf<-rownames(pass_fail_combined[pass_fail_combined[,"Pass ASDSF"]==FALSE,])
  
  
  
## List these out for ease of use:
  # fail_g_pass_common # Fail Geweke's but pass common diagnostics
  fail_topo_good_lnlESS # Fail topo but have good LnL ESS
  fail_topo_bad_lnlESS # Fail topo and LnL ESS
  fail_oneormore_ESS # Fail one of more parameter ESS
  fail_asdsf # fail asdsf 
  
  
  # ## Right now these above all refer to single chains, and so an analysis could be in each of these lists twice if both chains failed, to make sure that we have 1000 analyses to reanalyses, get the base analysis name and only the unique ones
  fail_topo_good_lnlESS<-gsub("\\.nex.*$|\\.run.*$", "", fail_topo_good_lnlESS)  # strip off the run info
  fail_topo_good_lnlESS<-unique(fail_topo_good_lnlESS)  # get only the unique entries here
  
  fail_topo_bad_lnlESS<-gsub("\\.nex.*$|\\.run.*$", "", fail_topo_bad_lnlESS)  # strip off the run info
  fail_topo_bad_lnlESS<-unique(fail_topo_bad_lnlESS)  # get only the unique entries here

  fail_oneormore_ESS<-gsub("\\.nex.*$|\\.run.*$", "", fail_oneormore_ESS)  # strip off the run info
  fail_oneormore_ESS<-unique(fail_oneormore_ESS)  # get only the unique entries here
  
  fail_asdsf<-gsub("\\.nex.*$|\\.run.*$", "", fail_asdsf)  # strip off the run info
  fail_asdsf<-unique(fail_asdsf)  # get only the unique entries here
  
  
  
    
  # # Sample 1,000 analyses out from each of these for each that have more than 1000 analyses - if not, then just use all
  if(length(fail_topo_good_lnlESS)>1000){
    fail_topo_good_lnlESS_LARGE<-sample(fail_topo_good_lnlESS, 1000)
  }else{
    fail_topo_good_lnlESS_LARGE<-fail_topo_good_lnlESS
  }
  
  if(length(fail_topo_bad_lnlESS)>1000){
    fail_topo_bad_lnlESS_LARGE<-sample(fail_topo_bad_lnlESS, 1000)
  }else{
    fail_topo_bad_lnlESS_LARGE<-fail_topo_bad_lnlESS
  }

  if(length(fail_oneormore_ESS)>1000){
    fail_oneormore_ESS_LARGE<-sample(fail_oneormore_ESS, 1000)    
  }else{
    fail_oneormore_ESS_LARGE<-fail_oneormore_ESS
  }
  
  if(length(fail_asdsf)>1000){
    fail_asdsf_LARGE<-sample(fail_asdsf, 1000)
  }else{
    fail_asdsf_LARGE<-fail_asdsf
  }
  

  # # Write a csv file containing the lines for each of these
  # # This has already been done, so commented out
  # write.csv(fail_topo_good_lnlESS_LARGE, file="fail_topo_good_LnL_ESS_LARGE.csv")
  # write.csv(fail_topo_bad_lnlESS, file="fail_topo_bad_LnL_ESS_LARGE.csv")
  # write.csv(fail_oneormore_ESS_LARGE, file="fail_oneormore_ESS_LARGE.csv")
  # write.csv(fail_asdsf_LARGE, file="fail_ASDSF_LARGE.csv")

    
  ## Here, we'll read in the csv's
  set_1<-read.csv("fail_ASDSF_LARGE.csv", stringsAsFactors = FALSE)
  set_2<-read.csv("fail_oneormore_ESS_LARGE.csv", stringsAsFactors = FALSE)
  set_3<-read.csv("fail_topo_bad_LnL_ESS_LARGE.csv", stringsAsFactors = FALSE)
  set_4<-read.csv("fail_topo_good_LnL_ESS_LARGE.csv", stringsAsFactors = FALSE)
  
  
  ## Write out an R data file that has the diagnostic info on these subsets of analyses
  all_to_reestimate<-c(set_1[,2], set_2[,2], set_3[,2], set_4[,2])
  orig_diags_subsets_LARGE<-diagnosticsAll[unique(all_to_reestimate)]  # the original diagnostics for the datasets that we reestimated trees for with different heating and with nst=mixed
  # get out the pass/fail table for these same sets of analyses
  to_pull_LARGE<-which(gsub("\\.nex.*$|\\.run.*$", "", rownames(pass_fail_combined)) %in% unique(all_to_reestimate))  # get out the indices of the rows to pull out of the combined pass/fail table
  pass_fail_orig_reestimated_LARGE<-pass_fail_combined[to_pull_LARGE,]
  # Identify parameters that failed ESS in each chain
  orig_diags_subsets_LARGE1<-orig_diags_subsets_LARGE # strip the names off the orig diags so the next step doesn't end up with analysis name pasted onto analysis+chain name 
  names(orig_diags_subsets_LARGE1)<-NULL
  failed_ESS_orig_reest_pre_LARGE<-unlist(lapply(orig_diags_subsets_LARGE1, function(x) x$Params_ESS_less_200), recursive=FALSE)
  failed_ESS_orig_reest_LARGE<-lapply(failed_ESS_orig_reest_pre_LARGE, names)   # get the actual names of these parameters as the value in the object
  
  # save these as a list
  orig_diags_LARGE<-list(pass_fail_orig_reestimated_LARGE, failed_ESS_orig_reest_LARGE, orig_diags_subsets_LARGE)
  names(orig_diags_LARGE)<-c("pass_fail", "Pars_failed_ESS", "raw_diags")
  save(orig_diags_LARGE, file="Orig_diags_for_reestimated_LARGE.RData")
  
    ## Get the chain swap acceptance for each of these so that we can adjust heating accordingly

  ## Need to get these out from the orig_diags_subsets_LARGE object - don't have these values at present, use otherwise same code as above for only analyses passing all diags  
  failed_cswapsSub<-which(sapply(orig_diags_subsets_LARGE, function(x) "No log file found" %in% x$chain_swaps[[2]]))
  if(length(failed_cswapsSub)>0){ 
    diag_cswapsSub<-orig_diags_subsets_LARGE[-failed_cswapsSub]  
  }else{
    diag_cswapsSub<-orig_diags_subsets_LARGE
  }  
  cswapsSub<-lapply(diag_cswapsSub, function(x) x$chain_swaps) # get out chain swaps
    cswaps12Sub<-unlist(lapply(cswapsSub, function(x) lapply(x, function(y) y[1,2])))  ## Get just the swap acceptance between chains 1 and 2
    cswaps_for_reest<-cbind.data.frame(names(cswaps12Sub), as.numeric(cswaps12Sub), stringsAsFactors = FALSE)
    colnames(cswaps_for_reest)<-c("analysis_name", "cswap_12")
    cswaps_for_reest$analysis_name<-gsub(".$", "", names(cswaps12Sub))
    # This so far includes a separate chain swap acceptance rate for each chain, just get the mean of these for each analysis
    cswap_mean_reest<-sapply(unique(cswaps_for_reest$analysis_name), function(x) mean(cswaps_for_reest[cswaps_for_reest$analysis_name==x, 2]))
    ## specify if we want to increase or decrease heating for each analysis - this is coarse, but let's decrease if less than or equal to 0.5 and increase otherwise
    heating_change<-cswap_mean_reest
    heating_change[cswap_mean_reest>0.5]<-"increase"
    heating_change[cswap_mean_reest<=0.5]<-"decrease"
    save(cswap_mean_reest, file="example_chains_cswaps_LARGE.RData")
    save(heating_change, file="example_chains_heatChange.RData")
  
    
## Same thing for anlayses with I+G and high correlations between I and G for reestimation with only G  
  ## Pull out studies with a high I+G correlation that either have bad ESS values or that have all good ESS - run these through without +I and see what changes
    ess_All<-unlist(lapply(diagnosticsAll, function(x) x$ESS_values), recursive=FALSE)  # Get ESS values for all analyses (not just those that pass all diagnostic thresholds)
    names(ess_All)<-unlist(lapply(diagnosticsAll, function(x) names(x$"ACT params"))) # assign names to the object
    ess_pinvar_alphaAll<-do.call(rbind, lapply(ess_All, function(x) x[c("alpha", "pinvar")]))   # Get out ESS values for only pinvar and alpha
    colnames(ess_pinvar_alphaAll)<-c("alpha", "pinvar")
    
    ## Get the correlations of I & G for all analyses, including those that fialed a diagnostic
    diagnosticsAll1<-diagnosticsAll # Do this to make an object we can clean names off of so that the subsequent unlisting gives us just the names of each chain, not analysisname.chainname 
    names(diagnosticsAll1)<-NULL
    param_corsAll<-unlist(lapply(diagnosticsAll1, function(x) lapply(x$Param_corrs, reshape2::melt)), recursive=FALSE)
    
    ## Get just Pinvar and Alpha correlations, these are all we need here
    has_pinvaralpha<-lapply(param_corsAll, function(x) "pinvar"  %in% c(as.character(x[,1]), as.character(x[,2])) && "alpha"  %in% c(as.character(x[,1]), as.character(x[,2])))
    param_corAllPA<-param_corsAll[unlist(has_pinvaralpha)]
    
    alpha_pinvar_cors<-vector()
    for(j in 1:length(param_corAllPA)){
      param_corAllPA[[j]][,1]<-as.character(param_corAllPA[[j]][,1])
      param_corAllPA[[j]][,2]<-as.character(param_corAllPA[[j]][,2])
      alpha_pinvar_cors[[j]]<-param_corAllPA[[j]][grep("pinvar alpha", paste(param_corAllPA[[j]][,1], param_corAllPA[[j]][,2])),3]
    }
    names(alpha_pinvar_cors)<-names(param_corAllPA)
    
    
    ess_high_ap_corAll<-ess_pinvar_alphaAll[names(which(alpha_pinvar_cors>0.6)),]  # Pull out the ESS values for analyses with a pinvar/alpha correlation of >0.6
    ess_low_ap_corAll<-ess_pinvar_alphaAll[names(which(alpha_pinvar_cors<0.6)),]  # Pull out the ESS values for analyses with a pinvar/alpha correlation of <0.6
    
    
  all_high_IG_cor_bad_I_ESSAll<-rownames(ess_high_ap_corAll[which(ess_high_ap_corAll[,"alpha"]<200),])  # get out the analyses that have high correlations between IG and also ESS values <200 for G
  all_high_IG_cor_bad_G_ESSAll<-rownames(ess_high_ap_corAll[which(ess_high_ap_corAll[,"pinvar"]<200),])  # get out the analyses that have high correlations between IG and also ESS values <200 for I
  all_high_IG_cor_bad_ESSAll<-c(all_high_IG_cor_bad_I_ESSAll, all_high_IG_cor_bad_G_ESSAll)
  
  pass_fail_common_diags<-pass_fail_combined[,c("Pass ESS", "Pass ASDSF", "Pass PSRF", "Pass Approx topo ESS")]  # make an object that contains the pass/fail stats for the common diagnostics (e.g., exclude pseudo ESS, SF Correlations) - leave out the multichain diagnostics that were calculated for only 2 chains in the 4 chain analyses, too
  pass_common_diags<-names(which(apply(pass_fail_common_diags, 1, function(x) !(FALSE %in% x))))  # Get a list of the analyses that pass all common diags 
  high_IG_cor_good_conv_both_chains<-rownames(ess_high_ap_corAll)[which(rownames(ess_high_ap_corAll) %in% pass_common_diags)] # which of the analyses that converge well (as far as we can tell) have high IG correlations
  high_IG_cor_good_conv_both_chains<-gsub("\\.nex.*$|\\.run.*$", "", high_IG_cor_good_conv_both_chains) # strip off file extension stuff
  high_IG_cor_good_conv<-high_IG_cor_good_conv_both_chains[which(!duplicated(high_IG_cor_good_conv_both_chains))] ## high_IG_cor_good_conv_both_chains looks at each chain independently, i.e., if one chain passes all ESS, but the other doesn't, then it can be included - instead remove any analyses that don't pass common diagnostics for both chains
  
  ## Strip off run/file extensions and get unique entries only for all_high_IG_cor_bad_ESS
  all_high_IG_cor_bad_ESSAll<-gsub("\\.nex.*$|\\.run.*$", "", all_high_IG_cor_bad_ESSAll)  # strip off the run info
  all_high_IG_cor_bad_ESSAll<-unique(all_high_IG_cor_bad_ESSAll)  # get only the unique entries here

## sample out 1000 analyses each that have high IG correlation and either good or bad ESS
  if(length(all_high_IG_cor_bad_ESSAll)>1000){
    high_IG_cor_bad_ESS_LARGE<-sample(all_high_IG_cor_bad_ESSAll, 1000)
  }else{
    high_IG_cor_bad_ESS_LARGE<-all_high_IG_cor_bad_ESSAll 
  }
  
  if(length(high_IG_cor_good_conv)>1000){
    high_IG_cor_good_conv_LARGE<-sample(high_IG_cor_good_conv, 1000)
  }else{
    high_IG_cor_good_conv_LARGE<-high_IG_cor_good_conv
  }
  
  
  # Write a csv file containing the lines for each of these
  # # I've done this already, so commented out for now
  # write.csv(high_IG_cor_bad_ESS_LARGE, file="high_IG_cor_bad_ESS_LARGE.csv")
  # write.csv(high_IG_cor_good_conv_LARGE, file="high_IG_cor_good_conv_LARGE.csv")
  
  ## Here, we'll read in the csv's
  IGset_1<-read.csv("high_IG_cor_bad_ESS_LARGE.csv", stringsAsFactors = FALSE)[,2]
  IGset_2<-read.csv("high_IG_cor_good_conv_LARGE.csv", stringsAsFactors = FALSE)[,2]
  
  ## Write out an R data file of the diagnostics for these IG subset analyses
  datasets_IG_reestimated_LARGE<-c(IGset_1, IGset_2)
  orig_diags_IG_subsets_LARGE<-diagnosticsAll[unique(datasets_IG_reestimated_LARGE)]  # the original diagnostics for the datasets that we reestimated trees for with different heating and with nst=mixed
  # get out the pass/fail table for these same sets of analyses
  to_pull_IG_LARGE<-which(gsub("\\.nex.*$|\\.run.*$", "", rownames(pass_fail_combined)) %in% unique(datasets_IG_reestimated_LARGE))  # get out the indices of the rows to pull out of the combined pass/fail table
  pass_fail_IG_orig_reestimated_LARGE<-pass_fail_combined[to_pull_IG_LARGE,]
  # Identify parameters that failed ESS in each chain
  failed_ESS_orig_IG_reest_pre_LARGE<-unlist(lapply(orig_diags_IG_subsets_LARGE, function(x) x$Params_ESS_less_200), recursive=FALSE)
  failed_ESS_orig_IG_reest_LARGE<-lapply(failed_ESS_orig_IG_reest_pre_LARGE, names)   # get the actual names of these parameters as the value in the object
  
  # save these as a list
  orig_IG_diags_LARGE<-list(pass_fail_IG_orig_reestimated_LARGE, failed_ESS_orig_IG_reest_LARGE, orig_diags_IG_subsets_LARGE)
  names(orig_IG_diags_LARGE)<-c("pass_fail", "Pars_failed_ESS", "raw_diags")
  save(orig_IG_diags_LARGE, file="Orig_diags_for_IG_reestimated.RData")


##########################################################################################  
##########################################################################################
##########################################################################################
## Manually tweak some analyses that perform poorly to try to improve convergence
##########################################################################################

  # of the subset of ESS analyses, pull out those that also fail asdsf
  fail_ess_and_asdsf<-fail_oneormore_ESS_LARGE[which(fail_oneormore_ESS_LARGE %in% fail_asdsf)]
  num_failed_ESS_orig<-sapply(failed_ESS_orig_reest_LARGE, length) # from failed_ESS_orig_reest_LARGE get the number of failed ESS parameters for analyses
  ana_names_failedESS<-gsub("\\.nex.*$|\\.run.*$", "", names(num_failed_ESS_orig)) # get names of the analyses (not individual chains) for each element of num_failed_ESS_orig
  num_failed_ESS_w_names<-cbind(ana_names_failedESS, num_failed_ESS_orig) # tack these on as an extra column
  num_fail_ESS_failESSasdaf<-num_failed_ESS_w_names[num_failed_ESS_w_names[,"ana_names_failedESS"] %in% fail_ess_and_asdsf,] # extract the info only for analyses in fail_ess_and_asdsf
  unique_analyses_forNumESS<-unique(num_fail_ESS_failESSasdaf[,"ana_names_failedESS"]) # unique analyses in this object
  total_failed_ESS<-vector(length=length(unique_analyses_forNumESS)) # set up a vector to drop total number of failed params into 
  names(total_failed_ESS)<-unique_analyses_forNumESS # name the vector with analysis names
  for(i in unique_analyses_forNumESS){ # for each of the unique analyses in the object, get the total number of params that fail ESS 
    index<-which(num_fail_ESS_failESSasdaf[,"ana_names_failedESS"]==i)
    total_failed_ESS[i]<-sum(as.numeric(num_fail_ESS_failESSasdaf[index,"num_failed_ESS_orig"]))
  }
  manual_reanalyze<-sort(total_failed_ESS, decreasing=TRUE)[1:20] # get out the 20 analyses that fail the most ESS values

# Write a csv file containing the names of these analyses
## I've done this already, so commented out for now
  # write.csv(names(manual_reanalyze), file="Manually_reanalyze.csv")
  
  
## Read in a set of analyses that were pulled out manually previously
manual_reanalyze<-read.csv("Manually_reanalyze.csv", stringsAsFactors=FALSE)[,2]

  # Get the original diagnostics
  orig_manual_raw_diags<-diagnosticsAll[manual_reanalyze]
  
  ## Write out the diags in the same way as done for the large sets of reestimates
  # Get pass/fail table
  to_pull_manual<-which(gsub("\\.nex.*$|\\.run.*$", "", rownames(pass_fail_combined)) %in% names(orig_manual_raw_diags))  # get out the indices of the rows to pull out of the combined pass/fail table
  pass_fail_orig_manual<-pass_fail_combined[to_pull_manual,]
  # Identify parameters that failed ESS in each chain
  failed_ESS_orig_manual_pre<-unlist(lapply(orig_manual_raw_diags, function(x) x$Params_ESS_less_200), recursive=FALSE)
  failed_ESS_orig_manual<-lapply(failed_ESS_orig_manual_pre, names)   # get the actual names of these parameters as the value in the object
  
  # save these as a list
  setwd(write_to)
  orig_manual_diags<-list(pass_fail_orig_manual, failed_ESS_orig_manual, orig_manual_raw_diags)
  names(orig_manual_diags)<-c("pass_fail", "Pars_failed_ESS", "raw_diags")
  save(orig_manual_diags, file="orig_manual_diags.RData")


####################################################################################################
####################################################################################################
############  PCA of convergence diagnostics to see what's up
####################################################################################################
####################################################################################################  

  # Get data all set up to run a PCA

  ### Add in autocorrelation times to object with ESS values, etc.
  act_table_PCA<-act_table
  colnames(act_table_PCA)<-paste(colnames(act_table), "ACT", sep=".")
  data_for_pca_1<-merge(ess_data_tl_branch, act_table_PCA, by.x="chain_name", by.y="row.names", all.x=TRUE) 


  ## Get out the PRSF values and merge with the object with other relevant info - use PSRF from 4 chains for analyses that have them
  psrf_2<-lapply(diagnostics[num_chains==2], function(x) x$PSRF$psrf[,1])
  psrf_4<-lapply(diagnostics[num_chains==4], function(x) x$PSRF$PSRF$psrf[,1])
  psrf_1<-c(psrf_2, psrf_4)
  psrf_names<-unique(unlist(lapply(psrf_1, names)))    ## Different analyses included different parameters (e.g., GTR vs HKY), need to find the full list of parameters that were used
  psrf<-do.call(rbind, lapply(psrf_1, function(x) x[match(psrf_names, names(x))]))
  colnames(psrf)<-paste(psrf_names, "psrf", sep=".")
  data_for_pca_2<-merge(data_for_pca_1, psrf, by.x="analysis_name", by.y="row.names", all.x=TRUE)  # merge together
  colnames(data_for_pca_2)[colnames(data_for_pca_2) %in% c("LnL", "TL", "r.A...C.", "r.A...G.", "r.A...T.", "r.C...G.", "r.C...T.", "r.G...T.", "pi.A.", "pi.C.", "pi.G.", "pi.T.", "alpha", "pinvar", "topo", "kappa", "LnPr")]<-paste(colnames(data_for_pca_2)[colnames(data_for_pca_2) %in% c("LnL", "TL", "r.A...C.", "r.A...G.", "r.A...T.", "r.C...G.", "r.C...T.", "r.G...T.", "pi.A.", "pi.C.", "pi.G.", "pi.T.", "alpha", "pinvar", "topo", "kappa", "LnPr")], "ESS", sep=".")
  colnames(data_for_pca_2)[colnames(data_for_pca_2)=="asdsf_for_TL"]<-"asdsf"
  
  ## add in move acceptance info
  ## Attach analysis name and chain info within each element of the list of acceptances for 2 chains
  acceptance_only_2_chains_prep<-list()
  for(i in 1:length(acceptance_only_2_chains)){
    analysis_name<-rep(names(acceptance_only_2_chains)[[i]], 2)
    chain_name<-paste0(rep(names(acceptance_only_2_chains)[[i]], 2), ".",c(1:2))
    acceptance_only_2_chains_prep[[i]]<-cbind(analysis_name, chain_name, acceptance_only_2_chains[[i]])
  }
  ## Attach analysis name and chain info within each element of the list of acceptances for 4 chains
  acceptance_4_chains_prep<-list()
  for(i in 1:length(acceptance_4_chains)){
    analysis_name<-rep(names(acceptance_4_chains)[[i]], 2)
    chain_name<-paste0(rep(names(acceptance_4_chains)[[i]], 2), ".",c(1:4))
    acceptance_4_chains_prep[[i]]<-cbind(analysis_name, chain_name, acceptance_4_chains[[i]])
  }
  
  acc_2_chains_only_table<-do.call(rbind.fill, lapply(acceptance_only_2_chains_prep, as.data.frame))
  acc_4_chains_table<-do.call(rbind.fill, lapply(acceptance_4_chains_prep, as.data.frame))
  acc_all_labelled<-rbind.fill(acc_2_chains_only_table, acc_4_chains_table)
  acc_all_labelled_comb<-acc_all_labelled[,which(!colnames(acc_all_labelled) %in% "analysis_name")]  # for easier combining, let's actually remove the column with analysis name
  acc_all_labelled_comb$chain_name<-as.character(acc_all_labelled_comb$chain_name)
              
  # Convert all to numeric
  indx <- sapply(acc_all_labelled_comb, is.factor)
  acc_all_labelled_comb[indx] <- lapply(acc_all_labelled_comb[indx], function(x) as.numeric(as.character(x)))
              
              
  # combine acceptance rates with other stuff for PCA & mult reg
  data_for_pca_3<-merge(data_for_pca_2, acc_all_labelled_comb, by.x="chain_name", by.y="chain_name", all.x=TRUE)
  
  
  ## Trim out the chain and analysis name stuff from this before feeding it into the PCA
  data_for_pca<-data_for_pca_3[,which(!names(data_for_pca_3) %in% c("chain_name", "analysis_name"))]
  data_for_pca$ntax<-as.numeric(data_for_pca$ntax) # fix up structure issues here
  data_for_pca$nchar<-as.numeric(data_for_pca$nchar) # fix up structure issues here
  data_for_pca_phylota<-data_for_pca[which(data_for_pca[,"Study"]=="PhyLoTA"),] # make an object that includes only the Phylota analyses
  
  ## PCA of only convergence diagnostics
  ## separate out core diagnostics of convergence from dataset properties and things
   ## that are properties of the MCMC itself (move acceptance, etc.)
  diags_to_include<-c("LnL.ESS", "TL.ESS", "r.A...C..ESS", "r.A...G..ESS",  "r.A...T..ESS", "r.C...G..ESS", 
                      "r.C...T..ESS", "r.G...T..ESS",  "pi.A..ESS",  "pi.C..ESS", "pi.G..ESS", "pi.T..ESS",
                      "alpha.ESS", "pinvar.ESS", "topo.ESS", "kappa.ESS", "LnPr.ESS",  "asdsf",
                      "LnL.psrf", "TL.psrf", "r.A...C..psrf", "r.A...G..psrf", "r.A...T..psrf", "r.C...G..psrf",
                      "r.C...T..psrf", "r.G...T..psrf", "pi.A..psrf", "pi.C..psrf", "pi.G..psrf", "pi.T..psrf",
                      "alpha.psrf", "pinvar.psrf", "kappa.psrf",  "LnPr.psrf", "topo.ACT", "LnL.ACT", "TL.ACT", 
                      "r.A...C..ACT", "r.A...G..ACT", "r.A...T..ACT", "r.C...G..ACT", "r.C...T..ACT", "r.G...T..ACT",
                      "pi.A..ACT", "pi.C..ACT", "pi.G..ACT", "pi.T..ACT", "alpha.ACT", "pinvar.ACT"
  )
  
  diags_PCA_no_ACT<-c("LnL.ESS", "TL.ESS", "r.A...C..ESS", "r.A...G..ESS",  "r.A...T..ESS", "r.C...G..ESS", 
                      "r.C...T..ESS", "r.G...T..ESS",  "pi.A..ESS",  "pi.C..ESS", "pi.G..ESS", "pi.T..ESS",
                      "alpha.ESS", "pinvar.ESS", "topo.ESS", "kappa.ESS", "LnPr.ESS",  "asdsf",
                      "LnL.psrf", "TL.psrf", "r.A...C..psrf", "r.A...G..psrf", "r.A...T..psrf", "r.C...G..psrf",
                      "r.C...T..psrf", "r.G...T..psrf", "pi.A..psrf", "pi.C..psrf", "pi.G..psrf", "pi.T..psrf",
                      "alpha.psrf", "pinvar.psrf", "kappa.psrf",  "LnPr.psrf"
  )
  
  
  
  diags_PCA_data<-data_for_pca[,diags_to_include]
  diags_PCA_data_noACT<-data_for_pca[,diags_PCA_no_ACT]
  # diags_PCA_data_ALL<-data_for_pca[,diags_ALL]
  # diags_PCA_data_ALL_noACT<-data_for_pca[,diags_ALL_noACT]
  

  ################################################################################################
  ################################################################################################                
  ## Run a pca on this 
  ################################################################################################


  resPCAdiags_noACT<-pca(diags_PCA_data_noACT, nPcs=13, method="nipals", scale="uv") 
  
  pdf(file="PCA_Core_diagsnoACT.pdf", width=25, height=25)
  slplot(resPCAdiags_noACT)
  slplot(resPCAdiags_noACT, pcs=c(2,3))
  slplot(resPCAdiags_noACT, pcs=c(3,4))
  plotPcs(resPCAdiags_noACT, pcs=1:4)
  dev.off()
  
  
  
########################################################################
###     Multiple regression
########################################################################
  PC1_core_diags<-resPCAdiags_noACT@scores[,1]  # Get out the scores for PC1
  
  ## Multiple regression with Ntax, Nchars, TL, ACT, cred size, cswaps as predictors
  multregPC1core<-lm(PC1_core_diags ~ data_for_pca[,"ntax"] + data_for_pca[,"nchar"] + data_for_pca[,"tl_branch"])
  summary(multregPC1core)
  
  multregCorePC1_summ<-rbind(summary(multregPC1core)[[4]], c(summary(multregPC1core)[[9]], NA, NA, NA))
  rownames(multregCorePC1_summ)[nrow(multregCorePC1_summ)]<-"adjusted R square"
  write.csv(multregCorePC1_summ, file="mr_Core_PC1.csv")
  
  # take a look into PC2
  PC2_core_diags<-resPCAdiags_noACT@scores[,2]
  mr_core_PC2<-lm(PC2_core_diags ~ data_for_pca[,"ntax"] + data_for_pca[,"nchar"] + data_for_pca[,"tl_branch"])
  summary(mr_core_PC2)

          mr_core_PC2<-rbind(summary(mr_core_PC2)[[4]], c(summary(mr_core_PC2)[[9]], NA, NA, NA))
          rownames(mr_core_PC2)[nrow(mr_core_PC2)]<-"adjusted R square"
          
          write.csv(mr_core_PC2, file="mr_Core_PC2.csv")


          
      


#########################################################################################################  
#########################################################################################################          
########    PCA of chain acceptances          
#########################################################################################################          
diags_acceptance<-c("Dirichlet(Revmat)", "Slider(Revmat)",
                   "Dirichlet(Pi)", "Slider(Pi)", "Multiplier(Alpha)", "ExtSPR(Tau,V)", "ExtTBR(Tau,V)",
                   "NNI(Tau,V)", "ParsSPR(Tau,V)", "Multiplier(V)", "Nodeslider(V)", "TLMultiplier(V)",
                   "Slider(Pinvar)", "cswap_12"
)
          
# acceptance_PCA_data<-data_for_pca[,c("Study", diags_acceptance)] # used this to then after trimming out NA rows find out which datasets have acceptance info
# acceptance_PCA_data <- acceptance_PCA_data[rowSums(!is.na(acceptance_PCA_data[,diags_acceptance])) > 0,]
# unique(acceptance_PCA_data$Study) # Amniotes does not have info on acceptances

acceptance_PCA_data<-data_for_pca[,diags_acceptance]
acceptance_PCA_data <- acceptance_PCA_data[rowSums(!is.na(acceptance_PCA_data)) > 0,]  # remove rows that are all NA


resPCA_acceptance<-pca(acceptance_PCA_data, nPcs=13, method="nipals", scale="uv") 

pdf(file="PCA_acceptances.pdf", width=25, height=25)
slplot(resPCA_acceptance)
slplot(resPCA_acceptance, pcs=c(2,3))
slplot(resPCA_acceptance, pcs=c(3,4))
plotPcs(resPCA_acceptance, pcs=1:4)
dev.off()


## Save PCA objects for fancier plotting later
save(resPCA_acceptance, resPCAdiags_noACT, acceptance_PCA_data, diags_PCA_data_noACT, file="PCA_results.RData")

### Multiple regression of dataset properties against acceptance rates

PC1_acc<-resPCA_acceptance@scores[,1]  # Get out the scores for PC1

## Multiple regression with Ntax, Nchars, TL, ACT, cred size, cswaps as predictors
multregAcc<-lm(PC1_acc ~ data_for_pca[rownames(acceptance_PCA_data),"ntax"] + data_for_pca[rownames(acceptance_PCA_data),"nchar"] + data_for_pca[rownames(acceptance_PCA_data),"tl_branch"])
summary(multregAcc)

mr_Acceptances_PC1<-rbind(summary(multregAcc)[[4]], c(summary(multregAcc)[[9]], NA, NA, NA))
rownames(mr_Acceptances_PC1)[nrow(mr_Acceptances_PC1)]<-"adjusted R square"
write.csv(mr_Acceptances_PC1, file="mr_Acceptances_PC1.csv")


############################################################
####   Make fancier PCA plots                      #########
############################################################
PCA_diags<-merge(diags_PCA_data_noACT,  scores(resPCAdiags_noACT), by=0)
loadings_diags<-as.data.frame(loadings(resPCAdiags_noACT))

PCA_acc<-merge(acceptance_PCA_data, scores(resPCA_acceptance), by=0)
loadings_acc<-as.data.frame(loadings(resPCA_acceptance))



pdf(file="Fig_4_PCA_diags.pdf", width=10, height=10)
ggplot()+
  geom_point(data=PCA_diags, aes(x=PC1, y=PC2), size=0.01, alpha=0.8, colour="gray")+
  geom_segment(data=loadings_diags*13, aes(x=0, y=0, xend=PC1, yend=PC2)
               , arrow=arrow(length=unit(0.2,"cm")), alpha=0.75) +
  geom_text(data=loadings_diags*13, aes(x=PC1, y=PC2, label=rownames(loadings_diags)), alpha=1, size=3)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous("PC 1 \n 6.5% of variance")+
  scale_y_continuous("PC 2 \n 5.8% of variance")
dev.off()




pdf(file="Fig_5_PCA_acceptances.pdf", width=10, height=10)
ggplot()+
  geom_point(data=PCA_acc, aes(x=PC1, y=PC2), size=0.01, alpha=0.8, colour="gray")+
  geom_segment(data=loadings_acc*13, aes(x=0, y=0, xend=PC1, yend=PC2)
               , arrow=arrow(length=unit(0.2,"cm")), alpha=0.75) +
  geom_text(data=loadings_acc*13, aes(x=PC1, y=PC2, label=rownames(loadings_acc)), alpha=1, size=3)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous("PC 1 \n 53% of variance")+
  scale_y_continuous("PC 2 \n 13% of variance")
dev.off()




####################################################################################################
## final thing to pull out

#### Get out total number of analyses successfully pulled out for each study, including failed analyses
study_and_analysis_pre<-comp_fail6[,c("analysis_name", "ntax", "nchar", "Study")]
study_and_analysis<-study_and_analysis_pre[-which(duplicated(study_and_analysis_pre)),]
length(which(study_and_analysis[,"Study"]=="PhyLoTA"))
length(which(study_and_analysis[,"Study"]=="Amniotes"))
length(which(study_and_analysis[,"Study"]=="Barcoding"))
length(which(study_and_analysis[,"Study"]=="mtDNA"))



  