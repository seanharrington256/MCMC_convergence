## This script contains functions necessary to extract the trees and parameter values for chains
## and calculate several convergence diagnostics and summaries of the trees/parameter values


# Read in chains from a .tar.gz compressed archive - include a part to thin the chain if it has too many samples
# Trim is the thinning interval to be passed to load.multi, target_samp is the number of total samples to shoot for at the end (if there are that many), and max_gens is last generation that we want to include - e.g., if we want to read in only the first 100,000 generations of an analysis, we set this to 100,000
read.compressed.chain<-function(x, trim=1, target_samp=1000, max_gens=10000000){
  untar(x, exdir=sub(".tar.gz", "",  x)) # unzip the folder
  unz_folder<-sub(".tar.gz", "",  x) # get the name of the unzipped folder
  if(list.files(path=unz_folder)==unz_folder){ # some analyses get compressed in a way that when unzipped the folder has another folder of the same name inside that actually has the results in it - handle that here
    chains<-load.multi(paste(unz_folder, unz_folder, sep="/"), trim=trim)
  }else{
    chains<-load.multi(unz_folder, trim=trim) # read in the chains
  }
  chain_length<-max(chains[[1]]$ptable$Gen) # get out the number of generations
  if(chain_length>max_gens){
    end_gen<-which(chains[[1]]$ptable$Gen==max_gens)   # find the row number in the ptable and tree number that corresponds to the maximum generation that we want
    for(i in 1:length(chains)){
      chains[[i]]$trees<-chains[[i]]$trees[1:end_gen]
      chains[[i]]$ptable<-chains[[i]]$ptable[1:end_gen,]
    }
  }
  if(length(chains[[1]]$trees)>target_samp){   # If the number of samples is higher than our target:
    new_thin<-((length(chains[[1]]$trees)-1)/target_samp)*trim  # get a new thinning interval by dividing the number of samples we have (minus 1) by what we want, and then multiply by the thinning interval that was already used
    for(i in 1:length(chains)){
      chains[[i]]$trees<-chains[[i]]$trees[seq(from = 1, to = length(chains[[i]]$trees), by = new_thin)]
      chains[[i]]$ptable<-chains[[i]]$ptable[seq(from = 1, to = nrow(chains[[i]]$ptable), by = new_thin),]
      ### Also need to edit the number of generations/tree in the chains object, since this has now been changed
      chains[[i]]$gens.per.tree<-chains[[i]]$gens.per.tree*new_thin
    }
  }
  unlink(unz_folder, recursive=TRUE)
  return(chains)
}


## Write a trycatch for read.compressed.chain so that everything continues if something goes wrong reading in some chains
read.comp.catch<-function(folder, trim=1, target_samp=1000, max_gens=10000000){
  return(tryCatch(read.compressed.chain(folder, trim=trim, target_samp=target_samp, max_gens=max_gens), error=function(e) {
    fail<-list(paste(folder, "FAILED"), e)
    return(fail)
  }))
}


### Function to get values from the nexus file that follow some specific characters and "=" (e.g., for getting nruns, etc)
### Text refers to the text lines reading from, values is the name of the value (as referenced in the text) to get, which is 
### followed by "=" then the number of interest (e.g., "text" arg could be a nexus file and "value" arg could be ntax)
nums_from_txt<-function(text, values)
{
  line<-grep(values, text, value=TRUE, ignore.case=TRUE)    ## Get out the line with the value of interest
  res<-as.numeric(strapplyc(line, paste(values, " *= *([[:digit:]]+)", sep=""), simplify = TRUE, ignore.case=TRUE))
  return(res)
}



# Get the burnin of a chain
# Similar to Beiko et al. 2006, this function estimates the burnin of a chain by finding the first generation
# at which the likelihood is higher than average of the last 10% of a chain
getburnin<-function(chain) {
  ## Get the LnL and name the values by iteration
  LnL<-chain$ptable[,"LnL"]
  names(LnL)<-1:length(LnL)
  # Find the average LnL of the last 10% of the chain
  mean_end_gens<-mean(LnL[ceiling(length(LnL)*0.9):length(LnL)])
  # Find the first iteration that passes burnin
  burn_iter<-min(which(LnL>mean_end_gens))
  # Find what percentage this equates to (rounding up to a whole percentage)
  burn_percent<-ceiling((burn_iter/length(LnL))*100)
  return(burn_percent)
}

## Modified calc_ess (tracerer package) function to also spit the autocorrelation time - already calculated, might as well spit it out here
calc.ess.act<-function (trace, sample_interval){
  if (!is.numeric(trace)) {
    stop("trace must be numeric")
  }
  if (sample_interval < 1) {
    stop("sample interval must be at least one")
  }
  act <- tracerer::calc_act(trace = trace, sample_interval = sample_interval)
  ess <- length(trace)/(act/sample_interval)
  res<-c(ess, act)
  names(res)<-c("ESS", "ACT")
  return(res)
}



###### Modified makeplot.param function from rwty to pull out just ESS and pararmeter values
# the cutoff argument specifies how far into the chain to include, i.e., if we want to check ess values from burnin to halfway into the chain
# Burnin in also modified to reflect % burnin
# ESS values are calculated using the calc_ess function of Tracerer, same method as in Tracer
# This is also used to return autocorrelation times
param_vals_ess<-function (chains, burnin = 0, cutoff=100, parameter = "LnL") 
{
  chains = check.chains(chains)
  ptable = combine.ptables(chains, burnin)
  # Get the specific values for the parameter of interest, want to output this
  vals<-lapply(chains, FUN = function(x) x$ptable[parameter][(ceiling(((burnin/100)*(length(x$ptable[[parameter]])-1)))+1):ceiling(((cutoff/100)*length(x$ptable[[parameter]]))),])
  if (parameter %in% names(ptable)) {
    ess_auto <- unlist(lapply(vals, FUN = function(x) calc.ess.act(x, 1)))  # Using a sampling interval of 1 here for the calculation of ESS - shouldn't affect the calculation of ESS
    res<-c(ess_auto[[1]], vals, ess_auto[[2]])
    names(res)<-c("ESS", "vals", "ACT")
    return(res)
  }
  else {
    stop(sprintf("The variable '%s' is not a column in the table of parameters you have supplied", 
                 parameter))
  }
}

## Get only parameter values
param_vals<-function(chains, burnin = 0, cutoff=100, parameter = "LnL") 
{
  chains = check.chains(chains)
  ptable = combine.ptables(chains, burnin)
  # Get the specific values for the parameter of interest, want to output this
  vals<-lapply(chains, FUN = function(x) x$ptable[parameter][(ceiling(((burnin/100)*(length(x$ptable[[parameter]])-1)))+1):ceiling(((cutoff/100)*length(x$ptable[[parameter]]))),])
  return(vals)
}



# This function gets the final ESS values (i.e, ESS when cutoff = 100) for all parameters in an analysis
# Takes a single chain as input
# Spits out the aoutocorrelation time and values of parameters at each generation through the chains, also
final_ess_ACT_all<-function(chains, burn){
  params<-colnames(chains$ptable)[2:length(colnames(chains$ptable))] # get the list of parameters, exlcude the first column because this is generation number
  ess_vals_act<-lapply(params, FUN=function(x) param_vals_ess(chains=chains, burn=burn, parameter=x, cutoff=100))
  ess_act_only<- sapply(ess_vals_act, function(x)  x[c(1,3)])
  colnames(ess_act_only)<-params
  par_vals_only<-sapply(ess_vals_act, function(x) x[2])
  names(par_vals_only)<-params
  res<-list(ess_act_only, par_vals_only)
  return(res)
}


## Get the generations to ESS value of 200 using multiples of 10% of chain length
## If at the end of the chain, the ESS value is still <200, returns Inf
gens_to_200<-function(chains, burn, parameter){
  start<-round(burn+10, -1) # start getting ESS at the burnin + 10%, rounded to the nearest 10%
  # Get ess values at each of these cutoffs until the ESS hits 200
  for(i in seq(start, 100, by=10)) {
    ess<-param_vals_ess(chains=chains, burnin=burn, cutoff=i, parameter=parameter)[[1]]
    if(ess>=200){break}
    if(i==100 & ess<200){i<-NA}
  }
  percent<-i
  num_samples<-length(chains[[1]])
  gens_sample<-chains$gens.per.tree
  gens.200<-ceiling(percent*0.01*(num_samples-1))*gens_sample
  return(gens.200)
}

## Modify the internal rwty tree.ess function to use the tracerer function for calculating ESS and autocorrelation time
tree.ess.tracer<-function (tree.list, treedist = "PD") 
{
  i <- sample(1:length(tree.list), 1)
  distances <- rwty:::tree.distances(tree.list, i, treedist = treedist)
  ESS_ACT <- apply(distances, 2, calc.ess.act, sample_interval=1)
  ESS_ACT<-as.numeric(ESS_ACT)
  names(ESS_ACT)<-c("ESS", "ACT")
  return(ESS_ACT)
}

## Modify the internal rwty tree.ess.multi function to use tree.ess.tracer
tree.ess.tracer.multi<-function (tree.list, n = 20, treedist = "PD") 
{
  print(sprintf("Calculating pseudo ESS for %s trees and %s replicates, please be patient", 
                length(tree.list), n))
  data <- replicate(n, tree.ess.tracer(tree.list, treedist))
  return(data)
}

## Edit topological.pseudo.ess to use burnin as a percentage and to include a cutoff, also include part of makeplot.pseudo.ess to summarize the multiple replicates
## also spit out the autocorrelation times
topo.pseudo.ess.ACT.cutoff<-function (chains, burn = 0, cutoff=100, n = 20, treedist = "PD") 
{
  chains <- check.chains(chains)
  chain <- chains[[1]]
  indices <- seq(from = (ceiling((burn*0.01*(length(chains[[1]]$trees)-1)))+1), to = ceiling(cutoff*0.01*length(chains[[1]]$trees)), 
                 by = 1)
  trees <- lapply(chains, function(x) x[["trees"]][indices])
  raw.ess.ACT <- lapply(trees, tree.ess.tracer.multi, n, treedist)
  raw.ess<-lapply(raw.ess.ACT, function(x) x[1,])
  raw.act<-lapply(raw.ess.ACT, function(x) x[2,])
  final.ess <- data.frame(t(matrix(unlist(raw.ess), nrow = length(chains), 
                                   byrow = T)))
  final.act <- data.frame(t(matrix(unlist(raw.act), nrow = length(chains), 
                                   byrow = T)))
  colnames(final.ess) <- names(chains)
  colnames(final.act) <- names(chains)
  dat.ess <- data.frame(median.ess = apply(final.ess, 2, FUN = median), 
                        ci.lower = apply(final.ess, 2, FUN = function(x) quantile(x, 
                                                                                  0.025)), ci.upper = apply(final.ess, 2, FUN = function(x) quantile(x, 
                                                                                                                                                     0.975)), chain = names(chains))
  dat.act <- data.frame(median.act = apply(final.act, 2, FUN = median), 
                        ci.lower = apply(final.act, 2, FUN = function(x) quantile(x, 
                                                                                  0.025)), ci.upper = apply(final.act, 2, FUN = function(x) quantile(x, 
                                                                                                                                                     0.975)), chain = names(chains))
  
  return(list(dat.ess, dat.act))
}


## Get the topo generations to psuedo ESS value of 200 using multiples of 10% of chain length
## If at the end of the chain, the ESS value is still <200, returns NA
## Uses the median estimate of topo ESS
topo_pseudo_gens_to_200<-function(chains, burn=0, n=20){
  start<-round(burn+10, -1) # start getting ESS at the burnin + 10%, rounded to the nearest 10%
  # Get topo ess values at each of these cutoffs until the ESS hits 200
  for(i in seq(start, 100, by=10)) {
    ess<-topo.pseudo.ess.ACT.cutoff(chains, burn, cutoff=i, n, treedist = "PD")[[1]]$median.ess
    if(ess>=200){break}
    if(i==100 & ess<200){i<-NA}
  }
  percent<-i
  num_samples<-length(chains[[1]])
  gens_sample<-chains$gens.per.tree
  gens.200<-ceiling(percent*0.01*(num_samples-1))*gens_sample
  return(gens.200)
}



## Edits to rwty topological.autocorr to add in cutoff
## and to fix the default max.sampling.interval, see comment on that line
## This makes a call to an internal rwty function at present
topological.autocorr.cutoff<-function (chains, burn = 0, cutoff=100, max.sampling.interval = NA, autocorr.intervals = 100,
                                       squared = FALSE, treedist = "PD", use.all.samples = FALSE)
{
  chains = check.chains(chains)
  chain = chains[[1]]
  N = ceiling(cutoff*0.01*length(chains[[1]]$trees))
  burnin = (ceiling((burn*0.01*(length(chains[[1]]$trees)-1))))
  if (is.na(max.sampling.interval)) {
    max.sampling.interval = floor((N - burnin)*0.1)  # Original code from rwty had this as floor((N - burnin) - (0.9 * N)), but that fails with any burnin 10% or above, and only returns the 10% of samples it should if burnin is 0
  }
  indices = seq(from = burnin + 1, to = N, by = 1)
  trees = lapply(chains, function(x) x[["trees"]][indices])
  raw.autocorr = lapply(trees, rwty:::tree.autocorr, max.sampling.interval,
                        autocorr.intervals, squared, treedist, use.all.samples)
  final.autocorr = do.call("rbind", raw.autocorr)
  final.autocorr$chain = unlist(lapply(names(chains), rep,
                                       nrow(raw.autocorr[[1]])))
  rownames(final.autocorr) = NULL
  return(final.autocorr)
}



## Modify topological.approx.ess to take burnin as a percentage and work from burnin to cutoff value
## Also make this return the topological autocorrelation - this is already calculated as 
## part of the function, might as well spit it out here
## This makes a call to an internal rwty function at present
topological.approx.ess.cutoff.ACT<-function (chains, burn = 0, cutoff=100, max.sampling.interval = 100, treedist = "PD",
                                             use.all.samples = TRUE)
{
  chains = check.chains(chains)
  N = ceiling(cutoff*0.01*length(chains[[1]]$trees))
  burnin = (ceiling((burn*0.01*(length(chains[[1]]$trees)-1))))
  if (N - burnin < max.sampling.interval) {
    warning("Not enough trees to use your chosen max.sampling.interval")
    warning("Setting it to 90% of the length of your post-burnin chain instead")
    max.sampling.interval = floor((N - burnin) - (0.1 * N))
  }
  autocorr.intervals = max.sampling.interval
  print(sprintf("Calculating approximate ESS with sampling intervals from 1 to %d",
                max.sampling.interval))
  autocorr.df = topological.autocorr.cutoff(chains, burn, cutoff, max.sampling.interval,
                                            autocorr.intervals, squared = TRUE, treedist = treedist,
                                            use.all.samples = use.all.samples)
  autocorr.m = estimate.autocorr.m(autocorr.df)
  approx.ess.df = rwty:::approx.ess.multi(autocorr.df, autocorr.m, (N-burnin))
  res<-list(approx.ess.df, autocorr.m)
  names(res)<-c("Approx.topo.ESS", "Autocorr.time")
  return(res)
}


## Create a trycatch version of topological.approx.ess.cutoff.ACT so that cases when there are not enough samples to evaluate approx ESS, the function will not stop
## This is currently used only in the context of getting approx ESS without autocorr time, and so autocorr time is not returned here
approx.ess.cutoff.catch<-function(chains, burn = 0, cutoff=100, max.sampling.interval = 100, treedist = "PD",
                                  use.all.samples = TRUE){
  tryCatch(topological.approx.ess.cutoff.ACT(chains, burn, cutoff, max.sampling.interval, treedist, use.all.samples)[[1]], error = function(e) {
    ess<-0
    return(list("FAIL", ess, paste("ESS could not be evaluated for cutoff = ", cutoff, sep=""), cutoff,e))
  }
  )
}


### Get the number of generations to an approximate ESS value of 200
## Uses a trycatch function to get approx ESS so that it doesn't fail when approx ESS cannot be calculated, happens for some small cutoff values
## If at the end of the chain, the ESS value is still <200, returns NA
topo_gens_to_200<-function(chains, burn=0){
  start<-round(burn+10, -1) # start getting ESS at the burnin + 10%, rounded to the nearest 10%
  # Get topo ess values at each of these cutoffs until the ESS hits 200
  failed_cuts<-list()
  for(i in seq(start, 100, by=10)) {
    res<-approx.ess.cutoff.catch(chains, burn, cutoff=i)   #Use the trycatch function for calculating approximate ESS - if approx ESS cannot be evaluated for a given cutoff, ess will be returned as zero, and messages of error and the cutoff that caused the problem will be returned
    ess<-res[[2]]
    if(is.na(ess)){  # for some cutoffs for some datasets, the tree distsances at all intervals are zero (I believe from there being a single tree in these sets) and approx ESS cannot be calculated at these cutoffs - in this case, break and return NA
      i<-NA
      break
    }
    if(ess>=200){break}
    if(i==100 & ess<200){i<-NA}  # if ESS never gets above 100, return NA
    if(res[[1]]=="FAIL"){failed_cuts<-c(failed_cuts, res[c(3,4,5)])}
  }
  percent<-i
  num_samples<-length(chains[[1]])
  gens_sample<-chains$gens.per.tree
  gens.200<-ceiling(percent*0.01*(num_samples-1))*gens_sample
  out<-list(gens.200, failed_cuts)
  return(out)
}




## Function to get the correlations in split frequencies among chains and ASDSF - takes a few pieces from makeplot.splitfreq.matrix and from https://stackoverflow.com/questions/45825685/correlations-for-pairs-of-combinations
get.sf.cors.asdsf<-function (chains, burn = 0) 
{
  chains<-check.chains(chains)
  burnin = (ceiling((burn*0.01*(length(chains[[1]]$trees)-1)))+1)
  dat = rwty:::get.comparison.table(chains, burnin, min.freq = 0)
  asdsf = dat$asdsf
  dat = dat$cladetable
  dat = dat[, (names(dat) %in% names(chains))]
  cr <- cor(dat)
  cr[upper.tri(cr, diag=TRUE)] <- NA
  cr_data_frame<-reshape2::melt(cr, na.rm=TRUE, value.name="cor")
  res<-list(cr_data_frame, asdsf)
  names(res)<-c("sf_corrs", "ASDSF")
  return(res)
}

# Function to convert the parameters of chains to an mcmc class - exclude the generation column - do not want to run diagnostics on this
rwty_to_coda<-function(chains, burn=0, cutoff=100){
  chains = check.chains(chains)
  gens_to_include<-(ceiling(((burn/100)*(length(chains[[1]]$ptable[,1])-1)))+1):ceiling(((cutoff/100)*length(chains[[1]]$ptable[,1])))   #generations between burnin and cutoff
  mcmc.list(lapply(chains, FUN=function(x) mcmc(x$ptable[gens_to_include,!colnames(x$ptable) %in% "Gen"]))) # exclude generation column from the object - not interested in using this to calculate
}

# Function to convert Geweke's z scores to Geweke's p values
# This uses code from https://github.com/mikeryanmay/bonsai/blob/master/package/R/NumericalParameterRefClass.R to convert Geweke's z scores to Geweke's p values
Gewekes.z.to.p<-function(z_scores){
  res<-2 * (1 - pnorm(abs(as.numeric(z_scores))))
  names(res)<-names(z_scores)
  return(res)
}



## Function that takes the lines of a log file (log arg) and the indices of the lines in this object that contain the acceptance rates for a chain (lines arg)
##    returns the acceptance rates for each parameter in the chain
get.acc.rates<-function(log, lines){
  rates<-list()
  for(i in 1:length(lines)){
    cleaned<-gsub("^\\s*", "",  log[lines[[i]]]) # take line in a chain and get rid of the white space at the beginning
    split_out<-strsplit(cleaned, " ") # split the string into different elements separated by spaces
    rate<-as.numeric(split_out[[1]][[1]]) # acceptance rate
    move<-split_out[[1]][[length(split_out[[1]])]]
    names(rate)<-move
    rates[[i]]<-rate
  }
  return(unlist(rates))
}


## Function to get the tree length (both the full list of TL through the chain and median) when supplied some chains and a burnin
get_TL<-function(chains, burnin){
  # Get out the tree length
  TL<-param_vals(chains, burnin = burnin, cutoff=100, parameter = "TL")
  med_TL<-lapply(TL, median)
  names(med_TL)<-rep("median", length(med_TL))
  res<-c(TL, med_TL)
  return(res)
}

## Trycatch version of get_TL
get_TL.catch<-function(chains, burnin){
  return(tryCatch(get_TL(chains, burnin), error=function(e) {
    fail<-list(paste(names(chains), "FAILED"), e)
    return(fail)
  }))
}




# Function takes a compressed analysis that includes chains, nexus file, etc, and reads in various information
# about the analysis as well as the chains
# Read in compressed chains - includes a part to thin the chain if it has too many samples
# Trim is the thinning interval to be passed to load.multi, target_samp is the number of total samples to shoot for at the end (if there are that many), and max_gens is last generation that we want to include - e.g., if we want to read in only the first 100,000 generations of an analysis, we set this to 100,000
get_info_comp<-function(x, trim=1, target_samp=1000, max_gens=10000000){
  untar(x, exdir=sub(".tar.gz", "",  x)) # unzip the folder
  unz_folder<-sub(".tar.gz", "",  x) # get the name of the unzipped folder
  #######################################
  ##### Read in the chains
  #######################################
  if(list.files(path=unz_folder)==unz_folder){ # some analyses get compressed in a way that when unzipped the folder has another folder of the same name inside that actually has the results in it - handle that here
    chains<-load.multi(paste(unz_folder, unz_folder, sep="/"), trim=trim)
    unz_path<-paste(unz_folder, unz_folder, sep="/")
    infile<-paste(unz_path, list.files(path=unz_path, pattern = "\\.nex$", ignore.case=TRUE), sep="/")  # Get the name of the nexus file for down below
    if(length(grep("mbb_", infile))>0){
      infile<-infile[-grep("mbb_", infile)]
    }
  }else{
    chains<-load.multi(unz_folder, trim=trim) # read in the chains
    unz_path<-unz_folder
    infile<-paste(unz_path, list.files(path=unz_path, pattern = "\\.nex$", ignore.case=TRUE), sep="/")  # Get the name of the nexus file for down below
    if(length(grep("mbb_", infile))>0){
      infile<-infile[-grep("mbb_", infile)]
    }
  }
  chain_length<-max(chains[[1]]$ptable$Gen) # get out the number of generations
  if(chain_length>max_gens){
    end_gen<-which(chains[[1]]$ptable$Gen==max_gens)   # find the row number in the ptable and tree number that corresponds to the maximum generation that we want
    for(i in 1:length(chains)){
      chains[[i]]$trees<-chains[[i]]$trees[1:end_gen]
      chains[[i]]$ptable<-chains[[i]]$ptable[1:end_gen,]
    }
  }
  if(length(chains[[1]]$trees)>target_samp){   # If the number of samples is higher than our target:
    new_thin<-((length(chains[[1]]$trees)-1)/target_samp)*trim  # get a new thinning interval by dividing the number of samples we have (minus 1) by what we want, and then multiply by the thinning interval that was already used
    for(i in 1:length(chains)){
      chains[[i]]$trees<-chains[[i]]$trees[seq(from = 1, to = length(chains[[i]]$trees), by = new_thin)]
      chains[[i]]$ptable<-chains[[i]]$ptable[seq(from = 1, to = nrow(chains[[i]]$ptable), by = new_thin),]
      ### Also need to edit the number of generations/tree in the chains object, since this has now been changed
      chains[[i]]$gens.per.tree<-chains[[i]]$gens.per.tree*new_thin
    }
  }
  #########################################################
  ##### Read in the number of characters and taxa 
  #########################################################
  # Identify if there is a bayes block as a separate file from the Nexus file
  bayes_block<-list.files(path=unz_path, pattern = "\\.bayesblock$|\\.bb$|^mbb_", ignore.case=TRUE)
  # If the folder has a bayes block file separate from the main nexus file, then extract the info sepatately from the bayes block and the nexus file. Else, pull everything out of the nexus file.
  if(length(bayes_block)>0){
    ## Read in the nexus file as lines of text
    nex<-readLines(infile)
    # values to be extracted from the nexus file
    from_nex<-c("ntax", "nchar")
    # Extract these values from the file
    ntax_nchar<-sapply(from_nex, FUN=function(x) nums_from_txt(text=nex, values=x))
    ## Read in the mrbayes file to pull out info
    bb<-readLines(paste(unz_path, bayes_block, sep="/"))
    # Values to be extracted from the bayes block
    from_bayes_block<-c("ngen", "nchains", "nruns", "samplefreq")
    # Extract these values from the bayes block
    bb_params<-sapply(from_bayes_block, FUN=function(x) nums_from_txt(text=bb, values=x))
    # Combine the analysis pars from the bayes block and nexus into a single object
    analysis_pars<-c(ntax_nchar, bb_params)
  }else{
    ## Read in the nexus file as lines of text
    nex<-readLines(infile)
    # values to be extracted from the nexus file
    from_nex<-c(" ntax", " nchar", "ngen", "nchains", "nruns", "samplefreq") # use a space on the beginning of ntax and nchar so that when these character strings show up in taxon names (and they do) those don't get mixed in and cause problems
    # Extract these values from the file
    analysis_pars<-sapply(from_nex, FUN=function(x) nums_from_txt(text=nex, values=x))
    names(analysis_pars)<-c("ntax", "nchar", "ngen", "nchains", "nruns", "samplefreq") # assign the names to not 
  }
  ################################################################################################
  ##### Get acceptance rates and chain swap rates for any analyses that have a log file
  ################################################################################################
  logfile<-list.files(path = unz_path,  pattern=".log$") # First check if there is a logfile
  if(length(logfile)>0){
    ### First get out acceptances
    log<-readLines(paste(unz_path, logfile,sep="/"))
    blank_lines<-grep("^\\s*$", log, value=FALSE) # find all of the blank lines
    line_start_chain_1<-grep('Acceptance rates for the moves in the "cold" chain of run 1:', log, value=FALSE, ignore.case=TRUE)+2  # Line starting the table for chain 1 - +2 is because there is this line searched for that designated the chain, then a header, THEN the table starts   
    line_end_chain_1<-min(blank_lines[blank_lines>line_start_chain_1])-1  # end line is the first line of white space only past this, -1
    lines_chain_1<-line_start_chain_1:line_end_chain_1 # Make these lines into an object
    line_start_chain_2<-grep('Acceptance rates for the moves in the "cold" chain of run 2:', log, value=FALSE, ignore.case=TRUE)+2  # same for chain 2 as for chain 1
    line_end_chain_2<-min(blank_lines[blank_lines>line_start_chain_2])-1 
    lines_chain_2<-line_start_chain_2:line_end_chain_2 # Make these lines into an object
    acc_rates_chain_1<-get.acc.rates(log=log, lines=lines_chain_1)  # use get.acc.rates function to get the acceptance rates given the log file and lines that contain the acceptances
    acc_rates_chain_2<-get.acc.rates(log=log, lines=lines_chain_2)
    acceptances<-rbind(acc_rates_chain_1, acc_rates_chain_2)
    ## Then get chains swaps ###################################
    nruns<-length(grep("Chain swap information", log))    ## Find how many independent runs were run by identifying the number of times "Chain swap information" shows up - this will show up once for each run
    ## Make a for loop to iterate over the number of analyses, which may differ across analyses
    chain_swaps<-list()
    for(i in 1:nruns){
      line_start<-grep(paste("Chain swap information for run", i, sep=" "), log, value=FALSE, ignore.case=TRUE)+2  # Line starting the table of chain swaps: +2 is because there is this line searched for that designated the analysis, then a blank space, THEN the table starts   
      chain_indices<-strsplit(log[line_start], " ")[[1]] # the first line of the table is the numbers of the different chains - we need to know how many chains were run before we proceed, so split out values by spaces here
      nchains<-max(chain_indices)
      ## make a for loop to iterate over the number of chains, which may vary across analyese
      rows_acceptance<-list()
      for(j in 1:nchains){
        row_vals<-gsub("^.*\\|\\s*", "",  log[line_start+1+j]) # indexing here is tricky: line start is the header of the table that has the numebr of the chains, +1 gets us to a line of underscores, +j gets us to the desired line in the table - clean out everything before we get to the actual table
        acc_rates<-strsplit(row_vals, " ")[[1]][strsplit(row_vals, " ")[[1]]!=""] # split these out into independent elements and get rid of the blanks
        ## We want to replace the jth element of the list with NA so that the diagonal of the matrix will be NA - first move any indices higher than j up by 1
        acc_rates[which(1:length(acc_rates)>=j)+1]<-acc_rates[1:length(acc_rates)>=j]
        acc_rates[j]<-NA
        rows_acceptance[[j]]<-acc_rates
      }
      chain_swaps[[i]]<-do.call(rbind, rows_acceptance)
    }
  }else{
    acceptances<-list(NA, "No log file found")
    chain_swaps<-list(NA, "No log file found")
  }
  ########################################################################
  ##### Get credible set size for any analyses that have a trprobs file
  ########################################################################
  trprobs_file<-list.files(path=unz_path, pattern="trprobs")
  if(length(trprobs_file)==0){
    cred_size<-list(NA, "No trprobs file found")
  }else{
    trprobs_lines<-readLines(paste(unz_path, trprobs_file, sep="/")) # read in the trprobs file as lines
    tree_lines<-trprobs_lines[grep("P = ", trprobs_lines)] # find the lines that contain "P =" - these are the lines of the distinct trees, and "P =" will be followed by the cumulative probability of all trees to that point 
    cum_probs<-strapplyc(tree_lines, "P = ([[:digit:]]+\\.*[[:digit:]]*)", simplify = TRUE, ignore.case=FALSE) ## extract out the cumulative probabilities from these lines
    cred_size<-which(cum_probs>=0.95)[[1]] ## Find the first tree that is equal to or greater than 0.95. The number of this tree is the number of trees in the credible set
  }
  unlink(unz_folder, recursive=TRUE) # delete the unzipped folder
  res<-list(chains, analysis_pars, acceptances, cred_size, chain_swaps) # prep the output
  names(res)<-c("chains", "analysis_pars", "acc_rates", "cred_size", "chain_swaps")
  return(res)
}

## Write a trycatch for get_info_comp so that everything continues if something goes wrong reading in some chains
read.comp.info.catch<-function(folder, trim=1, target_samp=1000, max_gens=10000000){
  return(tryCatch(get_info_comp(folder, trim=trim, target_samp=target_samp, max_gens=max_gens), error=function(e) {
    fail<-list(paste(folder, "FAILED"), e)
    return(fail)
  }))
}

