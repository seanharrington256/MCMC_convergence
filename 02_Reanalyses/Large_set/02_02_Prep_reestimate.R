### Script to pull out just the necessary files for each of the analyses that we want to reestimate trees for and get them ready to run through MrBayes with new heating or Nst=mixed
library(stringr)
library(gsubfn)

####################################################################################################
##### !!!!!! NOTE THAT directories here are specific to how we ran things on the UH cluster    #####
####################################################################################################


main_dir<-"~/thomson_lts/Sean_storage/MCMC_reanalyses"
chains_dir<-paste0(main_dir, "/all_original_analyses")  # where folders housing the files from the original chains are
processed_dir<-"/home/seanh256/thomson_lts/Sean_storage/MCMC_reanalyses/Prepped_Nex" # where we want to write the processed nexus files - this has to be in the full notation that terminal uses, rather than R, because of the way that this will be called
mb_ready_dir<-paste0(main_dir, "/Nex_new_heat_chains_READY")  # directory for files that are ready for MrBayes with the new heating
mb_ready_dir_nstMixed<-paste0(main_dir, "/Nex_nstMixed_chains_READY")  # directory for files that are ready for MrBayes with nst=mixed
## For IG analyses
chains_dir_IG<-paste0(main_dir, "/IG_original_analyses")  # where folders housing the files from the original chains are
processed_dir_IG<-"/home/seanh256/thomson_lts/Sean_storage/MCMC_reanalyses/Prepped_Nex_IG" # where we want to write the processed nexus files - this has to be in the full notation that terminal uses, rather than R, because of the way that this will be called
mb_ready_dir_IG<-paste0(main_dir, "/Nex_IG_chains_READY")  # directory for files that are ready for MrBayes with the new heating or number of chains


# read in RData file saying to increase or decrease heating
setwd(main_dir)
load("example_chains_heatChange.RData")


### This function is from 01_R_functions_calculate_diagnostics_V2.R, but just copying it here to avoid needing to source in all the functions in that file just to get this one
nums_from_txt<-function(text, values){
  line<-grep(values, text, value=TRUE, ignore.case=TRUE)    ## Get out the line with the value of interest
  res<-as.numeric(strapplyc(line, paste(values, " *= *([[:digit:]]+)", sep=""), simplify = TRUE, ignore.case=TRUE))
  return(res)
}


# Set working directory and list out the .tar.gz archives containing the chains
setwd(chains_dir)
to_extract<-list.files(pattern="*.tar.gz") # list out the .tar.gz folders that will need to be extracted

# Use a loop here to perform the extraction on each of the compressed archives in the to_extract object
for(j in 1:length(to_extract)){
  # list out the files within the archive
  files_in_archive<-untar(to_extract[[j]], list=TRUE)

  # we only want to extract the files ending in .nex, .bb, or .bayesblock
  files_to_extract<-files_in_archive[grep(pattern=".nex$|.bayesblock$|.bb", files_in_archive)]
  
  # Untar the files to the location of the processed nexus file
  if(length(grep("[[:alpha:]].*/.*[[:alpha:]]", files_to_extract)>0)){ # phylota analyses are all coming out in a directory - need to identify this, then go into the directpory and move the file to the containing folder - identify nex nested within its own directory by searching for letters separated from each other by /
    untar(to_extract[[j]], files=files_to_extract, exdir=processed_dir) # untar 
    files_to_extract<-unique(gsub("^\\./|^\\./._", "", files_to_extract)) # drop out ._ as for others
    file.copy(paste(processed_dir, files_to_extract, sep="/"), processed_dir) # copy the file from it's own directory into the directory we want as a sinle nexus file
    unlink(paste(processed_dir, strsplit(files_to_extract, "/")[[1]][[1]], sep="/"), recursive = TRUE) # delete the directory - no longer necessary
    files_to_extract<-strsplit(files_to_extract, "/")[[1]][[2]] # Need to change the files_to_extract object to no longer include the containing directory
  }else{
    untar(to_extract[[j]], files=files_to_extract, exdir=processed_dir)
    files_to_extract<-unique(gsub("^\\./|^\\./._", "", files_to_extract)) # get the names of the files without the ./ that was the prefix for untar-ing
  }
  
  
  # If there is a single file, it is a nexus with a Bayes Block, and we're done for this folder
  # If there is more than one file, then one of the files in the bayes block, and we need to append this to the nexus file
  if(length(files_to_extract)>1){
    setwd(processed_dir) # set the working directory to where the files were extracted
    nex<-files_to_extract[grep(".nex$",files_to_extract)] # get the name of the nexus file
    bayesblock<-files_to_extract[grep(".bayesblock$|.bb",files_to_extract)] # get the name of the bayes block
    nex_lines<-readLines(nex) # read in the lines of text from the nexus file
    bayesblock_lines<-readLines(bayesblock) # read in the lines of text from the bayes block
    nex_proc_name<-gsub(".nex", "_processed.nex" ,nex) # create a name for the processed Nexus file to be written
    cat(nex_lines, bayesblock_lines, file=nex_proc_name, fill=TRUE, sep="\n") # concatenate the lines of the nexus and bayes block into a new file
    unlink(c(nex, bayesblock)) # delete out the original nexus file and bayes block
    setwd(chains_dir) # set the working directory back to where it was before this if statement
  }
}

## Now we need to go into the Nexus files and reset the heating
setwd(processed_dir)

# List all of the Nexus files that we want to iterate over 
to_edit<-list.files(pattern=".nex$")
# for each of these files, change heating based on the chain acceptance for analyses that we know this for
for(i in 1:length(to_edit)){ 
  ## First identify if the analysis needs heating to be increased or decreased - if the anlaysis is not included in this object (i.e., didn't have a log file to be able to identify chain swaps), then just decrease heating--of the analyses that had this info, decreasing was more common
  analysis_base_name<-gsub("_processed|\\.nex","",to_edit[[i]])
  if(!analysis_base_name %in% names(heating_change)){
    change_direction<-"decrease"
  }else{
    change_direction<-heating_change[[analysis_base_name]]
  }
  nex_to_edit<-readLines(to_edit[[i]])# Read in the lines for the file
  if(length(grep("temp", nex_to_edit))>0){ # if the files already have a temperature specified, we need to replace the specified value
    temp_line<-nex_to_edit[grep("temp *= *([[:digit:]]+).([[:digit:]]+)", nex_to_edit)]
    old_temp<-str_extract(temp_line, "temp *= *([[:digit:]]+).([[:digit:]]+)")
    old_temp<-as.numeric(gsub("temp *= *", "", old_temp))
    if(change_direction=="increase"){    ## make the new temp either double or half of the old temp depending on if it should be increased or decreased
      new_temp<-old_temp*2
    }else{
      new_temp<-old_temp/2
    }
    nex_new_temp<-gsub("temp *= *([[:digit:]]+).([[:digit:]]+)", paste("temp=", new_temp, sep=""), nex_to_edit, ignore.case = TRUE) ### handle what to do with the files that have a temperature specified
  }else{  # otherwise, handle the nexus files that don't have a temp specified already
    old_temp<-0.1  ## if no temp specified, assume that it was the MrBayes default of 0.1
    if(change_direction=="increase"){    ## make the new temp either double or half of the old temp depending on if it should be increased or decreased
      new_temp<-old_temp*2
    }else{
      new_temp<-old_temp/2
    }
    nex_new_temp<-gsub("ngen", paste("temp=", new_temp ," ngen", sep=""), nex_to_edit, ignore.case = TRUE)  # simply find ngen (should be in all Nexus files) and replace it with a statement specifying temperature and ngen
  }
  nex_new_temp<-gsub("mcmcp", "mcmc", nex_new_temp, ignore.case = TRUE)
  nex_new_temp<-gsub("append=yes", "append=no", nex_new_temp, ignore.case = TRUE) # if any had append=yes, change that up
  setwd(mb_ready_dir) # change the directory to where we want the final MrBayes files to end up
  out_name<-gsub(".nex", "_mbReady.nex", to_edit[[i]]) # creat a name for the new Nexus file to be written
  out_folder_name<-gsub("_mbReady.nex", "", out_name)  # create new name for a folder to make and put file into - get rid of anything past _mbReady.nex
  out_folder_name<-gsub("_processed", "", out_folder_name)  # if there is a _processed, get rid of that
  dir.create(out_folder_name) # create a folder for each of these Nexus files
  setwd(out_folder_name) # move into the newly created folder
  writeLines(nex_new_temp, out_name) # write the new nexus file
  setwd(processed_dir) # set the working directory back to where the other nexus files to modify are kept
}

# Now, for each of these folders, go in and make a slurm file for running these on the cluster
setwd(mb_ready_dir) # set working directory to where the Nexus files ready for MrBayes are
folders<-list.files() # list out the folders with Nexus files in them
  
## Make a loop to generate a slurm file within each folder
for(i in 1:length(folders)){
  setwd(folders[[i]]) # move into each folder
  nexus_name<-list.files(pattern=".nex") # get the name of the nexus file
  nexlines<-readLines(nexus_name)
  nruns<-nums_from_txt(nexlines, "nruns")
  chains<-nums_from_txt(nexlines, "nchains")
  processors<-nruns*chains
  slurm<-c("#!/bin/bash", "", paste("#SBATCH -J", folders[[i]], sep=" "), paste0("#SBATCH -n ", (processors+1)),  # Set up an object with all of the slurm info
           "#SBATCH -c 1", "#SBATCH -N 1", "#SBATCH -p shared", "#SBATCH --mem=24GB",
           "#SBATCH -t 3-00:00", paste("#SBATCH -o ",folders[[i]], ".out", sep=""),
           paste("#SBATCH -e ", folders[[i]], ".err", sep=""), "source ~/.bash_profile",
           "module load bio/MrBayes/3.2.5-intel-2018.5.274-mpi", paste("mpirun -np", processors, "mb ", nexus_name, sep=" "))
  writeLines(slurm, paste(folders[[i]], ".slurm", sep="")) # write the slurm to file
  setwd(mb_ready_dir) # set the working directory back
}


#################################################################################
### Do the same for changing to use nst=mixed
setwd(processed_dir)

# List all of the Nexus files that we want to iterate over 
to_edit<-list.files(pattern=".nex$")
# for each of these files, change nst=X to nst=mixed
for(i in 1:length(to_edit)){ 
  nex_to_edit<-readLines(to_edit[[i]])# Read in the lines for the file
  nex_nst_mixed<-gsub("nst *= *([[:digit:]]+)", "nst=mixed", nex_to_edit, ignore.case = TRUE) ### change to nst=mixed
  nex_nst_mixed<-gsub("mcmcp", "mcmc", nex_nst_mixed, ignore.case = TRUE)
  nex_nst_mixed<-gsub("append=yes", "append=no", nex_nst_mixed, ignore.case = TRUE)
  setwd(mb_ready_dir_nstMixed) # change the directory to where we want the final MrBayes files to end up
  out_name<-gsub(".nex", "_mbReady.nex", to_edit[[i]]) # creat a name for the new Nexus file to be written
  out_folder_name<-gsub("_mbReady.nex", "", out_name)  # create new name for a folder to make and put file into - get rid of anything past _mbReady.nex
  out_folder_name<-gsub("_processed", "", out_folder_name)  # if there is a _processed, get rid of that
  dir.create(out_folder_name) # create a folder for each of these Nexus files
  setwd(out_folder_name) # move into the newly created folder
  writeLines(nex_nst_mixed, out_name) # write the new nexus file
  setwd(processed_dir) # set the working directory back to where the other nexus files to modify are kept
}


# Now, for each of these folders, go in and make a slurm file for running these on the cluster
setwd(mb_ready_dir_nstMixed) # set working directory to where the Nexus files ready for MrBayes are
folders<-list.files() # list out the folders with Nexus files in them

## Make a loop to generate a slurm file within each folder
for(i in 1:length(folders)){
  setwd(folders[[i]]) # move into each folder
  nexus_name<-list.files(pattern=".nex") # get the name of the nexus file
  nexlines<-readLines(nexus_name)
  nruns<-nums_from_txt(nexlines, "nruns")
  chains<-nums_from_txt(nexlines, "nchains")
  processors<-nruns*chains
  slurm<-c("#!/bin/bash", "", paste("#SBATCH -J", folders[[i]], sep=" "), paste0("#SBATCH -n ", (processors+1)),  # Set up an object with all of the slurm info
           "#SBATCH -c 1", "#SBATCH -N 1", "#SBATCH -p shared", "#SBATCH --mem=24GB",
           "#SBATCH -t 3-00:00", paste("#SBATCH -o ",folders[[i]], ".out", sep=""),
           paste("#SBATCH -e ", folders[[i]], ".err", sep=""), "source ~/.bash_profile",
           "module load bio/MrBayes/3.2.5-intel-2018.5.274-mpi", paste("mpirun -np", processors, "mb ", nexus_name, sep=" "))
  writeLines(slurm, paste(folders[[i]], ".slurm", sep="")) # write the slurm to file
  setwd(mb_ready_dir_nstMixed) # set the working directory back
}




#################################################################################
#################################################################################
### Now get the IG analyses ready to go
#################################################################################
#################################################################################
# Set working directory and list out the .tar.gz archives containing the chains
setwd(chains_dir_IG)
to_extract_IG<-list.files(pattern="*.tar.gz") # list out the .tar.gz folders that will need to be extracted

# Use a loop here to perform the extraction on each of the compressed archives in the to_extract object
for(j in 1:length(to_extract_IG)){
  # list out the files within the archive
  files_in_archive<-untar(to_extract_IG[[j]], list=TRUE)
  
  # we only want to extract the files ending in .nex, .bb, or .bayesblock
  files_to_extract<-files_in_archive[grep(pattern=".nex$|.bayesblock$|.bb", files_in_archive)]
  
  # Untar the files to the location of the processed nexus file
  if(length(grep("[[:alpha:]].*/.*[[:alpha:]]", files_to_extract)>0)){ # phylota analyses are all coming out in a directory - need to identify this, then go into the directpory and move the file to the containing folder - identify nex nested within its own directory by searching for letters separated from each other by /
    untar(to_extract_IG[[j]], files=files_to_extract, exdir=processed_dir_IG) # untar 
    files_to_extract<-unique(gsub("^\\./|^\\./._", "", files_to_extract)) # drop out ._ as for others
    file.copy(paste(processed_dir_IG, files_to_extract, sep="/"), processed_dir_IG) # copy the file from it's own directory into the directory we want as a sinle nexus file
    unlink(paste(processed_dir_IG, strsplit(files_to_extract, "/")[[1]][[1]], sep="/"), recursive = TRUE) # delete the directory - no longer necessary
    files_to_extract<-strsplit(files_to_extract, "/")[[1]][[2]] # Need to change the files_to_extract object to no longer include the containing directory
  }else{
    untar(to_extract_IG[[j]], files=files_to_extract, exdir=processed_dir_IG)
    files_to_extract<-unique(gsub("^\\./|^\\./._", "", files_to_extract)) # get the names of the files without the ./ that was the prefix for untar-ing
  }
  
  # If there is a single file, it is a nexus with a Bayes Block, and we're done for this folder
  # If there is more than one file, then one of the files in the bayes block, and we need to append this to the nexus file
  if(length(files_to_extract)>1){
    setwd(processed_dir_IG) # set the working directory to where the files were extracted
    nex<-files_to_extract[grep(".nex$",files_to_extract)] # get the name of the nexus file
    bayesblock<-files_to_extract[grep(".bayesblock$|.bb",files_to_extract)] # get the name of the bayes block
    nex_lines<-readLines(nex) # read in the lines of text from the nexus file
    bayesblock_lines<-readLines(bayesblock) # read in the lines of text from the bayes block
    nex_proc_name<-gsub(".nex", "_processed.nex" ,nex) # create a name for the processed Nexus file to be written
    cat(nex_lines, bayesblock_lines, file=nex_proc_name, fill=TRUE, sep="\n") # concatenate the lines of the nexus and bayes block into a new file
    unlink(c(nex, bayesblock)) # delete out the original nexus file and bayes block
    setwd(chains_dir_IG) # set the working directory back to where it was before this if statement
  }
}

## Now we need to go into the Nexus files and change from using I+G to only using gamma
setwd(processed_dir_IG)

# List all of the Nexus files that we want to iterate over 
to_edit<-list.files(pattern="\\.nex$")
# for each of these files, change IG to G
for(i in 1:length(to_edit)){ 
  nex_to_edit<-readLines(to_edit[[i]])# Read in the lines for the file
  nex_new_temp<-gsub("invgamma", "gamma", nex_to_edit, ignore.case = TRUE)  # replace invgamma with gammato switch from using I+G to just +G
  nex_new_temp<-gsub("mcmcp", "mcmc", nex_new_temp, ignore.case = TRUE)  # replace invgamma with gammato switch from using I+G to just +G
  nex_new_temp<-gsub("append=yes", "append=no", nex_new_temp, ignore.case = TRUE)
  setwd(mb_ready_dir_IG) # change the directory to where we want the final MrBayes files to end up
  out_name<-gsub("\\.nex", "_mbReady.nex", to_edit[[i]]) # creat a name for the new Nexus file to be written
  out_folder_name<-gsub("_mbReady.nex", "", out_name)  # create new name for a folder to make and put file into - get rid of anything past _mbReady.nex
  out_folder_name<-gsub("_processed", "", out_folder_name)  # if there is a _processed, get rid of that
  dir.create(out_folder_name) # create a folder for each of these Nexus files
  setwd(out_folder_name) # move into the newly created folder
  writeLines(nex_new_temp, out_name) # write the new nexus file
  setwd(processed_dir_IG) # set the working directory back to where the other nexus files to modify are kept
}

# Now, for each of these folders, go in and make a slurm file for running these on the cluster
setwd(mb_ready_dir_IG) # set working directory to where the Nexus files ready for MrBayes are
folders<-list.files() # list out the folders with Nexus files in them

## Make a loop to generate a slurm file within each folder
for(i in 1:length(folders)){
  setwd(folders[[i]]) # move into each folder
  nexus_name<-list.files(pattern=".nex") # get the name of the nexus file
  nexlines<-readLines(nexus_name)
  nruns<-nums_from_txt(nexlines, "nruns")
  chains<-nums_from_txt(nexlines, "nchains")
  processors<-nruns*chains
  slurm<-c("#!/bin/bash", "", paste("#SBATCH -J", folders[[i]], sep=" "), paste0("#SBATCH -n ", (processors+1)),  # Set up an object with all of the slurm info
           "#SBATCH -c 1", "#SBATCH -N 1", "#SBATCH -p shared", "#SBATCH --mem=24GB",
           "#SBATCH -t 3-00:00", paste("#SBATCH -o ",folders[[i]], ".out", sep=""),
           paste("#SBATCH -e ", folders[[i]], ".err", sep=""), "source ~/.bash_profile",
           "module load bio/MrBayes/3.2.5-intel-2018.5.274-mpi", paste("mpirun -np", processors, "mb ", nexus_name, sep=" "))
  writeLines(slurm, paste(folders[[i]], ".slurm", sep="")) # write the slurm to file
  setwd(mb_ready_dir_IG) # set the working directory back
}