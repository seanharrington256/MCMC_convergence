## Script to get out all of the original folders for the chains to be reanalyzed from their storage locations
  # this is specific to where these files are stored on the UH cluster for our analyses, change all files 
  # based on their location on your system

# Set up the locations where the data folders of the original analyses are stored
storage1<-"Sean_storage/All_Barcoding_w_BayesBlocks/Data"
storage2<-"Sean_storage/mtDNA_cleaned/Data"
storage3<-"Sean_storage/All_turtles/Data"
storage4<-"Van_storage/Phylota_v2.1_complete_all"

# Set up the location that we will mainly work from
main_dir<-"Sean_storage/MCMC_reanalyses"
# Set up the location that we will want to copy things over into
copy_to<-paste0(main_dir, "/all_original_analyses")
copy_to_IG<-paste0(main_dir, "/IG_original_analyses")
copy_to_manual<-paste0(main_dir, "/Manual_original_analyses")

## Read in the names of the analyses we want to reanalyze
setwd(main_dir)
set_1<-read.csv("fail_ASDSF_LARGE.csv", stringsAsFactors = FALSE)
set_2<-read.csv("fail_oneormore_ESS_LARGE.csv", stringsAsFactors = FALSE)
set_3<-read.csv("fail_topo_bad_LnL_ESS_LARGE.csv", stringsAsFactors = FALSE)
set_4<-read.csv("fail_topo_good_LnL_ESS_LARGE.csv", stringsAsFactors = FALSE)
all_to_reestimate<-c(set_1[,2], set_2[,2], set_3[,2], set_4[,2]) # combine these into a single object for ease

## Read in the analyses to rerun without +I
IG_get_1<-read.csv("high_IG_cor_bad_ESS_LARGE.csv", stringsAsFactors = FALSE)
IG_get_2<-read.csv("high_IG_cor_good_conv_LARGE.csv", stringsAsFactors = FALSE)
IG_to_reestimate<-c(IG_get_1[,2], IG_get_2[,2]) # combine these into a single object for ease

## Read in analyses to figure out manually
manual_to_reestimate<-read.csv("Manually_reanalyze.csv", stringsAsFactors = FALSE)[,2]



## For each of the storage locations,
# go to the storage location of the folders and copy over the folders we need
all_storage<-c(storage1, storage2, storage3, storage4)
for(j in all_storage){
  setwd(j)
  # List out all of the .tar.gz folders
  all_folders<-list.files(pattern=".tar.gz")
  base_names_all_folders<-cbind(all_folders, gsub("_READY|.tar.gz", "", all_folders)) #strip off extensions
  # Pull out the analyses we need and copy them over
  index<-which(base_names_all_folders[,2] %in% all_to_reestimate)
  if(length(index)>0){
    folders_to_copy<-base_names_all_folders[index,1]
    # copy over the folders
    for(i in 1:length(folders_to_copy)){
      folder_copy<-paste(j, folders_to_copy[[i]], sep="/")
      file.copy(folder_copy, copy_to, recursive=TRUE)
    }
  }
  ### Same thing for +IG
  # Pull out the analyses
  index_IG<-which(base_names_all_folders[,2] %in% IG_to_reestimate)
  if(length(index_IG)>0){
    folders_to_copy_IG<-base_names_all_folders[index_IG,1]
    # copy over the folders
    for(i in 1:length(folders_to_copy_IG)){
      folder_copy_IG<-paste(j, folders_to_copy_IG[[i]], sep="/")
      file.copy(folder_copy_IG, copy_to_IG, recursive=TRUE)
    }
  }
  ### Same thing for manual adjustment
  # Pull out the analyses
  index_manual<-which(base_names_all_folders[,2] %in% manual_to_reestimate)
  if(length(index_manual)>0){
    folders_to_copy_manual<-base_names_all_folders[index_manual,1]
    # copy over the folders
    for(i in 1:length(folders_to_copy_manual)){
      folder_copy_manual<-paste(j, folders_to_copy_manual[[i]], sep="/")
      file.copy(folder_copy_manual, copy_to_manual, recursive=TRUE)
    }
  }
}

