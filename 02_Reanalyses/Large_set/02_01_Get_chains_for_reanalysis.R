## Script to get out all of the original folders for the chains to be reanalyzed

# Set up the location where the data folders of the original analyses are stored
storage<-"~/apps/1_All_orig_chains"
# Set up the location that we' will want to copy things over into'll mainly work from
main_dir<-getwd()
# Set up the location that we will want to copy things over into
copy_to<-paste(getwd(), "/all_original_analyses", sep="")
copy_to_IG<-paste(getwd(), "/IG_original_analyses", sep="")

# go to the storage location of the folders:
setwd(storage)

# List out all of the .tar.gz folders
all_folders<-list.files(pattern=".tar.gz")
base_names_all_folders<-cbind(all_folders, gsub("_READY|.tar.gz", "", all_folders)) #strip off extensions

## Read in the names of the analyses we want to reanalyze
setwd(main_dir)
set_1<-read.csv("fail_Geweke_pass_common_LARGE.csv", stringsAsFactors = FALSE)
set_2<-read.csv("fail_oneormore_ESS_LARGE.csv", stringsAsFactors = FALSE)
set_3<-read.csv("fail_topo_bad_LnL_ESS_LARGE.csv", stringsAsFactors = FALSE)
set_4<-read.csv("fail_topo_good_LnL_ESS_LARGE.csv", stringsAsFactors = FALSE)
all_to_reestimate<-c(set_1[,2], set_2[,2], set_3[,2], set_4[,2]) # combine these into a single object for ease
# Pull out the analyses that match these
    index<-which(base_names_all_folders[,2] %in% all_to_reestimate)
    folders_to_copy<-base_names_all_folders[index,1]
    # Create a directory to drop the folders into
    dir.create("all_original_analyses")
    # copy over the folders
    for(i in 1:length(folders_to_copy)){
      folder_copy<-paste(storage, folders_to_copy[[i]], sep="/")
      file.copy(folder_copy, "all_original_analyses", recursive=TRUE)
    }

#######################################################################################################
#######################################################################################################
### Do the same as above for the analyses that we want to analyze with +G only instead of +IG
#######################################################################################################
#######################################################################################################
    
## Read in the analyses we want to rerun without +I
    setwd(main_dir)
    IG_get_1<-read.csv("high_IG_cor_bad_ESS_LARGE.csv", stringsAsFactors = FALSE)
    IG_get_2<-read.csv("high_IG_cor_good_conv_LARGE.csv", stringsAsFactors = FALSE)
    IG_to_reestimate<-c(IG_get_1[,2], IG_get_2[,2]) # combine these into a single object for ease
    # Pull out the analyses that match these
    index_IG<-which(base_names_all_folders[,2] %in% IG_to_reestimate)
    folders_to_copy_IG<-base_names_all_folders[index_IG,1]
    # Create a directory to drop the folders into
    dir.create("IG_original_analyses")
    # copy over the folders
    for(i in 1:length(folders_to_copy_IG)){
      folder_copy_IG<-paste(storage, folders_to_copy_IG[[i]], sep="/")
      file.copy(folder_copy_IG, "IG_original_analyses", recursive=TRUE)
    }
    
