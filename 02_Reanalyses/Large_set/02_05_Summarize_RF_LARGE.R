## Script to summarize the RF distances output by 02_04_RF_distances_LARGE.R


########################################################################################################################
########################################################################################################################
########################################################################################################################
##### Start with IG vs G
########################################################################################################################
########################################################################################################################


setwd("~/compare_reanalyses/RF_dists/IG")

# Read in all output - need to ensure that don't overwrite each time, as each file has same object names for summ
to_read<-list.files(pattern="rf_comp_IG")
out_all_pre<-list()    #make a list to dump the batches into
for(i in 1:length(to_read)){   # Dump the summ object from each batch into an element of a new list
  load(to_read[[i]])
  out_all_pre[[i]]<-rf_comparisons_IG
}

out_all<-unlist(out_all_pre, recursive=FALSE)

## Get out the RF scores only and plot them out
IGvsG_RF_only<-do.call(rbind, lapply(out_all, function(x) x[[1]]))

# Get the median to add to the title of each
med_RF_95<-round(median(IGvsG_RF_only[,1], na.rm=TRUE), 3)
med_RF_50<-round(median(IGvsG_RF_only[,2], na.rm=TRUE), 3)
med_RF_MCC<-round(median(IGvsG_RF_only[,3], na.rm=TRUE), 3)

### Plots of the change between new and old heating for 95% and 50% consensus trees and Maximum Clade Credibility (MCC) trees
setwd("~/compare_reanalyses")
pdf(file="RF IG_G.pdf", height=4, width=11)
par(mfrow=c(1,3))
# For chains that had bad Geweke's
hist(IGvsG_RF_only[,1], col="gray", main=paste("IG vs G 95% con \nMedian = ",med_RF_95,sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL)
hist(IGvsG_RF_only[,2], col="gray", main=paste("IG vs G 50% con \nMedian = ",med_RF_50,sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL)
hist(IGvsG_RF_only[,3], col="gray", main=paste("IG vs G MCC \nMedian = ",med_RF_MCC,sep=""), xlab="Norm RF", xlim=c(0,1), breaks = 20, ylab=NULL)
dev.off()# turn off the pdf plotting device


########################################################################################################################
########################################################################################################################
########################################################################################################################
##### Move on to Nst = mixed
########################################################################################################################
########################################################################################################################
setwd("~/compare_reanalyses/RF_dists/nstMixed")

# Read in all output - need to ensure that don't overwrite each time, as each file has same object names for summ
to_read_NstMixed<-list.files(pattern="rf_comp_NstMixed")
out_all_NstMixed_pre<-list()    #make a list to dump the batches into
for(i in 1:length(to_read_NstMixed)){   # Dump the summ object from each batch into an element of a new list
  load(to_read_NstMixed[[i]])
  out_all_NstMixed_pre[[i]]<-rf_comparisons_NstMixed
}

out_all_NstMixed<-unlist(out_all_NstMixed_pre, recursive=FALSE)

## Get out the RF scores only and plot them out
NstMixed_RF_only<-do.call(rbind, lapply(out_all_NstMixed, function(x) x[[1]]))

# Get the median to add to the title of each
med_RF_95_NstMixed<-round(median(NstMixed_RF_only[,1], na.rm=TRUE), 3)
med_RF_50_NstMixed<-round(median(NstMixed_RF_only[,2], na.rm=TRUE), 3)
med_RF_MCC_NstMixed<-round(median(NstMixed_RF_only[,3], na.rm=TRUE), 3)

### Plots of the change between new and old heating for 95% and 50% consensus trees and Maximum Clade Credibility (MCC) trees
setwd("~/compare_reanalyses")
pdf(file="RF_NstMixed.pdf", height=4, width=11)
par(mfrow=c(1,3))
# For chains that had bad Geweke's
hist(NstMixed_RF_only[,1], col="gray", main=paste("Nst Mixed 95% con \nMedian = ",med_RF_95_NstMixed,sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL)
hist(NstMixed_RF_only[,2], col="gray", main=paste("Nst Mixed 50% con \nMedian = ",med_RF_50_NstMixed,sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL)
hist(NstMixed_RF_only[,3], col="gray", main=paste("Nst Mixed MCC \nMedian = ",med_RF_MCC_NstMixed,sep=""), xlab="Norm RF", xlim=c(0,1), breaks = 20, ylab=NULL)
dev.off()# turn off the pdf plotting device



########################################################################################################################
########################################################################################################################
########################################################################################################################
##### Move on to new heating
########################################################################################################################
########################################################################################################################


setwd("~/compare_reanalyses/RF_dists/NewHeat")

# Read in all output - need to ensure that don't overwrite each time, as each file has same object names for summ
to_read_New_heat<-list.files(pattern="rf_comp_New_heat")
out_all_New_heat_pre<-list()    #make a list to dump the batches into
for(i in 1:length(to_read_New_heat)){   # Dump the summ object from each batch into an element of a new list
  load(to_read_New_heat[[i]])
  out_all_New_heat_pre[[i]]<-rf_comparisons_New_heat
}

out_all_New_heat<-unlist(out_all_New_heat_pre, recursive=FALSE)

## Get out the RF scores only and plot them out
New_heat_RF_only<-do.call(rbind, lapply(out_all_New_heat, function(x) x[[1]]))

# Get the median to add to the title of each
med_RF_95_New_heat<-round(median(New_heat_RF_only[,1], na.rm=TRUE), 3)
med_RF_50_New_heat<-round(median(New_heat_RF_only[,2], na.rm=TRUE), 3)
med_RF_MCC_New_heat<-round(median(New_heat_RF_only[,3], na.rm=TRUE), 3)

### Plots of the change between new and old heating for 95% and 50% consensus trees and Maximum Clade Credibility (MCC) trees
setwd("~/compare_reanalyses")
pdf(file="RF_New_heat.pdf", height=4, width=11)
par(mfrow=c(1,3))
# For chains that had bad Geweke's
hist(New_heat_RF_only[,1], col="gray", main=paste("New Heat 95% con \nMedian = ",med_RF_95_New_heat,sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL)
hist(New_heat_RF_only[,2], col="gray", main=paste("New Heat 50% con \nMedian = ",med_RF_50_New_heat,sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL)
hist(New_heat_RF_only[,3], col="gray", main=paste("New Heat MCC \nMedian = ",med_RF_MCC_New_heat,sep=""), xlab="Norm RF", xlim=c(0,1), breaks = 20, ylab=NULL)
dev.off()# turn off the pdf plotting device



### Plot all of the histograms together into 1 figure

### Plots of the change between new and old heating for 95% and 50% consensus trees and Maximum Clade Credibility (MCC) trees
setwd("~/compare_reanalyses")
pdf(file="RF_dists_Large.pdf", height=12, width=11)
par(mfrow=c(3,3))
# For chains that had bad Geweke's
hist(IGvsG_RF_only[,1], col="gray", main=paste("IG vs G 95% con \nMedian = ",med_RF_95,sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL)
hist(IGvsG_RF_only[,2], col="gray", main=paste("IG vs G 50% con \nMedian = ",med_RF_50,sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL)
hist(IGvsG_RF_only[,3], col="gray", main=paste("IG vs G MCC \nMedian = ",med_RF_MCC,sep=""), xlab="Norm RF", xlim=c(0,1), breaks = 20, ylab=NULL)

hist(NstMixed_RF_only[,1], col="gray", main=paste("Nst Mixed 95% con \nMedian = ",med_RF_95_NstMixed,sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL)
hist(NstMixed_RF_only[,2], col="gray", main=paste("Nst Mixed 50% con \nMedian = ",med_RF_50_NstMixed,sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL)
hist(NstMixed_RF_only[,3], col="gray", main=paste("Nst Mixed MCC \nMedian = ",med_RF_MCC_NstMixed,sep=""), xlab="Norm RF", xlim=c(0,1), breaks = 20, ylab=NULL)

hist(New_heat_RF_only[,1], col="gray", main=paste("New Heat 95% con \nMedian = ",med_RF_95_New_heat,sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL)
hist(New_heat_RF_only[,2], col="gray", main=paste("New Heat 50% con \nMedian = ",med_RF_50_New_heat,sep=""), xlab="Norm RF", xlim=c(0,1), ylab=NULL)
hist(New_heat_RF_only[,3], col="gray", main=paste("New Heat MCC \nMedian = ",med_RF_MCC_New_heat,sep=""), xlab="Norm RF", xlim=c(0,1), breaks = 20, ylab=NULL)
dev.off()# turn off the pdf plotting device

