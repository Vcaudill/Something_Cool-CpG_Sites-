# your floder where all pages are saved
setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-")

#How to load data
virusname = 'EnterovirusB_VP2.fasta_pruned.mu.trim05.txt'
splitname<-unlist(strsplit(virusname,".fasta"))
truename<-splitname[1]

truenameRda= paste(truename, ".Rda", sep="")
setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-/Rda_Files")
load(truenameRda)


setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-")
# Wilcox test
source("RyanWilcox.R")
Wilcox_test(DF, truename)

########
########


