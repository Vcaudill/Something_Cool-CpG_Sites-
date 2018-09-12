# your floder where all pages are saved
#setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")


# DataSet <-read.fasta("DengueVirus1.fasta_pruned.mu.trim05.txt")
source("R_scripts/BioInfo/Freq.R")
# must place your file as a txt takes a few minutes 
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/virus")
virusname = 'EnterovirusB_VP2.fasta.mu.trim08'
virusname = 'DengueVirus1.fasta_pruned.mu.trim05.txt'
DF<-meanFreq(virusname)
splitname<-unlist(strsplit(virusname,".fasta"))
truename<-splitname[1]
DF$wtnt<-as.character(DF$wtnt)
DF$Virus<-(truename)
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")

# for Ryan 
# make sure that M is the first amino acid check all virus
start = 1
end = 0
DF<-DF[c(start:(nrow(DF)-end)),]
DF$num<-(1:nrow(DF))


source("WTAA.R")
DF<-getWTAA(DF)

source("MUTAA.R")
DF<-getMUTAA(DF)

source("Drastic_AA_Change.R")
DF<-big_aa_change(DF)

source("SynNonSyn.R")
DF<-synFunction(DF)

source("CPG_Function.R")
DF<-CPG_site(DF)
#How to save data
truenameRda= paste(truename, ".Rda", sep="")
save("virus" = DF,file=truenameRda)


#How to load data
virusname = 'EnterovirusB_VP2.fasta_pruned.mu.trim05.txt'
splitname<-unlist(strsplit(virusname,".fasta"))
truename<-splitname[1]

truenameRda= paste(truename, ".Rda", sep="")
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/Rda_Files")
load(truenameRda)

# graphs/ tables 
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")
source("Synonymous and nonsynonymous graph.R")
SynNonsynAT(DF)
SynNonsynCG(DF)


source("MakesCPG Graph.R")
CPGNoCPGAT(DF)

source("CPG_Syn_Nonsyn_graph.R")
comparing_CpG_Syn_Nonsyn (DF)

# Wilcox test
#source("stats.R")
#Wilcox_test(DF)

source("RyanWilcox.R")
Wilcox_test(DF, truename)

source("redoplot.R")
comparing_CpG_Syn_Nonsyn_new(DF, truename)

#########
########


