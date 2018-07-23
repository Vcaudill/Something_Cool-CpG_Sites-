# your floder where all pages are saved
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")


# DataSet <-read.fasta("DengueVirus1.fasta_pruned.mu.trim05.txt")
source("MeaFreq.R")
# must place your file as a txt takes a few minutes 
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/virus")
virusname = 'HumanBocavirus1_NS1.fasta_pruned.mu.trim05'
virusname = 'DengueVirus1.fasta_pruned.mu.trim05.txt'
DF<-meanFreq(virusname)
splitname<-unlist(strsplit(virusname,".fasta"))
truename<-splitname[1]
DF$wtnt<-as.character(DF$wtnt)
DF$Virus<-(truename)
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")

# for Ryan 
# make sure that M is the first amino acid check all virus
start = 178
end = 3060
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
virusname = 'DengueVirus1.fasta_pruned.mu.trim05.txt'
splitname<-unlist(strsplit(virusname,".fasta"))
truename<-splitname[1]

truenameRda= paste(truename, ".Rda", sep="")
load(truenameRda)

# graphs 

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


