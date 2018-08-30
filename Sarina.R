# your floder where all pages are saved
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")


# DataSet <-read.fasta("DengueVirus1.fasta_pruned.mu.trim05.txt")
source("MeaFreq.R")
# must place your file as a txt takes a few minutes 
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/virus")
virusname = 'BKpolyomavirus_VP1.fasta.mu'
DF<-meanFreq(virusname)
splitname<-unlist(strsplit(virusname,".fasta"))
truename<-splitname[1]
DF$wtnt<-as.character(DF$wtnt)
DF$Virus<-(truename)


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
virusname = 'BKpolyomavirus_VP1.fasta.mu'
splitname<-unlist(strsplit(virusname,".fasta"))
truename<-splitname[1]

truenameRda= paste(truename, ".Rda", sep="")
load(truenameRda)

# graphs 

truenamepdf = paste(truename,".pdf",sep="")
pdf(truenamepdf, width = 6.75, height = 6.75)
source("Synonymous and nonsynonymous graph.R")
dev.off()
SynNonsynAT(DF, truename)
SynNonsynCG(DF, truename)

truenamepdf = paste(truename,"CpG.pdf",sep="")
pdf(truenamepdf, width = 6.75, height = 3.75)
source("MakesCPG Graph.R")
dev.off()
CPGNoCPGAT(DF, truename)

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




# graphs 
truenamepng = paste(truename,".df",sep="")
pdf(truenamepng, width = 6.75, height = 6.75)
source("Synonymous and nonsynonymous graph.R")
dev.off()
SynNonsynAT(DF,truename)
SynNonsynCG(DF,truename)

truenamepdf = paste(truename,".pdf",sep="")
pdf(truenamepdf, width = 6.75, height = 3.75)
source("MakesCPG Graph.R")
dev.off()
CPGNoCPGAT(DF)

source("CPG_Syn_Nonsyn_graph.R")
pdf('comparing_CpG_Syn_EnterovirusA_VP.pdf',width = 7, height = 5)
comparing_CpG_Syn_Nonsyn (DF)
dev.off()

# Wilcox test
source("stats.R")
Wilcox_test(DF)

#########
########
source("MeaFreq.R")
# must place your file as a txt takes a few minutes 
DF<-meanFreq('EnterovirusB_VP1.fasta_pruned.mu.trim05.txt')
DF$wtnt<-as.character(DF$wtnt)
start = 142
end= 0
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

# graphs 
pdf('Synon_Nonsyn_EnteroB_Graph.pdf',width = 7, height = 5)
source("Synonymous and nonsynonymous graph.R")
dev.off()
SynNonsynAT(DF)
SynNonsynCG(DF)

pdf('CPG_NonCPG_AT_EnteroB_Graph.pdf',width = 7, height = 5)
source("MakesCPG Graph.R")
dev.off()
CPGNoCPGAT(DF)

source("CPG_Syn_Nonsyn_graph.R")
pdf('comparing_CpG_Syn_EnteroB_Graph.pdf',width = 7, height = 5)
comparing_CpG_Syn_Nonsyn (DF)
dev.off()

# Wilcox test
source("stats.R")
Wilcox_test(DF)
################################################
source("MeaFreq.R")
# must place your file as a txt takes a few minutes 
DF<-meanFreq('EnterovirusB_VP2.fasta.mu.trim08')
DF$wtnt<-as.character(DF$wtnt)
start = 21
end= 0
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

# graphs 
pdf('Synon_Nonsyn_EnteroA_VP2.pdf',width = 7, height = 5)
source("Synonymous and nonsynonymous graph.R")
dev.off()
SynNonsynAT(DF)
SynNonsynCG(DF)

source("CPG_Syn_Nonsyn_graph.R")
pdf('comparing_CpG_Syn_EnteroA_VP2.pdf',width = 7, height = 5)
comparing_CpG_Syn_Nonsyn (DF)
dev.off()

# Wilcox test
source("stats.R")
Wilcox_test(DF)

##########################################
source("MeaFreq.R")
# must place your file as a txt takes a few minutes 
DF<-meanFreq('HumanBocavirus1_VP1.fasta_pruned.mu.trim05')
DF$wtnt<-as.character(DF$wtnt)
start = 1
end= 0
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

# graphs 
pdf('Synon_Nonsyn_Bova_VP1.pdf',width = 7, height = 5)
source("Synonymous and nonsynonymous graph.R")
dev.off()
SynNonsynAT(DF)
SynNonsynCG(DF)

source("CPG_Syn_Nonsyn_graph.R")
pdf('comparing_CpG_Syn_Boca_VP1.pdf',width = 7, height = 5)
comparing_CpG_Syn_Nonsyn (DF)
dev.off()

# Wilcox test
source("stats.R")
Wilcox_test(DF)

###########################
source("MeaFreq.R")
# must place your file as a txt takes a few minutes 
DF<-meanFreq('EnterovirusB_VP2.fasta.mu.trim08')
DF$wtnt<-as.character(DF$wtnt)
start = 42
end= 14
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

# graphs 
pdf('Synon_Nonsyn_InfB_HA.pdf',width = 7, height = 5)
source("Synonymous and nonsynonymous graph.R")
dev.off()
SynNonsynAT(DF)
SynNonsynCG(DF)

source("CPG_Syn_Nonsyn_graph.R")
pdf('comparing_CpG_Syn_Inf_HA.pdf',width = 7, height = 5)
comparing_CpG_Syn_Nonsyn (DF)
dev.off()

# Wilcox test
source("stats.R")
Wilcox_test(DF)



