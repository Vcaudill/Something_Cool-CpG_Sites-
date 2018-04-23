# your floder where all pages are saved
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")


# DataSet <-read.fasta("DengueVirus1.fasta_pruned.mu.trim05.txt")
source("MeaFreq.R")
# must place your file as a txt takes a few minutes 
DF<-meanFreq('EnterovirusA_VP1.fasta_pruned.mu.trim05.txt')
DF$wtnt<-as.character(DF$wtnt)

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
pdf('Synon_Nonsyn_Graph.pdf',width = 7, height = 5)
source("Synonymous and nonsynonymous graph.R")
dev.off()
SynNonsynAT(DF)
SynNonsynCG(DF)

pdf('CPG_NonCPG_AT_Graph.pdf',width = 7, height = 5)
source("MakesCPG Graph.R")
dev.off()
CPGNoCPGAT(DF)

source("CPG_Syn_Nonsyn_graph.R")
pdf('comparing_CpG_Syn_Graph.pdf',width = 7, height = 5)
comparing_CpG_Syn_Nonsyn (DF)
dev.off()

# Wilcox test
source("stats.R")
Wilcox_test(DF)

#########
########


