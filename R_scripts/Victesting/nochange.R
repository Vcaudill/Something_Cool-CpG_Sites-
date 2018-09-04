# your floder where all pages are saved
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")

# DataSet <-read.fasta("DengueVirus1.fasta_pruned.mu.trim05.txt")
source("MeaFreq.R")
# must place your file as a txt takes a few minutes 
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/virus")
DF<-meanFreq("EnterovirusB_VP2.fasta.mu.trim08")
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")
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
pdf('~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/PDF/Synon_Nonsyn_AT_Graph_EnterovirusB_VP2.pdf',width = 7, height = 5)
source("Synonymous and nonsynonymous graph.R")
dev.off()
#SynNonsynAT(DF)
#SynNonsynCG(DF)

pdf('~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/PDF/Synon_Nonsyn_CG_Graph_EnterovirusB_VP2.pdf',width = 7, height = 5)
#source("MakesCPG Graph.R")
#CPGNoCPGAT(DF)
# pdf('Compare_CpG.pdf',width = 7, height = 5)
source("CPG_Syn_Nonsyn_graph.R")
comparing_CpG_Syn_Nonsyn (DF)
dev.off()

# Wilcox test
source("stats.R")
Wilcox_test(DF)
# pdf('Synon_Nonsyn_AT_Graph.pdf',width = 7, height = 5)
# SynNonsynAT(DF)
# # dev.off()
# pdf('Synon_Nonsyn_CG_Graph.pdf',width = 7, height = 5)
# SynNonsynCG(DF)
# dev.off()
#########
########

print(toString(DF$wtnt))


