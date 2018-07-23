# your floder where all pages are saved
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")
#This is by virus must change for new virus
#start = 1
#end = 0 
#Virus = "DengueVirus1.fasta_pruned.mu.trim05.txt"

#start= 75 #D2
#end = 327 
#Virus= "DengueVirus2.fasta_pruned.mu.trim05.txt"

#start= 71 #D3
#end = 414 
#Virus= "DengueVirus3.fasta_pruned.mu.trim05.txt"

#start= 81 #D4
#end = 363 
#Virus= "DengueVirus4.fasta_pruned.mu.trim05.txt"

#start= 187 
# Virus= "EnterovirusA_VP1.fasta_pruned.mu.trim05"
#end = 0

start = 0
end = 0
#Virus = "InfluenzaBvirus_HA.fasta.mu.trim05"
start = 178
end = 3060
Virus = "HumanBocavirus1_NS1.fasta_pruned.mu.trim05"
start = 2981
end = 159
Virus = "HumanBocavirus1_NS1.fasta_pruned.mu.trim05"
# DataSet <-read.fasta("DengueVirus1.fasta_pruned.mu.trim05.txt")
source("startstop.R") # either of these works I keeping this one becuase someone might have changed the other
#source("MeaFreq.R") # everything is good
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/virus")
# File can be a txt or fasta file and takes a few minutes 
# this created the inital data frame where there are wildtype nucleotides/ frequency/ nucelotide number
DF<-meanFreq(Virus)
#this will shorten the data fram to the gene Must use traslate tool for start/end
DF<-DF[c(start:(nrow(DF)-end)), ]
DF$num<-(1:nrow(DF))
DF$wtnt<-as.character(DF$wtnt)

setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")
#finds the wildtype amino acid
source("WTAA.R")
DF<-getWTAA(DF)

#finds the mutated amino acid
source("MUTAA.R")
DF<-getMUTAA(DF)

#see if there is a drastic Amino Acid change. like from positive to negative
source("Drastic_AA_Change.R")
DF<-big_aa_change(DF)

#determine if AA change is synonymous or nonsynonymous
source("SynNonSyn.R")
DF<-synFunction(DF)

#Marks if there is a new CpG mutation
source("CPG_Function.R")
DF<-CPG_site(DF)

# graphs made into pdfs
pdf('Synon_Nonsyn_AT_Graph_HumanBocavirus1_VP1.pdf',width = 7, height = 5)
source("Synonymous and nonsynonymous graph.R")
dev.off()
#SynNonsynAT(DF)
#SynNonsynCG(DF)

pdf('Synon_Nonsyn_CG_HumanBocavirus1_VP1.pdf',width = 7, height = 5)
#source("MakesCPG Graph.R")
#CPGNoCPGAT(DF)
# pdf('Compare_CpG.pdf',width = 7, height = 5)
source("CPG_Syn_Nonsyn_graph.R")
comparing_CpG_Syn_Nonsyn (DF)
dev.off()

# Wilcox test
source("stats.R")
Wilcox_test(DF)



