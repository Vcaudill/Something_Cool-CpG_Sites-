# your floder where all pages are saved




setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-")
source("MeaFreq.R")
library(knitr)
#makeTable <- function(wilcox_test){
 # kable(wilcox_test)
  #table <- matrix(c(wilcox_test),ncol=3,byrow=TRUE)
  #colnames(smoke) <- c("A-G","T-C")
  #rownames(smoke) <- c("Syn v. Nonsyn","SynCpG v. NonCpG", "NonSynCpG v. NonCpG")
#}
# must place your file as a txt takes a few minutes 
DF<-meanFreq('EnterovirusD_VP1.fasta_pruned.mu.trim05.txt')
DF$wtnt<-as.character(DF$wtnt)
start = 86
end = 979


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
pdf("VP1EnteroGraphsns.pdf", width = 7, height= 5)
source("Synonymous and nonsynonymous graph.R")
dev.off()
SynNonsynAT(DF)
SynNonsynCG(DF)



#source("MakesCPG Graph.R")
#CPGNoCPGAT(DF)

source("CPG_Syn_Nonsyn_graph.R")
pdf("VP1EnteroGraph2.pdf", width = 7, height= 5)
comparing_CpG_Syn_Nonsyn (DF)
#dev.off()

b # Wilcox test
source("RyanWilcox.R")
Wilcox_test(DF)
dev.off()
########
########


