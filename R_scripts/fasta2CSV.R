# your floder where all pages are saved
#setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")
setwd("~ryanwinstead/Something_Cool-CpG_Sites-")

#loop using a csv to find the file name, strating and stoping points, and readign frame

Virus_info<- read.csv("new_data/CpG_list_Newdata.csv")
#RW changed file from MEME_CpG_list.csv to Newdata.csv
i=1
# DataSet <-read.fasta("DengueVirus1.fasta_pruned.mu.trim05.txt")
for(i in 1:nrow(Virus_info)){
  print(i, Virus_info$name[i])
  setwd("~ryanwinstead/Something_Cool-CpG_Sites-")
  source("R_scripts/BioInfo/Freq.R")
  # must place your file as a txt takes a few minutes 
  #setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/virus")
  fasta_file = paste('new_data/good_to_go/', Virus_info$name[i], ".fasta",sep="")#trimmed
  DF<-Freq(fasta_file)
  nucleotideNumber <- nrow(DF)
  print(nucleotideNumber)
  name = as.character(Virus_info$name[i])
  splitname<-unlist(strsplit(as.character(Virus_info$name[i]),".fasta"))
  truename<-splitname[1]
  DF$wtnt_consensus<-as.character(DF$wtnt_consensus)
  DF$Virus<-(truename)
  #setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")
  
  source("R_scripts/BioInfo/WTAA_consensus.R")
  DF<-getWTAA(DF)
  
  source("R_scripts/BioInfo/MUTAA.R")
  DF<-getMUTAA(DF)
  
  source("R_scripts/BioInfo/Drastic_AA_Change.R")
  DF<-big_aa_change(DF)
  
  source("R_scripts/BioInfo/CPG_Function.R")
  DF<-CPG_site(DF)
  
  source("R_scripts/BioInfo/SynNonSyn.R")
  DF<-synFunction(DF)
  
  #How to save data
  truenameCSV= paste('new_data/Consensus/', truename, ".csv", sep="")
  write.csv(DF, file = truenameCSV)

  #add in the graphs/tables
  }

source("R_scripts/BioInfo/SynNonSyn.R")
DF<-synFunction(DF)

source("R_scripts/BioInfo/CPG_Function.R")
DF<-CPG_site(DF)
#How to load data
virusname = 'EnterovirusB_VP2.fasta_pruned.mu.trim05.txt'
splitname<-unlist(strsplit(virusname,".fasta"))
truename<-splitname[1]

truenameRda= paste('new_data/Rda/', truename, ".Rda", sep="")
truenameCSV= paste('new_data/Csv/', truename, ".csv", sep="")
write.csv(DF, file = "MyData.csv")

#setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/Rda_Files")
# load(truenameRda)
###############################################
# graphs/ tables 
#setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")
# source("R_scripts/BioInfo/Synonymous and nonsynonymous graph.R")
# SynNonsynAT(DF)
# SynNonsynCG(DF)
# 
# 
# source("R_scripts/BioInfo/MakesCPG Graph.R")
# CPGNoCPGAT(DF)
# 
# source("R_scripts/BioInfo/CPG_Syn_Nonsyn_graph.R")
# comparing_CpG_Syn_Nonsyn (DF)

# Wilcox test
#source("stats.R")
#Wilcox_test(DF)
# that tables
#source("R_scripts/BioInfo/RyanWilcox.R")
#Wilcox_test(DF, truename)

#source("R_scripts/graphs/redoplot.R")
#comparing_CpG_Syn_Nonsyn_new(DF, truename)

#########
########


