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
virusname = 'HumanBocavirus1_VP1.fasta_pruned.mu.trim05'
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

source("R_scripts/graphs/redoplot.R")
comparing_CpG_Syn_Nonsyn_new(truename)

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


ListofFastaFiles <-list.files("data/fasta/", pattern = c("trim05"), all.files = TRUE)
#my.list <- list('DengueVirus1', 'DengueVirus2', 'DengueVirus3', 'DengueVirus4', 'humanparainfluenzavirus1_F', 'humanparainfluenzavirus1_HN', 'humanparainfluenzavirus3_HN', 'InfluenzaAvirus_HA_H1N1','InfluenzaAvirus_HA_H3N2', 'InfluenzaAvirus_NA_H1N1', 'InfluenzaAvirus_NA_H3N2','InfluenzaBvirus_HA', 'InfluenzaBvirus_NA', 'EnterovirusA_VP1', 'EnterovirusA_VP2','EnterovirusB_VP1', 'EnterovirusB_VP2','EnterovirusC_VP1','EnterovirusC_VP2','EnterovirusD_VP1', 'BKpolyomavirus_VP1', 'HumanBocavirus1_NS1', 'HumanBocavirus1_VP1')
fasta<-c("BKpolyomavirus_VP1.fasta.mu.trim05","DengueVirus1.fasta_pruned.mu.trim05","DengueVirus2.fasta_pruned.mu.trim05","DengueVirus3.fasta_pruned.mu.trim05","DengueVirus4.fasta_pruned.mu.trim05","EnterovirusA_VP1.fasta_pruned.mu.trim05","EnterovirusA_VP2.fasta.mu.trim08","EnterovirusB_VP1.fasta_pruned.mu.trim05"
      ,"EnterovirusB_VP2.fasta.mu.trim08","EnterovirusC_VP1.fasta_pruned.mu.trim05","EnterovirusC_VP2.fasta_pruned.mu.trim05","EnterovirusD_VP1.fasta_pruned.mu.trim05","HepatitisB.fasta_pruned.mu.trim05"
      ,"HumanBocavirus1_NS1.fasta_pruned.mu.trim05","HumanBocavirus1_VP1.fasta_pruned.mu.trim05","Humanherpesvirus2_gD.fasta_pruned.mu.trim05","Humanherpesvirus2_gD.fasta_pruned.mu.trim05","Humanherpesvirus2_glycoprotein_G.fasta_pruned.mu.trim05"
      ,"Humanpapillomavirus16.fasta_pruned.mu.trim05","Humanpapillomavirus16_L1.fasta_pruned.mu.trim05","humanparainfluenzavirus1.fasta_pruned.mu.trim05","humanparainfluenzavirus1_F.fasta_pruned.mu.trim05","humanparainfluenzavirus1_HN.fasta_pruned.mu.trim05"
      ,"humanparainfluenzavirus3.fasta_pruned.mu.trim05","humanparainfluenzavirus3_HN.fasta_pruned.mu.trim05","Humanrespiratorysyncytialvirus.fasta.mu.trim05","Humanrespiratorysyncytialvirus_G.fasta.mu.trim05","Humanrespiratorysyncytialvirus_G.fasta_pruned.mu.trim05"
      ,"InfluenzaAvirus_HA_H1N1.fasta.mu.trim05","InfluenzaAvirus_HA_H3N2.fasta.mu.trim05","InfluenzaAvirus_NA_H1N1.fasta.mu.trim05","InfluenzaAvirus_NA_H3N2.fasta.mu.trim05","InfluenzaBvirus_HA.fasta.mu.trim05","InfluenzaBvirus_NA.fasta.mu.trim05","JCpolyomavirus_VP1.fasta_pruned.mu.trim05"
      ,"Measles.fasta_pruned.mu.trim05","Measles_hemagglutinin_OR_haemagglutinin.fasta_pruned.mu.trim05","ParvovirusB19_NS1.fasta_pruned.mu.trim05","ParvovirusB19_VP1.fasta_pruned.mu.trim05","RhinovirusB.fasta_pruned.mu.trim05","RhinovirusB_polyprotein.fasta_pruned.mu.trim05"
      ,"RhinovirusC.fasta_pruned.mu.trim05","RotavirusA_VP6.fasta_pruned.mu.trim05")
source("R_scripts/sort.R")
#source("R_scriptssort_CSV.R")
fasta2<-c("DengueVirus1.fasta_pruned.mu.trim05","DengueVirus2.fasta_pruned.mu.trim05","DengueVirus3.fasta_pruned.mu.trim05","DengueVirus4.fasta_pruned.mu.trim05")
for (i in fasta2){
  if (i == "~$fluenzaAvirus_NA_H1N1.fasta.mu.trim05") {
    next
   }
  # if (i == "DengueVirus1.fasta_pruned.mu.trim05") {
  #   next
  # }
  # if (i == "DengueVirus2.fasta_pruned.mu.trim05") {
  #   next
  # }
  # if (i == "DengueVirus3.fasta_pruned.mu.trim05") {
  #   next
  # }
  # if (i == "DengueVirus4.fasta_pruned.mu.trim05") {
  #   next
  # }
  print(i)
  sort_CSV(i)
  
}
source("R_scripts/graphs/redoplot.R")
truename <- list('DengueVirus1', 'DengueVirus2', 'DengueVirus3', 'DengueVirus4', 'humanparainfluenzavirus1_F', 'humanparainfluenzavirus1_HN', 'humanparainfluenzavirus3_HN', 'InfluenzaAvirus_HA_H1N1','InfluenzaAvirus_HA_H3N2', 'InfluenzaAvirus_NA_H1N1', 'InfluenzaAvirus_NA_H3N2','InfluenzaBvirus_HA', 'InfluenzaBvirus_NA', 'EnterovirusA_VP1', 'EnterovirusA_VP2','EnterovirusB_VP1', 'EnterovirusB_VP2','EnterovirusC_VP1','EnterovirusC_VP2','EnterovirusD_VP1', 'HumanBocavirus1_NS1', 'HumanBocavirus1_VP1', 'HepatitisB', 'Humanherpesvirus2_glycoprotein_G', 'Humanpapillomavirus16', 'Humanpapillomavirus16_L1', 'Humanrespiratorysyncytialvirus', 'Humanrespiratorysyncytialvirus_G', 'JCpolyomavirus_VP1', 'Measles','Measles_hemagglutinin_OR_haemagglutinin','ParvovirusB19_NS1', 'ParvovirusB19_VP1', 'RhinovirusB', 'RhinovirusB_polyprotein', 'RhinovirusC', 'RotavirusA_VP6', 'humanparainfluenzavirus1', 'humanparainfluenzavirus3','BKpolyomavirus_VP1','HepatitisB_core','HepatitisB_s','HepatitisB_precore','HepatitisB_pre_S','HepatitisB_polymerase_truncated_precore','HepatitisB_polymerase')
for (i in truename){
  
  print(i)
  comparing_CpG_Syn_Nonsyn_new(i)
  
}

