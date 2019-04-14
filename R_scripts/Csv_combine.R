library(seqinr)

hyphy_virus<-read.csv("data/MEME_CpG_list1.csv")

for (i in 1: nrow(hyphy_virus)){
#these viruses dont have datamonkey files yet, except humanpap not using
if(hyphy_virus$name[i]== "DengueVirus1.fasta_pruned.mu.trim05"){
  next
}
if(hyphy_virus$name[i]== "DengueVirus2.fasta_pruned.mu.trim05"){
  next
}
if(hyphy_virus$name[i]== "DengueVirus3.fasta_pruned.mu.trim05"){
  next
}
if(hyphy_virus$name[i]== "InfluenzaAvirus_HA_H3N2.fasta.mu.trim05"){
    next
}
if(hyphy_virus$name[i]== "Humanpapillomavirus16.fasta_pruned.mu.trim05"){
    next
}
if(hyphy_virus$name[i]== "InfluenzaBvirus_HA.fasta.mu.trim05"){
    next
  }
  
  
  

name = as.character(hyphy_virus$name[i])
splitname<-unlist(strsplit(as.character(hyphy_virus$name[i]),".fasta"))
truename<-splitname[1]  
path_file<-paste("Hyphy/Consensus_Hyphy/",truename,".csv", sep="")
al1<-read.csv(path_file)
folder<-paste(hyphy_virus$Hyphy_folders[i],sep='')
path<-paste("Hyphy/",folder,"/datamonkey-table.csv",sep='')
datamonkey<-read.csv(path)
file<-al1$wtnt_consensus
nuc<-splitseq(file, frame = 0, word = 3)
makesCpG<-al1$makesCpG
CpG<-splitseq(makesCpG, frame = 0, word = 3)

datamonkey$consensus_codon<-nuc
datamonkey$consensus_AA<-translate(as.character(file))
datamonkey$potential_CpG<-"no"
datamonkey$makesCpG<-0
for(j in 1:nrow(datamonkey)){
  letter<-datamonkey$consensus_AA[j]
  if(letter=="H"||letter=="E"||letter=="Y"||letter=="C"||letter=="F"||letter=="N"||letter=="K"||letter=="Q"||letter=="D"){
    datamonkey$potential_CpG[j]="yes"
  }
  if(sum(as.numeric(CpG[j])) > 0){
    datamonkey$makesCpG[j] = 1
  }
}

write.csv(datamonkey,path)

}

  
