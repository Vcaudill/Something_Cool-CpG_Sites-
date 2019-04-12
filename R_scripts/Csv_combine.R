library(seqinr)
Virus_info <- read.csv("data/MEME_CpG_list.csv")
k=2
i=1

for(k in 1:nrow(Virus_info)){
name = as.character(Virus_info$name[k])
splitname<-unlist(strsplit(as.character(Virus_info$name[k]),".fasta"))
truename<-splitname[1]  
path_file<-paste("Hyphy/Consensus_Hyphy/",truename,".csv", sep="")
al<-read.csv(path_file)
#for(j in 1:length(Virus_info)){
  
  #list.files(path = ".", pattern = NULL, all.files = FALSE,
  #           full.names = FALSE, recursive = FALSE,
  #           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  #dir(path = ".", pattern = NULL, all.files = FALSE,
  #    full.names = FALSE, recursive = FALSE,
  #    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  #list.dirs(path = ".", full.names = TRUE, recursive = TRUE
  foldername<-list.files("Hyphy/", pattern = "stuff", all.files = TRUE)
  Virus_info$Hyphy_folders<-c(foldername)
  write.csv(Virus_info, file = "data/MEME_CpG_list1.csv")
}  #folder HepB_s, InfluenzaBvirus_Ha
i=1
hyphy_virus<-read.csv("data/MEME_CpG_list1.csv")
for (i in 1: length(hyphy_virus)){
if(hyphy_virus$name[i]== "DengueVirus1.fasta_pruned.mu.trim05"){
  next
}
if(hyphy_virus$name[i]== "DengueVirus2.fasta_pruned.mu.trim05"){
  next
}
if(hyphy_virus$name[i]== "DengueVirus3.fasta_pruned.mu.trim05"){
  next
}
if(hyphy_virus$name[i]== ".fasta_pruned.mu.trim05"){
    next
}
name = as.character(hyphy_virus$name[i])
splitname<-unlist(strsplit(as.character(Virus_info$name[i]),".fasta"))
truename<-splitname[1]  
path_file<-paste("Hyphy/Consensus_Hyphy/",truename,".csv", sep="")
al1<-read.csv(path_file)
folder<-paste(hyphy_virus$Hyphy_folders[i],sep='')
path<-paste("Hyphy/",folder,"/datamonkey-table.csv",sep='')
datamonkey<-read.csv(path)
file<-al1$wtnt_consensus
nuc<-splitseq(file, frame = 0, word = 3)

datamonkey$wtnt_codon<-nuc

write.csv(datamonkey,path)

}

  
