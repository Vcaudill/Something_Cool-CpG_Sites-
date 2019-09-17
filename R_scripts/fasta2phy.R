library(ape)
library(seqinr) 
#install.packages("phylotools")
library(phylotools) 

Virus_info<- read.csv("data/list/Final_CpG_list.csv")

for(i in 1:nrow(Virus_info)){
 
  viruplace = paste(Virus_info$Fasta_File_Path[i],"/fasta/", Virus_info$File[i], sep="")
  name = as.character(Virus_info$File[i])
  splitname<-unlist(strsplit(as.character(Virus_info$File[i]),".fasta"))
  truename<-splitname[1]
  output <- paste("data/tree_data/phy/", truename, ".phy",  sep="")

  dat2phylip(read.fasta(viruplace), outfile= output)
}

