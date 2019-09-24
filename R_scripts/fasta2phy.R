library(ape)
library(seqinr) 
#install.packages("phylotools")
library(phylotools) 
library(dplyr)
Virus_info<- read.csv("data/list/Final_CpG_list.csv")

for(i in 1:nrow(Virus_info)){
  viruplace = paste(Virus_info$Fasta_File_Path[i],"/fasta/", Virus_info$File[i], sep="")
  name = as.character(Virus_info$File[i])
  splitname<-unlist(strsplit(as.character(Virus_info$File[i]),".fasta"))
  truename<-splitname[1]
  
  sample_size= 200 #sequences
  DataSet <- read.fasta(viruplace)

  if (Virus_info$Number_of_Sequences[i] > sample_size) {
    DS <- dplyr::sample_n(DataSet, size = sample_size, replace=F)
    size = sample_size
  }else{
    DS <- dplyr::sample_n(DataSet, size = as.numeric(nrow(DataSet)))
  }
  
  
  Fasta_trimmed_sample<-paste("data/tree_data/samlpe/",truename,"_sample.fasta",sep="")
  dat2fasta(DS, outfile = Fasta_trimmed_sample)
  
  output <- paste("data/tree_data/phy/", truename, ".phy",  sep="")
  sampleviruplace = paste("data/tree_data/samlpe/",truename,"_sample.fasta",sep="")
  dat2phylip(read.fasta(sampleviruplace), outfile= output)
}

