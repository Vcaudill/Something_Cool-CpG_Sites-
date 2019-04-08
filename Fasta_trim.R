#library("ape")
#vp1<-read.dna("data/fasta/BKpolyomavirus_VP1.fasta_pruned.mu.trim05", format="fasta", as.character=TRUE)
#vp1.trimmed<-vp1[1:113]
#write.fasta(vp1.trimmed,open="w",names=BKpolyomavirus_VP1,)
#vp1[1]
Virus_info <- read.csv("data/CpG_List.csv")
library(ape)
library(seqinr)

for(k in 1:nrow(Virus_info)){
file<-Virus_info[k,1]
path_file<-paste("data/no_mid_stop/",file,".fasta",sep="")
DataSet <- read.fasta(path_file)
sample_size= 499
if (Virus_info$SeqNumber[k] > sample_size) {
  DS <- sample(DataSet, size = sample_size)
  size = sample_size
  }else{
    DS <- sample(DataSet, size = length(DataSet))
    size = length(DataSet)
}
Fasta_trimmed_sample<-paste("output/FastaSample/",size,Virus_info$name[k],".fasta",sep="")
write.fasta(DS,names(DS), file = Fasta_trimmed_sample)
}
#call file in
#find number of sequences 

