#library("ape")
#vp1<-read.dna("data/fasta/BKpolyomavirus_VP1.fasta_pruned.mu.trim05", format="fasta", as.character=TRUE)
#vp1.trimmed<-vp1[1:113]
#write.fasta(vp1.trimmed,open="w",names=BKpolyomavirus_VP1,)
#vp1[1]
Virus_info <- read.csv("data/CpG_List.csv")
library(ape)
library(seqinr)

for(j in 1:nrow(Virus_info)){
file<-Virus_info[j,1]
path_file<-paste("data/fasta/",file, sep="")
al<-read.fasta(path_file,as.string=TRUE)
al2<-read.fasta(path_file,as.string=FALSE)
number_of_nucs <-length(al2[[1]])
end<-as.numeric(number_of_nucs)-as.numeric(Virus_info$stop[j])
al.trimmed<-al[as.numeric(Virus_info$start[j]):end]


pos<-c()
for (i in 1: length(al)){
  seq<-paste(al[i])
  #split into 3 bases from the start (if frame1)
  seq<-substring(seq,seq(1,nchar(seq),3), seq(3, nchar(seq),3) )
  p<-which(seq =="taa"|seq=="tga"|seq=="tag")
  pos<-c(pos,p)
  
}

poss<-unique(pos)
pos.remove<-c(poss*3-2, poss*3-1,poss*3) 

algn1<-read.dna(path_file, format = "fasta")
algn.trimmed<-algn1[,-pos.remove]
output_file<-paste("output/NoStopFasta/",file,"_No_Stop.fasta",sep="")
write.FASTA(algn.trimmed, file=output_file)
}

k=1
for(k in 1:nrow(Virus_info)){
file<-Virus_info[j,1]
path_file<-paste("output/NoStopFasta/",file,"_No_Stop.fasta",sep="")
DataSet <- read.fasta(path_file)
sample_size= 499
if (Virus_info$SeqNumber[k] > sample_size) {
  DS <- sample(DataSet, size = sample_size)
  size = sample_size
}else{
  DS <- sample(DataSet, size = Virus_info$SeqNumber[k])
  size = Virus_info$SeqNumber[k]
}
Fasta_trimmed_sample<-paste("output/FastaSample/",size,Virus_info$name[k],".fasta",sep="")
write.fasta(DS,names(DS), file = Fasta_trimmed_sample)
}
#call file in
#find number of sequences 

