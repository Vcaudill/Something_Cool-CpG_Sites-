#library("ape")
#vp1<-read.dna("data/fasta/BKpolyomavirus_VP1.fasta_pruned.mu.trim05", format="fasta", as.character=TRUE)
#vp1.trimmed<-vp1[1:113]
#write.fasta(vp1.trimmed,open="w",names=BKpolyomavirus_VP1,)
#vp1[1]
Virus_info <- read.csv("data/CpG_List.csv")
library(ape)
library(seqinr)

for(i in 1:length(Virus_info)){
file<-Virus_info[i,1]
path_file<-paste("data/fasta/",file, sep="")
al<-read.fasta(path_file,as.string=TRUE)
end<-as.numeric(Virus_info$NucleotideNumber[i])-as.numeric(Virus_info$stop[i])
al.trimmed<-al[as.numeric(Virus_info$start[i]):end]

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
