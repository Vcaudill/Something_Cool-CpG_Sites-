#library("ape")
#vp1<-read.dna("data/fasta/BKpolyomavirus_VP1.fasta_pruned.mu.trim05", format="fasta", as.character=TRUE)
#vp1.trimmed<-vp1[1:113]
#write.fasta(vp1.trimmed,open="w",names=BKpolyomavirus_VP1,)
#vp1[1]

library(ape)
library(seqinr)

al<-read.fasta("data/fasta/DengueVirus1.fasta_pruned.mu.trim05",as.string=TRUE)

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

algn1<-read.dna("data/fasta/DengueVirus1.fasta_pruned.mu.trim05", format = "fasta")
algn.trimmed<-algn1[,-pos.remove]
write.FASTA(algn.trimmed, file="output/NoStopFasta/DengueVirus1.fasta")
