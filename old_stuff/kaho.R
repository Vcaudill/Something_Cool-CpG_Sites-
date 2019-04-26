setwd("~/Desktop/genecut")
library(ape)
library(seqinr)

al<-read.fasta("ParvovirusB19_VP1_genecut.fasta",as.string=TRUE)
al<-read.fasta("499HepatitisB_core.fasta_pruned.mu.trim05.fasta",as.string=TRUE)
al<-read.fasta(path_file,as.string=TRUE)
al<-read.fasta("data/fasta/HIV1_FLT_2017_pol_DNA.fasta",as.string=TRUE)
al<-read.fasta("data/fasta/EnterovirusC_VP1.fasta_pruned.mu.trim05",as.string=TRUE)
al<-read.fasta("output/length_testing/EnterovirusC_VP1.fasta_pruned.mu_more_aligned.fasta",as.string=TRUE)
i=29
pos<-c()
whichseq<-c()
for (i in 1: length(al)){
  seq<-paste(al[i])
  #split into 3 bases from the start (if frame1)
  seq<-substring(seq,seq(1,nchar(seq),3), seq(3, nchar(seq),3) )
  p<-which(seq =="taa"|seq=="tga"|seq=="tag")
  print(c("sequence",i,"stop points",p))
  if (length(p)!=0){
    whichseq<-c(whichseq,c("sequence",i,"stop points",p))}
  pos<-c(pos,p)
  
}

poss<-unique(pos)
whichseq
pos.remove<-c(poss*3-2, poss*3-1,poss*3) 
+
algn1<-read.dna("BKpolyomavirus_VP1.fasta", format = "fasta")
algn.trimmed<-algn1[,-pos.remove]
write.FASTA(algn.trimmed, file="Trimmed.fasta" )

strsplit(al.removed[[1]][1])

######################################
#### this removes the full sequence####
pos<-c() 
stop.seq.no<-c()
for (i in 1: length(al)){
  seq<-paste(al[i])
  #split into 3 bases from the start (if frame1)
  seq<-substring(seq,seq(1,nchar(seq),3), seq(3, nchar(seq),3) )
  p<-which(seq =="taa"|seq=="tga"|seq=="tag")
  
  for (j in 1:length(p)){
    if (p[j]==length(seq)) pos<-c(pos,p[j])
    else  stop.seq.no<-c(stop.seq.no,i)
  }
}

al.removed<-al[-stop.seq.no]
write.fasta(al.removed,names(al.removed),"BKpolyomavirus.trimmed.fasta" )


