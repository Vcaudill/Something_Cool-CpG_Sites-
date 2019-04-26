#library("ape")
#vp1<-read.dna("data/fasta/BKpolyomavirus_VP1.fasta_pruned.mu.trim05", format="fasta", as.character=TRUE)
#vp1.trimmed<-vp1[1:113]
#write.fasta(vp1.trimmed,open="w",names=BKpolyomavirus_VP1,)
#vp1[1]
Virus_info <- read.csv("data/CpG_List.csv")
library(ape)
library(seqinr)
k=5
for(k in 1:nrow(Virus_info)){
  file<-Virus_info[k,1]
  not_same_len=c()
  path_file<-paste("output/FastaSample/",'123',file,".fasta",sep="")
  al2<-read.fasta(path_file, as.string = TRUE)
  number_of_nucs <-nchar(al2[[1]][1])
  for(i in 1:length(al2)){
    print(al2[[2]][1])
    }
    print(nchar(al2[[2]][1]))
    if(nchar(al2[[i]][1])!=number_of_nucs)
      not_same_len<-c(not_same_len,i)}
  
Virus_info<- read.csv("data/CpG_List_SarVersion.csv")
# k 31 Humanrespiratorysyncytialvirus_G.fasta.mu.trim05, 32 Humanrespiratorysyncytialvirus.fasta.mu.trim05
# 41 Measles.fasta_pruned.mu.trim05.fasta no data, 45 RhinovirusB.fasta_pruned.mu.trim05.fasta"
for(k in 51:52){
file<-Virus_info[k,1]
path_file<-paste("output/length_testing/",file,".fasta",sep="")
#path_file<-paste("data/trimmed_no_stop/",file,".fasta",sep="")
DataSet <- read.fasta(path_file)


sample_size= 200
if (Virus_info$SeqNumber[k] > sample_size) {
  DS <- sample(DataSet, size = sample_size)
  size = sample_size
  }else{
    DS <- sample(DataSet, size = length(DataSet))
    size = length(DataSet)
}
Fasta_trimmed_sample<-paste("output/Length_Sample/",Virus_info$name[k],size,".fasta",sep="")
write.fasta(DS,names(DS), file = Fasta_trimmed_sample)
}
#call file in
#find number of sequences 

