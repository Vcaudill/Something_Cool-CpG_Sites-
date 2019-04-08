
Virus_info <- read.csv("data/CpG_List.csv")
library(ape)
library(seqinr)

k=10

for(k in 1:nrow(Virus_info)){
  file<-Virus_info[k,1]
  path_file<-paste("data/fasta/",file, sep="")
  al2<-read.dna(path_file,format="fasta")
  
  number_of_nucs <-length(as.list(al2)[[1]])
  end<-as.numeric(number_of_nucs)-as.numeric(Virus_info$stop[k])
  al.trimmed<-al2[,c(as.numeric(Virus_info$start[k]):(end-3))]
  
  Fasta_FF<-paste0("data/trimmed_no_stop/",Virus_info$name[k],".fasta")
  write.FASTA(al.trimmed,file=Fasta_FF )
  ####

  al<-read.fasta(Fasta_FF,as.string=TRUE)
  
  pos<-c() 
  stop.seq.no<-c()
  for (i in 1: length(al)){
    seq<-paste(al[i])
    #split into 3 bases from the start (if frame1)
    seq<-substring(seq,seq(1,nchar(seq),3), seq(3, nchar(seq),3) )
    p<-which(seq =="taa"|seq=="tga"|seq=="tag")
    if (length(p)!= 0){
      for (j in 1:length(p)){
        if (p[j]==length(seq)) pos<-c(pos,p[j])
        else  stop.seq.no<-c(stop.seq.no,i)
          
      }
      }
  }
  if (length(stop.seq.no)!=0) al.removed<-al[-unique(stop.seq.no)]
  Fasta_final<-paste("data/no_mid_stop/",Virus_info$name[k],".fasta",sep="")
  write.fasta(al.removed,names(al.removed),file=Fasta_final )
}
