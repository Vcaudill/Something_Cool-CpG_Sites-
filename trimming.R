
Virus_info <- read.csv("data/CpG_List_SarVersion.csv")
library(ape)
library(seqinr)


for(k in 1:nrow(Virus_info)){
  file<-Virus_info[k,1]
  print(file)
  path_file<-paste("data/fasta/",file, sep="")
  al2<-read.dna(path_file,format="fasta")
  
  number_of_nucs <-length(as.list(al2)[[1]])
  end<-as.numeric(Virus_info$start[k])+as.numeric(Virus_info$length[k])
  al.trimmed<-al2[,c(as.numeric(Virus_info$start[k]):(end))]
  
  Fasta_FF<-paste0("data/trimmed_no_stop/",Virus_info$name[k],".fasta")
  write.FASTA(al.trimmed,file=Fasta_FF )
  ####

  al<-read.fasta(Fasta_FF,as.string=TRUE)
  
  Fasta_final<-paste("data/no_mid_stop/",Virus_info$name[k],".fasta",sep="")
  
  pos<-c()
  stop.seq.no<-c()
  for (i in 1: length(al)){
    seq<-paste(al[i])
    #split into 3 bases from the start (if frame1)
    seq<-substring(seq,seq(1,nchar(seq),3), seq(3, nchar(seq),3) )
    p<-which(seq =="taa"|seq=="tga"|seq=="tag")
    if (length(p)!=0){
      for (j in 1:length(p)){
        if (p[j]==length(seq)) pos<-c(pos,p[j])
        else  stop.seq.no<-c(stop.seq.no,i)
      } 
    }
    
  }
  poss<-unique(pos)
  
  # if stop.seq.no is empty, al[-stop.seq.no] will erase the whole file. So let it go throug the if loops. At the end, 
  # you will have "xxxxx.trimmed.fasta" regardless whether you trimmed/removed seqeunce(s)
  length(unique(stop.seq.no))
  if (length(stop.seq.no)!=0 ){
    al<-al[-stop.seq.no]
    write.fasta(al,names(al),file=Fasta_final)
    
    if (length(poss)!=0){
      al.filtered<-read.dna(Fasta_FF, format = "fasta")
      al.trimmed<-al.filtered[,-c(poss*3-2, poss*3-1,poss*3)]
      write.FASTA(al.trimmed, file=Fasta_final)} 
  }
  
  if (length(stop.seq.no)==0){
    if (length(poss)!=0){
      al.filtered<-read.dna(Fasta_FF, format = "fasta")
      al.trimmed<-al.filtered[,-c(poss*3-2, poss*3-1,poss*3)]
      write.FASTA(al.trimmed, file= Fasta_final )}
    else file.copy(Fasta_FF, Fasta_final )
  }
  
}
