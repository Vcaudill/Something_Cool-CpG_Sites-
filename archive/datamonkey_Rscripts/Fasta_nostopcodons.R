library(seqinr)

Virus_info<- read.csv("data/CpG_List.csv")

fasta<-"BKpolyomavirus_VP1.fasta.mu.trim05"
i = 1
No_STOP_fasta<-function(fasta, i){
  fasta_file <-paste("data/fasta/",fasta,sep = "")
  virus_basic <- read.fasta(fasta_file)
  number_of_seqs <- length(virus_basic)
  #cat("number_of_seqs", number_of_seqs)
  virus_align <- read.alignment(fasta_file, format = "fasta", forceToLower = T)
  virus_consensus <- seqinr :: consensus(virus_align, method = "majority")
  virus_consensus_matrix <- seqinr :: consensus(virus_align, method = "profile")
  consensus_length <- length(virus_consensus)
  number_column <- seq(1, consensus_length)
  viruplace = paste('data/fasta/', Virus_info$name[i], sep="")
  
  DF<-Freq(viruplace)
  name = as.character(Virus_info$name[i])
  splitname<-unlist(strsplit(as.character(Virus_info$name[i]),".fasta"))
  truename<-splitname[1]
  DF$wtnt_consensus<-as.character(DF$wtnt_consensus)
  DF$Virus<-(truename)
  #setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")
  
  #  
  # make sure that M is the first amino acid check all virus
  start = Virus_info$start[i]
  end = Virus_info$stop[i] # add three or sub
  DF<-DF[c(start:(nrow(DF)-end)),]
  DF$num<-(1:nrow(DF))
  
  virus_basic1 <- read.fasta(fasta_file, as.string = TRUE)
  write.fasta("", open = "w", names = "test", file = paste("output/NoStopFasta/",truename,"_short.fasta",sep=""))
  splitname<-unlist(strsplit(fasta,".fasta"))
  truename<-splitname[1]
  for(i in 1:number_of_seqs){
    print(i/number_of_seqs)
    virus_DF$Differences[i]=adist(consensus,virus_basic1[[i]][1])
    if(virus_DF$Differences[i] > 99){
      next
    }
    write.fasta(virus_basic[i], open = "a", names = names(virus_basic[i]), file = paste("output/ShortFasta/",truename,"_short.fasta",sep=""))
    
  }
}


#####

# read in the file individulally a a g t a     not a sting aagta
# cut off the begining if need and the stop 
# translate * make sure there are no stop codons 
# wrint in new file 
