library(seqinr)
fasta<-"InfluenzaAvirus_HA_H1N1.fasta.mu.trim05"
Short_fasta<-function(fasta){
  fasta_file <-paste("data/fasta/",fasta,sep = "")
  virus_basic <- read.fasta(fasta_file, as.string = TRUE)
  number_of_seqs <- length(virus_basic)
  truename = "InfluenzaAvirus_HA_H1N1"
  virus_align <- read.alignment(fasta_file, format = "fasta", forceToLower = T)
  virus_consensus <- seqinr :: consensus(virus_align, method = "majority")
  virus_consensus_matrix <- seqinr :: consensus(virus_align, method = "profile")
  consensus_length <- length(virus_consensus)
  number_column <- seq(1, consensus_length)
  consensus<-paste(virus_consensus,collapse = '')
  virus_DF <- data.frame("num" = 1:number_of_seqs, "Differences" = 0, sequencename= row.names(summary(virus_basic)))
  
  virus_basic1 <- read.fasta(fasta_file, as.string = TRUE)
  write.fasta("", open = "w", names = "test", file = paste("output/ShortFasta/",truename,"_short.fasta",sep=""))
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
Short_fasta(fasta)
