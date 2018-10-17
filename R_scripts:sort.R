library(seqinr) 
fasta_file <-"data/fasta/BKpolyomavirus_VP1.fasta.mu.trim05"
  virus_basic <- read.fasta(fasta_file, as.string = TRUE)
  number_of_seqs <- length(virus_basic)
  
  virus_align <- read.alignment(fasta_file, format = "fasta", forceToLower = T)
  virus_consensus <- seqinr :: consensus(virus_align, method = "majority")
  virus_consensus_matrix <- seqinr :: consensus(virus_align, method = "profile")
  consensus_length <- length(virus_consensus)
  number_column <- seq(1, consensus_length)
  consensus<-paste(virus_consensus,collapse = '')
  virus_DF <- data.frame("num" = 1:number_of_seqs, "Differences" = 0, sequencename= row.names(summary(virus_basic)))
  
  virus_basic1 <- read.fasta(fasta_file, as.string = TRUE)
for(i in 1:number_of_seqs){
  virus_DF$Differences[i]=adist(consensus,virus_basic1[[i]][1])
}

  
