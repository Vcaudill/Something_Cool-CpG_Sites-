library(seqinr)
#fasta<-"Measles.fasta_pruned.mu.trim05"
sort_CSV<-function(fasta){
fasta_file <-paste("data/fasta/",fasta,sep = "")
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
  splitname<-unlist(strsplit(fasta,".fasta"))
  truename<-splitname[1]
  write.csv(virus_DF, file = paste("data/sort_CSV/",truename,".csv",sep=""))
}
  
