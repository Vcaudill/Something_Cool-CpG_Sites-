library(seqinr)
fasta<-"EnterovirusC_VP1.fasta_pruned.mu.trim05"
fasta_file <-paste("data/fasta/",fasta,sep = "")
virus_basic <- read.fasta(fasta_file, as.string = TRUE)
IDs <-names(virus_basic)
for (i in 1:length(virus_basic)){  
  enteroC<-grepl("poliovirus_3", format(IDs[i]))
  if(enteroC==FALSE){
    next
  }
  write.fasta(virus_basic[i], open = "a", names = names(virus_basic[i]), file = paste("output/","poliovirus_3","_only.fasta",sep=""))
  print(enteroC)
}
#||Enterovirus_C||Human_enterovirus
#Human_poliovirus_1||Human_poliovirus_2||Human_poliovirus_3
#Human_coxsackievirus