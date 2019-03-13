library(seqinr)
fasta<-"BKpolyomavirus_VP1.fasta.mu.trim05"
Short_fasta<-function(fasta){
  fasta_file <-paste("data/fasta/",fasta,sep = "")
  virus_basic <- read.fasta(fasta_file, as.string = TRUE)
  number_of_seqs <- length(virus_basic)
  virus_align <- read.alignment(fasta_file, format = "fasta", forceToLower = T)
  DF<-Freq(fasta)
  splitname<-unlist(strsplit(fasta,".fasta"))
  truename<-splitname[1]
  DF$wtnt<-as.character(DF$wtnt)
  DF$Virus<-(truename)
  
  start = 1
  end = 0
  DF<-DF[c(start:(nrow(DF)-end)),]
  DF$num<-(1:nrow(DF))
  source("R_scripts/BioInfo/WTAA_consensus.R")
  DF<-getWTAA(DF)
  
  source("R_scripts/BioInfo/MUTAA.R")
  DF<-getMUTAA(DF)
  write.fasta(virus_basic, open = "a", names = truename, file = paste("output/ShortFasta/",truename,"_short.fasta",sep=""))
    
  }
}
Short_fasta(fasta)
