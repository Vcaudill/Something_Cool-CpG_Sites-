ListofFastaFiles<-list.files(pattern = "trim05")  #do not take if csv
library(seqinr) #package phylo for making trees
samplesize=50
for(filename in ListofFastaFiles){
  print(filename)
  DataSet <-read.fasta(filename)
  DS<-sample(DataSet,size=samplesize)
  write.fasta(DS, names(DS), file = paste("DS",samplesize,filename,sep = "_"))
}

#DataSet <-read.fasta("InfluenzaAvirus_HA_H1N1.fasta.mu.trim05")
DataSet <-read.fasta("DengueVirus4.fasta_pruned.mu.trim05")
#DataSet <-read.fasta("InfluenzaBvirus_HA.fasta.mu.trim05")
#View(DataSet)
DS<-sample(DataSet,size=50)
#View(DS)
write.fasta(DS, names(DS), file = "Den4")
ListofFastaFiles<-list.files(pattern = "trim05")  
