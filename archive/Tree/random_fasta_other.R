ListofFastaFiles <-list.files("data/fasta/", pattern = "_", all.files = TRUE)
library(seqinr) #package phylo for making trees
library(ggtree)
library(ape)
#print(ListofFastaFiles[26])
#filename = ListofFastaFiles[26]

Treegraph<- function(DSpng, NJ){
  png(DSpng, width = 6.75, height = 6.75, units = "in", res= 300)
  plot(ggtree(NJ) +geom_tiplab(size=1)+ labs(title=filename))
  dev.off()
  
  DS <- sample(DataSet, size = 100)
  write.fasta(DS, names(DS), file = paste("data/sample/",samplesize,filename,sep = "_"))
  file<-paste("data/sample/",samplesize,filename,sep = "_")
  NAseqs<-ape::read.dna(file, format = "fasta")
  D<-dist.dna(NAseqs) #create a distance matrix
  # 
  NJ<-njs(D)
  plot(ggtree(NJ) +geom_tiplab(color='purple', size=3) + labs(title=filename, caption="powered by Victoria"))
  
  plot(ggtree(NJ, layout = "circular") +geom_tiplab(aes(angle=angle), color='purple', size=2))
  print(DSpng)
  dev.off()
}
samplesize = 100
for (filename in ListofFastaFiles) {
 
  print(filename)
  if (filename == ".DS_Store") {
    next
  }
  if (filename == "IfluenzaAvirus_NA_H1N1.fasta.mu.trim05") {
    next
  }
  if (filename == "EnterovirusC_VP2.fasta_pruned.mu.trim05") {
    next
  }
  if (filename == "Humanrespiratorysyncytialvirus_G.fasta_pruned.mu.trim05") {
    next
  }
  
  if (filename == "Humanrespiratorysyncytialvirus_G.fasta.mu.trim05") {
    next
  }
  if (filename == "EnterovirusB_VP2.fasta.mu.trim08") {
    next
  }
  DataSet <- read.fasta(paste("data/fasta/", filename, sep = ""))
  DSpng = paste("output/tree/",filename,".png",sep="")
  png(DSpng, width = 6.75, height = 6.75, units = "in", res= 300)
  if (length(DataSet) > samplesize) {
    DS <- sample(DataSet, size = samplesize)
    size = samplesize
  }else{
    DS <- sample(DataSet, size = length(DataSet))
    size = length(DataSet)
  }
  write.fasta(DS, names(DS), file = paste("data/sample/",samplesize,filename,sep = "_"))
  file<-paste("data/sample/",samplesize,filename,sep = "_")
  NAseqs<-ape::read.dna(file, format = "fasta")
  D<-dist.dna(NAseqs) #create a distance matrix
  # 
  NJ<-njs(D) #create a neighbor joining tree
  numseq = paste("number of sequences",size)
  plot(ggtree(NJ) +geom_tiplab(size=1)+ labs(title=filename, caption= numseq))
  dev.off()
  # 
  # DSpng2 = paste("output/tree/small_",filename,".png",sep="")
  # png(DSpng2, width = 6.75, height = 6.75, units = "in", res= 300)
  # DS <- sample(DataSet, size = 15)
  # write.fasta(DS, names(DS), file = paste("data/sample/small_",samplesize,filename,sep = "_"))
  # file<-paste("data/sample/small_",samplesize,filename,sep = "_")
  # NAseqs<-ape::read.dna(file, format = "fasta")
  # D<-dist.dna(NAseqs) #create a distance matrix
  # # 
  # NJ<-njs(D)
  # plot(ggtree(NJ) +geom_tiplab(color='purple', size=3) + labs(title=filename, caption="powered by Victoria"))
  # dev.off()
  # DSpng3 = paste("output/tree/circle_",filename,".png",sep="")
  # png(DSpng3, width = 6.75, height = 6.75, units = "in", res= 300)
  # plot(ggtree(NJ, layout = "circular") +geom_tiplab(aes(angle=angle), color='purple', size=2)+ labs(title=filename, caption="powered by Victoria & Sarina"))
  # dev.off()
  
}
DSpng = paste("output/sample/",filename,".png",sep="")
png(DSpng, width = 6.75, height = 6.75, units = "in", res= 300)

#DataSet <-read.fasta("InfluenzaAvirus_HA_H1N1.fasta.mu.trim05")
DataSet <- read.fasta("DengueVirus4.fasta_pruned.mu.trim05")
#DataSet <-read.fasta("InfluenzaBvirus_HA.fasta.mu.trim05")
#View(DataSet)
DS <- sample(DataSet, size = 50)
#View(DS)
write.fasta(DS, names(DS), file = "Den4")
ListofFastaFiles <- list.files(pattern = "trim05")

file <- paste("data/fasta", filename, sep = "")
NAseqs <- ape::read.dna(file, format = "fasta")
D <- dist.dna(NAseqs) #create a distance matrix

NJ <- nj(D) #create a neighbor joining tree
ggtree(NJ)
