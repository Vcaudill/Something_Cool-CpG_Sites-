ListofFastaFiles <-list.files("data/fasta/", pattern = "trim05", all.files = TRUE)
library(seqinr) #package phylo for making trees
library(ggtree)
library(ape)
#print(ListofFastaFiles[26])
#filename = ListofFastaFiles[26]
samplesize = 100
for (filename in ListofFastaFiles) {
  print(filename)
  if (filename == "~$fluenzaAvirus_NA_H1N1.fasta.mu.trim05") {
    next
  }
  DataSet <- read.fasta(paste("data/fasta/", filename, sep = ""))
  if (length(DataSet) > samplesize) {
    DS <- sample(DataSet, size = samplesize)
    write.fasta(DS, names(DS), file = paste("data/sample/",samplesize,filename,sep = "_"))
    file<-paste("data/sample/",samplesize,filename,sep = "_")
    NAseqs<-ape::read.dna(file, format = "fasta")
    D<-dist.dna(NAseqs) #create a distance matrix
    # 
    NJ<-njs(D) #create a neighbor joining tree
    DSpng = paste("data/sample/",filename,".png",sep="")
    png(DSpng, width = 6.75, height = 6.75, units = "in", res= 300)
    ggtree(NJ)
    dev.off()
    # 
  }else{
    DS <- sample(DataSet, size = length(DataSet))
    write.fasta(DS, names(DS), file = paste("data/sample/",samplesize,filename,sep = "_"))
    file<-paste("data/sample/",samplesize,filename,sep = "_")
    NAseqs<-ape::read.dna(file, format = "fasta")
    D<-dist.dna(NAseqs) #create a distance matrix
    # 
    NJ<-njs(D) #create a neighbor joining tree
    DSpng = paste("data/sample/",filename,".png",sep="")
    png(DSpng, width = 6.75, height = 6.75, units = "in", res= 300)
    ggtree(NJ)
    dev.off()

  }
  
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
