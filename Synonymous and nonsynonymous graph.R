#Read in CSV file
Dengue<-read.csv("DengueVirus1.fasta_pruned.mu.trim05_DF.csv")
SynNonsynAT<-function(Virus){
#select for only synonymous sites in Virus file
synonymous<-which(Virus$TypeOfSite=='syn')
#subset data into synonymous sites and designated wildtype nucleotide
synonymousA<-subset(Virus, wtnt=='a' & TypeOfSite=='syn')
synonymousT<-subset(Virus, wtnt=='t' & TypeOfSite=='syn') 
#binds synonymous subset data into a new dataset
synonymous<-rbind(synonymousA, synonymousT)
#select for only nonsynonymous sites in Virus file
nonsynonymous<-which(Virus$TypeOfSite=='nonsyn')
#subset data into nonsynonymous sites and designated wildtypes
nonsynonymousA<-subset(Virus, wtnt=='a' & TypeOfSite=='nonsyn')
nonsynonymousT<-subset(Virus, wtnt=='t' & TypeOfSite=='nonsyn')

#A and T synonymous and nonsynonymous graphs

#splits the graph into 4 sections, 2 rows 2 columns 
layout(matrix(c(1,2,3,4), nrow=2, byrow = TRUE))
#log meanfreq to get rid of 0's and make graph easier to read
plot(synonymousA$num,synonymousA$MeanFreq+.000001,log='y', col = "red", 
     #pch creates the shape of the plotted points
     ,xlab="Position Number", ylab="MeanFreq",main="Synonymous A", pch=2)
#abline creates a line where the mean of all the points are located
abline(h=mean(synonymousA$MeanFreq+.000001))

#plots continue 
plot(nonsynonymousA$num,nonsynonymousA$MeanFreq+.000001,log='y', col = "blue",xlab="Position Number", ylab="MeanFreq",main="NonSynonymous A", pch = 5)
abline(h=mean(nonsynonymousA$MeanFreq+.000001))

plot(synonymousT$num,synonymousT$MeanFreq+.000001,log='y', col = "green",xlab="Position Number", ylab="MeanFreq",main="Synonymous T")
abline(h=mean(synonymousT$MeanFreq+.000001))

plot(nonsynonymousT$num,nonsynonymousT$MeanFreq+.000001,log='y', col = "purple",xlab="Position Number", ylab="MeanFreq",main="NonSynonymous T", pch=0)
abline(h=mean(nonsynonymousT$MeanFreq+.000001))

}
SynNonsynCG<-function(Virus2){
  
#subset data into synonymous sites and designated wildtypes
synonymousC<-subset(Virus2, wtnt=='c' & TypeOfSite=='syn')
synonymousG<-subset(Virus2, wtnt=='g' & TypeOfSite=='syn')
#subset data into nonsynonymous sites and designated wildtypes
nonsynonymousC<-subset(Virus2, wtnt=='c' & TypeOfSite=='nonsyn')
nonsynonymousG<-subset(Virus2, wtnt=='g' & TypeOfSite=='nonsyn')
  
#Plot C and G synonymous and nonsynonymous graphs
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
plot(synonymousC$num,synonymousC$MeanFreq+.000001,log='y', col = "deeppink4",xlab="Position Number", ylab="MeanFreq",main="Synonymous C", pch=7)
abline(h=mean(synonymousC$MeanFreq+.000001))

plot(nonsynonymousC$num,nonsynonymousC$MeanFreq+.000001,log='y', col = "yellow4",xlab="Position Number", ylab="MeanFreq",main="NonSynonymous C", pch = 9)
abline(h=mean(nonsynonymousC$MeanFreq+.000001))

plot(synonymousG$num,synonymousG$MeanFreq+.000001,log='y', col = "palegreen4",xlab="Position Number", ylab="MeanFreq",main="Synonymous G",pch=13)
abline(h=mean(synonymousG$MeanFreq+.000001))

plot(nonsynonymousG$num,nonsynonymousG$MeanFreq+.000001,log='y', col = "dodgerblue4" ,xlab="Position Number", ylab="MeanFreq",main="NonSynonymous G", pch=14)
abline(h=mean(nonsynonymousG$MeanFreq+.000001))

}  
SynNonsynAT(Dengue)
SynNonsynCG(Dengue)

plot(1, type="n", axes=FALSE, xlab="", ylab="") 

legend("center", legend=c
       ("synonymous A","nonsynonymous A","synonymous T","nonsynonymous T","synonymous C","nonsynonymous C","synonymous G","nonsynonymous G"),
       col= c("red","blue","green","purple","deeppink4","yellow4","palegreen4","dodgerblue4"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1,xpd=TRUE,bty='n', pch = c(2,5,1,0,7,9,13,14))
