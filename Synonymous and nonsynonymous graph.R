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
#binds nonsynonymous subset data into a new dataset
nonsynonymous<-rbind(nonsynonymousA, nonsynonymousT)
#subset data into synonymous sites and designated wildtypes
synonymousC<-subset(Virus, wtnt=='c' & TypeOfSite=='syn')
synonymousG<-subset(Virus, wtnt=='g' & TypeOfSite=='syn')
#subset data into nonsynonymous sites and designated wildtypes
nonsynonymousC<-subset(Virus, wtnt=='c' & TypeOfSite=='nonsyn')
nonsynonymousG<-subset(Virus, wtnt=='g' & TypeOfSite=='nonsyn')

#A and T synonymous and nonsynonymous graphs

#splits the graph into 4 sections, 2 rows 2 columns 
layout(matrix(c(1,2,3,4), nrow=2, byrow = TRUE))
#log meanfreq to get rid of 0's and make graph easier to read
plot(synonymousA$num,log10(synonymousA$MeanFreq+.000001), col = "red", 
     #pch creates the shape of the plotted points
     ,xlab="Position Number", ylab="MeanFreq",main="Synonymous A", pch=2)
#abline creates a line where the mean of all the points are located
abline(h=mean(log10(synonymousA$MeanFreq+.000001)))
#creates a legend:inset is coordinates where the legend lies

legend("topleft",inset=c(0,-0.45), legend=c
       ("synonymous A"),col="red", horiz=TRUE, 
#lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       lty=1, cex=0.8,xpd=TRUE,bty='n')
#plots continue 
plot(nonsynonymousA$num,log10(nonsynonymousA$MeanFreq+.000001), col = "blue",xlab="Position Number", ylab="MeanFreq",main="NonSynonymous A", pch = 5)
abline(h=mean(log10(nonsynonymousA$MeanFreq+.000001)))
legend("topleft",inset=c(0,-0.45), legend=c
       ("nonsynonymous A"),col="blue", horiz=TRUE, lty=1, cex=0.8,xpd=TRUE,bty='n')

plot(synonymousT$num,log10(synonymousT$MeanFreq+.000001), col = "green",xlab="Position Number", ylab="MeanFreq",main="Synonymous T")
abline(h=mean(log10(synonymousT$MeanFreq+.000001)))
legend("topleft",inset=c(0,-0.45), legend=c
       ("synonymous T"),col="green", horiz=TRUE, lty=1, cex=0.8,xpd=TRUE,bty='n')

plot(nonsynonymousT$num,log10(nonsynonymousT$MeanFreq+.000001), col = "purple",xlab="Position Number", ylab="MeanFreq",main="NonSynonymous T", pch=0)
abline(h=mean(log10(nonsynonymousT$MeanFreq+.000001)))
legend("topleft",inset=c(0,-0.45), legend=c
       ("nonsynonymous T"),col="purple", horiz=TRUE, lty=1, cex=0.8,xpd=TRUE,bty='n')
}
SynNonsynCG<-function(Virus){
#Plot C and G synonymous and nonsynonymous graphs
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
plot(synonymousC$num,log10(synonymousC$MeanFreq+.000001), col = "magenta",xlab="Position Number", ylab="MeanFreq",main="Synonymous C", pch=7)
abline(h=mean(log10(synonymousC$MeanFreq+.000001)))
legend("topleft",inset=c(0,-0.45), legend=c
       ("synonymous C"),col="magenta", horiz=TRUE, lty=1, cex=0.8,xpd=TRUE,bty='n')

plot(nonsynonymousC$num,log10(nonsynonymousC$MeanFreq+.000001), col = "cyan",xlab="Position Number", ylab="MeanFreq",main="NonSynonymous C", pch = 9)
abline(h=mean(log10(nonsynonymousC$MeanFreq+.000001)))
legend("topleft",inset=c(0,-0.45), legend=c
       ("nonsynonymous C"),col="cyan", horiz=TRUE, lty=1, cex=0.8,xpd=TRUE,bty='n')

plot(synonymousG$num,log10(synonymousG$MeanFreq+.000001), col = "black",xlab="Position Number", ylab="MeanFreq",main="Synonymous G",pch=13)
abline(h=mean(log10(synonymousG$MeanFreq+.000001)))
legend("topleft",inset=c(0,-0.45), legend=c
       ("synonymous G"),col="black", horiz=TRUE, lty=1, cex=0.8,xpd=TRUE,bty='n')

plot(nonsynonymousT$num,log10(nonsynonymousT$MeanFreq+.000001), col = "yellow",xlab="Position Number", ylab="MeanFreq",main="NonSynonymous G", pch=14)
abline(h=mean(log10(nonsynonymousT$MeanFreq+.000001)))
legend("topleft",inset=c(0,-0.45), legend=c
       ("nonsynonymous G"),col="yellow", horiz=TRUE, lty=1, cex=0.8,xpd=TRUE,bty='n')
}
SynNonsynAT(Dengue)
SynNonsynCG(Dengue)
