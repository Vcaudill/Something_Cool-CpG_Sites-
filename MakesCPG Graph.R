Dengue<-read.csv("DengueVirus1.fasta_pruned.mu.trim05_DF.csv")

CPGNoCPGAT<-function(Virus){
#subset wildtype and if a CPG site is made
CPGA<-subset(Virus, wtnt=='a' & makesCpG==1)
CPGT<-subset(Virus, wtnt=='t' & makesCpG==1) 

#subsets wildtype and if no CPG sites are made
NoCPGA<-subset(Virus, wtnt=='a' & makesCpG==0)
NoCPGT<-subset(Virus, wtnt=='t' & makesCpG==0)

#layout graphs take same amount of space and 2 columns
layout(matrix(c(1,2),ncol = 2,byrow = TRUE))
#plot position number vs mean freq
plot(CPGA$num,log10(CPGA$MeanFreq+.000001), col = "red",xlab="Position Number", ylab="MeanFreq",main="CPG & No CPG A", pch=2)
#creates a line with mean of all points
abline(h=mean(log10(CPGA$MeanFreq+.000001)),col = "red")
#creates legend with red for CPG A and blue for No CPG A
legend("topleft",inset=c(0,1.01), legend=c
       ("CPG A","No CPG A"),col=c("red","blue"), horiz=TRUE, lty=1, cex=0.8,xpd=TRUE,bty='n')
#overlaps second graph on top of first graph
par(new=TRUE)
#plots second graph No CPG A. xaxt gets rid of x axis label yaxt gets rid of y
plot(NoCPGA$num,log10(NoCPGA$MeanFreq+.000001), col=c("red","blue"),xaxt="n", yaxt="n",xlab="Position Number", ylab="MeanFreq", pch=5)
abline(h=mean(log10(NoCPGA$MeanFreq+.000001)),col="blue")

plot(CPGT$num,log10(CPGT$MeanFreq+.000001), col = "green",xlab="Position Number", ylab="MeanFreq",main="CPG & No CPG T")
abline(h=mean(log10(CPGT$MeanFreq+.000001)),col = "green")
legend("topleft",inset=c(0,1.01), legend=c
       ("CPG T","No CPG T"),col=c("green","purple"), horiz=TRUE, lty=1, cex=0.8,xpd=TRUE,bty='n')
par(new=TRUE)
plot(NoCPGT$num,log10(NoCPGT$MeanFreq+.000001), col="purple",xaxt="n", yaxt="n",xlab="Position Number", ylab="MeanFreq", pch=0)
abline(h=mean(log10(NoCPGT$MeanFreq+.000001)),col="purple")
}

CPGNoCPGAT(Dengue)
