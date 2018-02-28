Dengue<-read.csv("DengueVirus1.fasta_pruned.mu.trim05_DF.csv")


CPGA<-subset(Dengue, wtnt=='a' & makesCpG==1)
CPGT<-subset(Dengue, wtnt=='t' & makesCpG==1) 


NoCPGA<-subset(Dengue, wtnt=='a' & makesCpG==0)
NoCPGT<-subset(Dengue, wtnt=='t' & makesCpG==0)



layout(matrix(c(1,2),2,2,byrow = TRUE))
plot(CPGA$num,log10(CPGA$MeanFreq+.000001), col = "red",xaxt="n", yaxt="n",xlab="Position Number", ylab="MeanFreq",main="CPG & No CPG A", pch=2)
abline(h=mean(log10(CPGA$MeanFreq+.000001)),col = "red")
legend("topleft",inset=c(0,1), legend=c
       ("CPG A","No CPG A"),col=c("red","blue"), horiz=TRUE, lty=1:1, cex=0.8,xpd=TRUE,bty='n')
par(new=TRUE)
plot(NoCPGA$num,log10(NoCPGA$MeanFreq+.000001), col=c("red","blue"),xaxt="n", yaxt="n",xlab="Position Number", ylab="MeanFreq", pch=5)
abline(h=mean(log10(NoCPGA$MeanFreq+.000001)),col="blue")

plot(CPGT$num,log10(CPGT$MeanFreq+.000001), col = "green",xaxt="n", yaxt="n",xlab="Position Number", ylab="MeanFreq",main="CPG & No CPG T")
abline(h=mean(log10(CPGT$MeanFreq+.000001)),col = "green")
legend("topleft",inset=c(0,1), legend=c
       ("CPG T","No CPG T"),col=c("green","purple"), horiz=TRUE, lty=1:1, cex=0.8,xpd=TRUE,bty='n')
par(new=TRUE)
plot(NoCPGT$num,log10(NoCPGT$MeanFreq+.000001), col="purple",xaxt="n", yaxt="n",xlab="Position Number", ylab="MeanFreq", pch=0)
abline(h=mean(log10(NoCPGT$MeanFreq+.000001)),col="purple")