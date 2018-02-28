Dengue<-read.csv("DengueVirus1.fasta_pruned.mu.trim05_DF.csv")

synonymous<-which(Dengue$TypeOfSite=='syn')

synonymousA<-subset(Dengue, wtnt=='a' & TypeOfSite=='syn')
synonymousT<-subset(Dengue, wtnt=='t' & TypeOfSite=='syn') 
synonymous<-rbind(synonymousA, synonymousT)
nonsynonymous<-which(Dengue$TypeOfSite=='nonsyn')

nonsynonymousA<-subset(Dengue, wtnt=='a' & TypeOfSite=='nonsyn')
nonsynonymousT<-subset(Dengue, wtnt=='t' & TypeOfSite=='nonsyn')
nonsynonymous<-rbind(nonsynonymousA, nonsynonymousT)

synonymousC<-subset(Dengue, wtnt=='c' & TypeOfSite=='syn')
synonymousG<-subset(Dengue, wtnt=='g' & TypeOfSite=='syn')

nonsynonymousC<-subset(Dengue, wtnt=='c' & TypeOfSite=='nonsyn')
nonsynonymousG<-subset(Dengue, wtnt=='g' & TypeOfSite=='nonsyn')

# hist(Dengue$makesCpG[synonymous],prob=TRUE,breaks=20,main="Synoymous Histogram", xlab= "CpG",ylab= "Number of CpG Sites",col= "black")
# hist(Dengue$makesCpG[nonsynonymous],prob=TRUE,breaks=20,main="Non-Synoymous Histogram", xlab= "CpG",ylab= "Number of CpG Sites",col= "red")
# 
# 
# hist(Dengue$MeanFreq[synonymous],
#      prob=TRUE,breaks=20,main="Synoymous Nonsynoymous Histogram", 
#      xlab= "Mean Frequency",ylab= "Mutation Frequency",col= rgb(1,0,0,0.5))
# hist(Dengue$MeanFreq[nonsynonymous],
#      prob=TRUE,breaks=20,main="Non-Synoymous Histogram", 
#      xlab= "Mean Frequnecy",ylab= "Mutation Frequency",
#      col= rgb(0,0,1,0.5), add=TRUE)
# lines(density(Dengue$MeanFreq[synonymous]),col="red",lwd=2)
# lines(density(Dengue$MeanFreq[nonsynonymous]),col="blue",lwd=2)
# 
# 
# hist(nonsynonymousA$MeanFreq,
#      prob=TRUE,breaks=20,main="Synoymous Nonsynoymous Histogram", 
#      xlab= "Mean Frequency",ylab= "Mutation Frequency",col= rgb(1,0,0,0.5))
# hist(synonymousA$MeanFreq,
#      prob=TRUE,breaks=20,main="Non-Synoymous Histogram", 
#      xlab= "Mean Frequnecy",ylab= "Mutation Frequency",
#      col= rgb(0,0,1,0.5), add=TRUE)
# 
# hist(nonsynonymousT$MeanFreq,
#      prob=TRUE,breaks=20,main="Synoymous Nonsynoymous Histogram", 
#      xlab= "Mean Frequency",ylab= "Mutation Frequency",col= rgb(1,0,0,0.5))
# hist(synonymousT$MeanFreq,
#      prob=TRUE,breaks=20,main="Non-Synoymous Histogram", 
#      xlab= "Mean Frequnecy",ylab= "Mutation Frequency",
#      col= rgb(0,0,1,0.5), add=TRUE)
# plot(synonymousA$num,synonymousA$MeanFreq, col = "red")
# par(new=TRUE)
# plot(nonsynonymousA$num,nonsynonymousA$MeanFreq, col = "blue")
# par(new=TRUE)
# plot(synonymousT$num,synonymousT$MeanFreq, col = "green")
# par(new=TRUE)
# plot(nonsynonymousT$num,nonsynonymousT$MeanFreq, col = "purple")
# par(new=TRUE)
# plot(synonymous$num,synonymous$MeanFreq, col = "red")
# par(new=TRUE)
# plot(nonsynonymous$num,nonsynonymous$MeanFreq)

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
plot(synonymousA$num,log10(synonymousA$MeanFreq+.000001), col = "red",xaxt="n", yaxt="n",xlab="Position Number", ylab="MeanFreq",main="Synonymous A", pch=2)
abline(h=mean(log10(synonymousA$MeanFreq+.000001)))
legend("topleft",inset=c(0,-0.45), legend=c
       ("nonsynonymousT"),col="red", horiz=TRUE, lty=1:1, cex=0.8,xpd=TRUE,bty='n')

plot(nonsynonymousA$num,log10(nonsynonymousA$MeanFreq+.000001), col = "blue",xaxt="n", yaxt="n",xlab="Position Number", ylab="MeanFreq",main="NonSynonymous A", pch = 5)
abline(h=mean(log10(nonsynonymousA$MeanFreq+.000001)))
legend("topleft",inset=c(0,-0.45), legend=c
       ("nonsynonymousA"),col="blue", horiz=TRUE, lty=1:1, cex=0.8,xpd=TRUE,bty='n')

plot(synonymousT$num,log10(synonymousT$MeanFreq+.000001), col = "green",xaxt="n", yaxt="n",xlab="Position Number", ylab="MeanFreq",main="Synonymous T")
abline(h=mean(log10(synonymousT$MeanFreq+.000001)))
legend("topleft",inset=c(0,-0.45), legend=c
       ("synonymousT"),col="green", horiz=TRUE, lty=1:1, cex=0.8,xpd=TRUE,bty='n')

plot(nonsynonymousT$num,log10(nonsynonymousT$MeanFreq+.000001), col = "purple",xaxt="n", yaxt="n",xlab="Position Number", ylab="MeanFreq",main="NonSynonymous T", pch=0)
abline(h=mean(log10(nonsynonymousT$MeanFreq+.000001)))
legend("topleft",inset=c(0,-0.45), legend=c
       ("nonsynonymousT"),col="purple", horiz=TRUE, lty=1:1, cex=0.8,xpd=TRUE,bty='n')




# plot(synonymousA$num,log10(synonymousA$MeanFreq+.000001), 
#      col = "red",xaxt="n", yaxt="n",xlab="", ylab="", pch=2)
# abline(h=mean(log10(synonymousA$MeanFreq+.000001)))
# par(new=TRUE)
# plot(nonsynonymousA$num,log10(nonsynonymousA$MeanFreq+.000001), col = "blue",xaxt="n", yaxt="n",xlab="Position Number", ylab="MeanFreq",main="Synonymous and Nonsynonymous A", pch=5)
# abline(h=mean(log10(nonsynonymousA$MeanFreq+.000001)))
# legend("topleft",inset=c(0,-0.45), legend=c("synonymousA", "nonsynonymousA"),col=c("red", "blue"), horiz=TRUE, lty=1:1, cex=0.8,xpd=TRUE,bty='n')
# plot(synonymousT$num,synonymousT$MeanFreq, col = "green",,xaxt="n", yaxt="n",xlab="", ylab="")
# par(new=TRUE)
# plot(nonsynonymousT$num,nonsynonymousT$MeanFreq, col = "purple",xaxt="n", yaxt="n",xlab="Position Number", ylab="MeanFreq",main="Synonymous and Nonsynonymous T", pch=0)
# legend("topleft",inset=c(0,-0.45), legend=c("synonymousT", "nonsynonymousT"),col=c("green", "purple"), horiz=TRUE, lty=1:1, cex=0.8,xpd=TRUE,bty='n')
# 
# syn_Cpg<-Dengue$MeanFreq[synonymous]
# nonsyn_Cpg<-Dengue$MeanFreq[nonsynonymous]
# wilcox.test(syn_Cpg,nonsyn_Cpg)

