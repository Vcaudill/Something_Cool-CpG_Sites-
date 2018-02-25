
Dvirus<-read.csv("DengueVirus1.fasta_pruned.mu.trim05_DF.csv", header = T)
Wilcox_test = function(data){
library(graphics)
library(dplyr)
library(plyr)


array1 = data$MeanFreq[data$wtnt =="a" & data$TypeOfSite == 'syn' & data$makesCpG == 1]
array2 = data$MeanFreq[data$wtnt =="a" & data$TypeOfSite == 'syn' & data$makesCpG == 0]
array3 = data$MeanFreq[data$wtnt =="a" & data$TypeOfSite == 'nonsyn' & data$makesCpG == 1]
array4 = data$MeanFreq[data$wtnt =="a" & data$TypeOfSite == 'nonsyn' & data$makesCpG == 0]


print(wilcox.test(array1, array2, alternative='greater'))
print(wilcox.test(array1, array3, alternative='greater'))
print(wilcox.test(array2, array4, alternative='greater'))


array5 = data$MeanFreq[data$wtnt =="t" & data$TypeOfSite == 'syn' & data$makesCpG == 1]
array6 = data$MeanFreq[data$wtnt =="t" & data$TypeOfSite == 'syn' & data$makesCpG == 0]
array7 = data$MeanFreq[data$wtnt =="t" & data$TypeOfSite == 'nonsyn' & data$makesCpG == 1]
array8 = data$MeanFreq[data$wtnt =="t" & data$TypeOfSite == 'nonsyn' & data$makesCpG == 0]

print(wilcox.test(array5, array6, alternative='greater'))
print(wilcox.test(array5, array7, alternative='greater'))
print(wilcox.test(array5, array8, alternative='greater'))


}
Wilcox_test(data)

# #subset into two groups yes makes cpg and no cpg
# cpg.y<-subset(Dvirus, makesCpG==1)
# cpg.n<-subset(Dvirus, makesCpG==0)
# #subset further into letters 
# AC<-subset(cpg.y, wtnt=='a')
# ANC<-subset(cpg.n, wtnt=='a')
# GC<-subset(cpg.y, wtnt=='g')
# GNC<-subset(cpg.n, wtnt=='g')
# TC<-subset(cpg.y, wtnt=='t')
# TNC<-subset(cpg.n, wtnt=='t')
# CC<-subset(cpg.y, wtnt=='c')
# CNC<-subset(cpg.n, wtnt=='c')
# 
# #for wilcox test
# CpGdataA<-data.frame(CpG=c(rep("Y",times=length(AC$MeanFreq)),rep("N",times=length(ANC$MeanFreq))), MeanFreq=c(AC$MeanFreq,ANC$MeanFreq))
# CpGdataANC<-data.frame(CpG=c(rep("Y",times=length(ANC$MeanFreq)),rep("N",times=length(AC$MeanFreq))), MeanFreq=c(ANC$MeanFreq,AC$MeanFreq))
# 
# CpGdataG<-data.frame(CpG=c(rep("Y",times=length(GC$MeanFreq)),rep("N",times=length(GNC$MeanFreq))), MeanFreq=c(GC$MeanFreq,GNC$MeanFreq))
# CpGdataGNC<-data.frame(CpG=c(rep("Y",times=length(GNC$MeanFreq)),rep("N",times=length(GC$MeanFreq))), MeanFreq=c(GNC$MeanFreq,GC$MeanFreq))
# 
# CpGdataT<-data.frame(CpG=c(rep("Y",times=length(TC$MeanFreq)),rep("N",times=length(TNC$MeanFreq))), MeanFreq=c(TC$MeanFreq,TNC$MeanFreq))
# CpGdataTNC<-data.frame(CpG=c(rep("Y",times=length(TNC$MeanFreq)),rep("N",times=length(TC$MeanFreq))), MeanFreq=c(TNC$MeanFreq,TC$MeanFreq))
# 
# CpGdataC<-data.frame(CpG=c(rep("Y",times=length(CC$MeanFreq)),rep("N",times=length(CNC$MeanFreq))), MeanFreq=c(CC$MeanFreq,CNC$MeanFreq))
# CpGdataCNC<-data.frame(CpG=c(rep("Y",times=length(CNC$MeanFreq)),rep("N",times=length(CC$MeanFreq))), MeanFreq=c(CNC$MeanFreq,CC$MeanFreq))
# 
# 
# # test for A and T what are the differences? 
# # will not work if everthing is no or yes, must have two levels?
# wilcox.test(MeanFreq ~ CpG, data=CpGdataA, alternative='greater', paired = FALSE)
# wilcox.test(MeanFreq ~ CpG, data=CpGdataANC, alternative='greater', paired = FALSE)
# 
# wilcox.test(MeanFreq ~ CpG, data=CpGdataT, alternative='greater', paired = FALSE)
# wilcox.test(MeanFreq ~ CpG, data=CpGdataTNC, alternative='greater', paired = FALSE)














