#RYAN WILcOXTEST FILE
#Dvirus<-read.csv("DengueVirus1.fasta_pruned.mu.trim05_DF.csv", header = T)
#data<-read.csv("DengueVirus1.fasta_pruned.mu.trim05_DF.csv", header = T)
#place things in table

Wilcox_test = function(data, truename){
  truenamepdf= paste(truename, ".pdf", sep="")
  library(graphics)
  library(dplyr)
  library(plyr)
  pVals = c()
  shrtval = 0
  options(scipen=999)
  
  array1 = data$MeanFreq[data$wtnt =="a" & data$TypeOfSite == 'syn' & data$makesCpG == 1]
  array2 = data$MeanFreq[data$wtnt =="a" & data$TypeOfSite == 'syn' & data$makesCpG == 0]
  array3 = data$MeanFreq[data$wtnt =="a" & data$TypeOfSite == 'nonsyn' & data$makesCpG == 1]
  array4 = data$MeanFreq[data$wtnt =="a" & data$TypeOfSite == 'nonsyn' & data$makesCpG == 0]
  syna = data$MeanFreq[data$wtnt =="a" & data$TypeOfSite == 'syn']
  nonsyna = data$MeanFreq[data$wtnt =="a" & data$TypeOfSite == 'nonsyn']
  CpGa = data$MeanFreq[data$wtnt =="a" & data$makesCpG == 1]
  nonCpGa = data$MeanFreq[data$wtnt =="a" &  data$makesCpG == 0]
  
  print("For a: Comparing makes CpG with noCpG (syn). Wilcox test less: red/blue")
  print(wilcox.test(array1, array2, alternative='less'))
  pVals = c(pVals,format(wilcox.test(array1, array2, alternative='less')$p.value, nsmall = 6))
  print(pVals)
  print("For a: Comparing CpG with noCpG (nonsyn). Wilcox test less: yellow/green")
  print(wilcox.test(array3, array4, alternative='less'))
  pVals = c(pVals,format(wilcox.test(array3, array4, alternative='less')$p.value, nsmall = 6))
  print(pVals)
  print("For a: Comparing  syn to nonsyn. Wilcox test greater red&blue vs yellow&green")
  print(wilcox.test(syna, nonsyna, alternative='greater'))
  print(wilcox.test(syna, nonsyna, alternative='greater')$p.value)
  pVals = c(pVals,format(wilcox.test(syna, nonsyna, alternative='greater')$p.value, nsmall = 6))
  print(pVals)
  
  array5 = data$MeanFreq[data$wtnt =="t" & data$TypeOfSite == 'syn' & data$makesCpG == 1]
  array6 = data$MeanFreq[data$wtnt =="t" & data$TypeOfSite == 'syn' & data$makesCpG == 0]
  array7 = data$MeanFreq[data$wtnt =="t" & data$TypeOfSite == 'nonsyn' & data$makesCpG == 1]
  array8 = data$MeanFreq[data$wtnt =="t" & data$TypeOfSite == 'nonsyn' & data$makesCpG == 0]
  synt = data$MeanFreq[data$wtnt =="t" & data$TypeOfSite == 'syn']
  nonsynt = data$MeanFreq[data$wtnt =="t" & data$TypeOfSite == 'nonsyn']
  CpGt = data$MeanFreq[data$wtnt =="t" & data$makesCpG == 1]
  nonCpGt = data$MeanFreq[data$wtnt =="t" &  data$makesCpG == 0]
  
  print("For t: Comparing makes CpG with noCpG (syn). Wilcox test less: red/blue")
  print(wilcox.test(array5, array6, alternative='less'))
  pVals = c(pVals,format(wilcox.test(array5, array6, alternative='less')$p.value, nsmall = 6))
  print("For t: Comparing CpG with noCpG (nonsyn). Wilcox test less: yellow/green")
  print(wilcox.test(array7, array8, alternative='less'))
  pVals = c(pVals,format(wilcox.test(array7, array8, alternative='less')$p.value, nsmall = 6))
  print("For t: Comparing  syn to nonsyn. Wilcox test greater red&blue vs yellow&green")
  print(wilcox.test(synt, nonsynt, alternative='greater'))
  pVals = c(pVals,format(wilcox.test(synt, nonsynt, alternative='greater')$p.value, nsmall = 6))
  
  Pvalues= c(pVals)
  setwd("~/Desktop/Something_Cool-CpG_Sites-")
  source("maketable.R")
  makeTable(Pvalues, truenamepdf, truename)
  
}







