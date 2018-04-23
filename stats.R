
#Dvirus<-read.csv("DengueVirus1.fasta_pruned.mu.trim05_DF.csv", header = T)
#data<-read.csv("DengueVirus1.fasta_pruned.mu.trim05_DF.csv", header = T)
#place things in table
Wilcox_test = function(data){
library(graphics)
library(dplyr)
library(plyr)


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
print("For a: Comparing CpG with noCpG (nonsyn). Wilcox test less: yellow/green")
print(wilcox.test(array3, array4, alternative='less'))
# print("For a: Comparing makes syn mutation to nonsyn mutation. Both makes CpG, Wilcox test greater")
# print(wilcox.test(array1, array3, alternative='greater'))
# print("For a: Comparing  makes syn mutation to nonsyn mutation. No CpG made, Wilcox test greater")
# print(wilcox.test(array2, array4, alternative='greater'))
print("For a: Comparing  syn to nonsyn. Wilcox test greater red&blue vs yellow&green")
print(wilcox.test(syna, nonsyna, alternative='greater'))
# print("For a: Comparing  CpG to no CpG. Wilcox test less")
# print(wilcox.test(CpGa, nonCpGa, alternative='less'))


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
print("For t: Comparing CpG with noCpG (nonsyn). Wilcox test less: yellow/green")
print(wilcox.test(array7, array8, alternative='less'))
# print("For t: Comparing  makes syn mutation to nonsyn mutation. Both makes CpG, Wilcox test greater")
# print(wilcox.test(array5, array7, alternative='greater'))
# print("For t: Comparing  makes syn mutation to nonsyn mutation. No CpG made, Wilcox test greater")
# print(wilcox.test(array5, array8, alternative='greater'))
print("For t: Comparing  syn to nonsyn. Wilcox test greater red&blue vs yellow&green")
print(wilcox.test(synt, nonsynt, alternative='greater'))
# print("For t: Comparing  CpG to no CpG. Wilcox test less")
# print(wilcox.test(CpGt, nonCpGt, alternative='less'))


# sync = data$MeanFreq[data$wtnt =="c" & data$TypeOfSite == 'syn']
# nonsync = data$MeanFreq[data$wtnt =="c" & data$TypeOfSite == 'nonsyn']
# print("For c: Comparing  syn to nonsyn. Wilcox test greater")
# print(wilcox.test(sync, nonsync, alternative='greater'))
# 
# syng = data$MeanFreq[data$wtnt =="g" & data$TypeOfSite == 'syn']
# nonsyng = data$MeanFreq[data$wtnt =="g" & data$TypeOfSite == 'nonsyn']
# print("For g: Comparing  syn to nonsyn. Wilcox test greater")
# print(wilcox.test(syng, nonsyng, alternative='greater'))

}
#Wilcox_test(data)

# low <- c(1,1,2,2,3,3,4,4,5,5)
# high <- c(6,6,7,7,8,8,9,9,10,10, 12, 78)
# 
# wilcox.test(low, high, alternative='less')
# wilcox.test(low, high, alternative='greater')
# wilcox.test(low, high)
# # if p<.05 is significant
# num = 3
# low2= rep(1,num)
# high2 = rep (3, num)
# wilcox.test(low2, high2, alternative='less')








