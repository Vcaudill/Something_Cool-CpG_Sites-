virusname = 'DengueVirus1.fasta_pruned.mu.trim05.txt'
 
splitname<-unlist(strsplit(virusname,".fasta"))
truename<-splitname[1]

truenameRda= paste(truename, ".Rda", sep="")
load(truenameRda)
DengueVirus1 = DF
DengueVirusTest1= DF
DengueVirusTest2= DF
DengueVirusTest3= DF
DengueVirusTest4= DF

DengueVirusTest1$Virus<-('DVT1')
DengueVirusTest2$Virus<-('DVT2')
DengueVirusTest3$Virus<-('DVT3')
DengueVirusTest4$Virus<-('DVT4')

my.list <- list(DengueVirus1, DengueVirusTest1, DengueVirusTest2, DengueVirusTest3, DengueVirusTest4)
data_points = data.frame("Count"= 1:length(my.list), "Virus"= 1:length(my.list))
count = 1

for (data in my.list){  
  cpg.y<-subset(data, makesCpG==1)
  cpg.n<-subset(data, makesCpG==0)
  #subset further into letters nuclotideCpgforming or nucotideNonGpg
  AC<-subset(cpg.y, wtnt=='a')
  ANC<-subset(cpg.n, wtnt=='a')
  TC<-subset(cpg.y, wtnt=='t')
  TNC<-subset(cpg.n, wtnt=='t')
 
  
  #Function to help create errorbars
  sem<-function(x){
    return(sd(x,na.rm = FALSE)/sqrt(length(x)))
  }
  
  #making the data frames with all information about a, t, c, g 
  AllA = rbind(AC, ANC)
  AllT = rbind(TC, TNC)
 
  # for loops to caculate mean ans errorbars
  for (i in 1:length(AllA$makesCpG)) {
    if (AllA$makesCpG[i] == 1 && AllA$TypeOfSite[i] == "syn") {
      AllA_mean_value_syn_CpG <- mean(AllA$MeanFreq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") )])
      AllA_sem_vals_syn_CpG<-sem(AllA$MeanFreq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") )])
    }
    if (AllA$makesCpG[i] == 1 && AllA$TypeOfSite[i] == "nonsyn") {
      AllA_mean_value_nonsyn_CpG<- mean(AllA$MeanFreq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "nonsyn") )])
      AllA_sem_vals_nonsyn_CpG<-sem(AllA$MeanFreq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "nonsyn") )])
    }
    if (AllA$makesCpG[i] == 0 && AllA$TypeOfSite[i] == "syn") {
      AllA_mean_value_syn_nCpG <- mean(AllA$MeanFreq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") )])
      AllA_sem_vals_syn_nCpG<-sem(AllA$MeanFreq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") )])
    }
    if (AllA$makesCpG[i] == 0 && AllA$TypeOfSite[i] == "nonsyn") {
      AllA_mean_value_nonsyn_nCpG<- mean(AllA$MeanFreq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "nonsyn") )])
      AllA_sem_vals_nonsyn_nCpG<-sem(AllA$MeanFreq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "nonsyn") )])
    }
  }
  
  for (i in 1:length(AllT$makesCpG)) {
    if (AllT$makesCpG[i] == 1 && AllT$TypeOfSite[i] == "syn") {
      AllT_mean_value_syn_CpG <- mean(AllT$MeanFreq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )])
      AllT_sem_vals_syn_CpG<-sem(AllT$MeanFreq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )])
    }
    if (AllT$makesCpG[i] == 1 && AllT$TypeOfSite[i] == "nonsyn") {
      AllT_mean_value_nonsyn_CpG <- mean(AllT$MeanFreq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "nonsyn") )])
      AllT_sem_vals_nonsyn_CpG<-sem(AllT$MeanFreq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "nonsyn") )])
    }
    if (AllT$makesCpG[i] == 0 && AllT$TypeOfSite[i] == "syn") {
      AllT_mean_value_syn_nCpG<- mean(AllT$MeanFreq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )])
      AllT_sem_vals_syn_nCpG<-sem(AllT$MeanFreq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )])
    }
    if (AllT$makesCpG[i] == 0 && AllT$TypeOfSite[i] == "nonsyn") {
      AllT_mean_value_nonsyn_nCpG <- mean(AllT$MeanFreq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "nonsyn") )])
      AllT_sem_vals_nonsyn_nCpG<-sem(AllT$MeanFreq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "nonsyn") )])
    }
  }
  
  # There are the upper and lower limits of the error bar
  # AllA_LCLS = AllA$mean_value - AllA$sem_vals
  # AllA_UCLS = AllA$mean_value + AllA$sem_vals
  # 
  # AllT_LCLS = AllT$mean_value - AllT$sem_vals
  # AllT_UCLS = AllT$mean_value + AllT$sem_vals
  # 
  
  data_points$Virus[count]= data$Virus[1]
  data_points$AsynNC_C[count]= AllA_mean_value_syn_nCpG/AllA_mean_value_syn_CpG
  data_points$AnonsynNC_C[count]= AllA_mean_value_nonsyn_nCpG/AllA_mean_value_nonsyn_CpG
  
  data_points$TsynNC_C[count] = AllT_mean_value_syn_nCpG/AllT_mean_value_syn_CpG
  data_points$meanTnonsynNC_C[count] = AllT_mean_value_nonsyn_nCpG/AllT_mean_value_nonsyn_CpG
 

  data_points$AsynC_LCLS[count]= AllA_mean_value_syn_CpG - AllA_sem_vals_syn_CpG
  data_points$AnonsynC_LCLS[count]= AllA_mean_value_nonsyn_CpG - AllA_sem_vals_nonsyn_CpG
  data_points$AsynNC_LCLS[count]= AllA_mean_value_syn_nCpG - AllA_sem_vals_syn_nCpG
  data_points$AnonsynNC_LCLS[count]= AllA_mean_value_nonsyn_nCpG - AllA_sem_vals_nonsyn_nCpG
  data_points$TsynC_LCLS[count] = AllT_mean_value_syn_CpG - AllT_sem_vals_syn_CpG
  data_points$TnonsynC_LCLS[count] =AllT_mean_value_nonsyn_CpG - AllT_sem_vals_nonsyn_CpG
  data_points$TsynNC_LCLS[count] = AllT_mean_value_syn_nCpG - AllT_sem_vals_syn_nCpG 
  data_points$TnonsynNC_LCLS[count] =AllT_mean_value_nonsyn_nCpG - AllT_sem_vals_nonsyn_nCpG
  
  data_points$AsynC_UCLS[count]= AllA_mean_value_syn_CpG + AllA_sem_vals_syn_CpG
  data_points$AnonsynC_UCLS[count]= AllA_mean_value_nonsyn_CpG + AllA_sem_vals_nonsyn_CpG
  data_points$AsynNC_UCLS[count]= AllA_mean_value_syn_nCpG + AllA_sem_vals_syn_nCpG
  data_points$AnonsynNC_UCLS[count]= AllA_mean_value_nonsyn_nCpG + AllA_sem_vals_nonsyn_nCpG
  data_points$TsynC_UCLS[count] = AllT_mean_value_syn_CpG + AllT_sem_vals_syn_CpG
  data_points$TnonsynC_UCLS[count] =AllT_mean_value_nonsyn_CpG + AllT_sem_vals_nonsyn_CpG
  data_points$TsynNC_UCLS[count] = AllT_mean_value_syn_nCpG + AllT_sem_vals_syn_nCpG 
  data_points$TnonsynNC_UCLS[count] =AllT_mean_value_nonsyn_nCpG + AllT_sem_vals_nonsyn_nCpG
  count = count +1
}

# graphing 

plot(data_points$Count-.3, data_points$AsynNC_C, main="Scatterplot Example", 
     xlab="Virus ", ylab="Rate", pch=19, col= "red", log = 'y', yaxt="n", xaxt="n", ylim=c(1,10), xlim=c(.5, length(my.list) +1.5))
aty <- axTicks(2)
labels <- sapply(aty,function(i)
  as.expression(bquote(10^ .(i)))
)
axis(2,at=aty,labels=labels)
points(data_points$Count -.1, data_points$AnonsynNC_C, col= "green", pch=19)
points(data_points$Count+.1, data_points$TsynNC_C, col= "blue", pch=19)
points(data_points$Count+.3, data_points$TnonsynNC_C, col= "purple", pch=19)
par(new=TRUE)

# hack: we draw arrows but with very special "arrowheads"
arrows(data_points$Count-.3, data_points$AsynNC_LCLS/data_points$AsynC_LCLS, data_points$Count -.3, data_points$AsynNC_UCLS/data_points$AsynC_UCLS, length=0.05, angle=90, code=3, col= "red")
arrows(data_points$Count-.1, data_points$AnonsynNC_LCLS/data_points$AnonsynC_LCLS, data_points$Count -.1, data_points$AnonsynNC_UCLS/data_points$AnonsynC_UCLS, length=0.05, angle=90, code=3, col= "green")
arrows(data_points$Count+.1, data_points$TsynNC_LCLS/data_points$TsynC_LCLS, data_points$Count+.1, data_points$TsynNC_UCLS/data_points$TsynC_UCLS, length=0.05, angle=90, code=3, col= "blue")
arrows(data_points$Count+.3, data_points$TnonsynNC_LCLS/data_points$TnonsynC_LCLS, data_points$Count +.3, data_points$TnonsynNC_UCLS/data_points$TnonsynC_UCLS, length=0.05, angle=90, code=3, col= "purple")


axis(1, at=1:length(my.list), labels=data_points$Virus, las= 2)
legend((length(my.list) + .4), 7, legend=c("A Syn", "A NonSyn", "T Syn", "T NonSyn"),
       col=c("red", "green", "blue", "purple"), lty=1, cex=1)





