
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/Rda_Files")

library(scales)
library(plotrix)
library(sfsmisc)

DengueVirus1 <- read.csv("data/csv/DengueVirus1.CSV")
BK_polyomavirus_VP1<-read.csv("data/csv/BKpolyomavirus_VP1.CSV")

load("DengueVirus2.Rda")
DengueVirus2= DF
load("DengueVirus3.Rda")
DengueVirus3= DF
load("DengueVirus4.Rda")
DengueVirus4= DF

load("EnterovirusA_VP1.Rda")
EnterovirusA_VP1= DF
load("EnterovirusA_VP2.Rda")
EnterovirusA_VP2= DF
load("EnterovirusB_VP1.Rda")
EnterovirusB_VP1= DF
load("EnterovirusB_VP2.Rda")
EnterovirusB_VP2= DF

load("EnterovirusC_VP1.Rda")
EnterovirusC_VP1= DF
load("EnterovirusC_VP2.Rda")
EnterovirusC_VP2= DF
load("EnterovirusD_VP1.Rda")
EnterovirusD_VP1= DF
load("HumanBocavirus1_NS1.Rda")
HumanBocavirus1_NS1= DF
load("HumanBocavirus1_VP1.Rda")
HumanBocavirus1_VP1= DF

load("humanparainfluenzavirus1_F.Rda")
humanparainfluenzavirus1_F= DF
load("humanparainfluenzavirus1_HN.Rda")
humanparainfluenzavirus1_HN= DF
load("humanparainfluenzavirus3_HN.Rda")
humanparainfluenzavirus3_HN= DF
load("InfluenzaAvirus_HA_H1N1.Rda")
InfluenzaAvirus_HA_H1N1= DF

load("InfluenzaAvirus_HA_H3N2.Rda")
InfluenzaAvirus_HA_H3N2= DF
load("InfluenzaAvirus_NA_H1N1.Rda")
InfluenzaAvirus_NA_H1N1= DF
load("InfluenzaAvirus_NA_H3N2.Rda")
InfluenzaAvirus_NA_H3N2= DF
load("InfluenzaBvirus_NA.Rda")
InfluenzaBvirus_NA= DF

load("InfluenzaBvirus_HA.Rda")
InfluenzaBvirus_HA= DF



DengueVirus1$Virus<-('DengueVirus1')
BK_polyomavirus_VP1$Virus<-('BK_polyomavirus_VP1')
DengueVirus3$Virus<-('DengueVirus3')
DengueVirus4$Virus<-('DengueVirus4')
DengueVirus2$Virus<-('DengueVirus2')

EnterovirusA_VP1$Virus<-('EnterovirusA_VP1')
EnterovirusA_VP2$Virus<-('EnterovirusA_VP2')
EnterovirusB_VP1$Virus<-('EnterovirusB_VP1')
EnterovirusB_VP2$Virus<-('EnterovirusB_VP2')

EnterovirusC_VP1$Virus<-('EnterovirusC_VP1')
EnterovirusD_VP1$Virus<-('EnterovirusD_VP1')
HumanBocavirus1_NS1$Virus<-('HumanBocavirus1_NS1')
HumanBocavirus1_VP1$Virus<-('HumanBocavirus1_VP1')

humanparainfluenzavirus1_F$Virus<-('humanparainfluenzavirus1_F')
humanparainfluenzavirus1_HN$Virus<-('humanparainfluenzavirus1_HN')
humanparainfluenzavirus3_HN$Virus<-('humanparainfluenzavirus3_HN')
InfluenzaAvirus_HA_H1N1$Virus<-('InfluenzaAvirus_HA_H1N1')

InfluenzaAvirus_HA_H3N2$Virus<-('InfluenzaAvirus_HA_H3N2')
InfluenzaAvirus_NA_H1N1$Virus<-('InfluenzaAvirus_NA_H1N1')
InfluenzaAvirus_NA_H3N2$Virus<-('InfluenzaAvirus_NA_H3N2')
InfluenzaBvirus_HA$Virus<-('InfluenzaBvirus_HA')

InfluenzaBvirus_NA$Virus<-('InfluenzaBvirus_NA')


my.list <- list(DengueVirus1, DengueVirus2, DengueVirus3, DengueVirus4, humanparainfluenzavirus1_F, humanparainfluenzavirus1_HN, humanparainfluenzavirus3_HN, InfluenzaAvirus_HA_H1N1,InfluenzaAvirus_HA_H3N2, InfluenzaAvirus_NA_H1N1, InfluenzaAvirus_NA_H3N2,InfluenzaBvirus_HA, InfluenzaBvirus_NA, EnterovirusA_VP1, EnterovirusA_VP2,EnterovirusB_VP1, EnterovirusB_VP2,EnterovirusC_VP1,EnterovirusC_VP2,EnterovirusD_VP1, BK_polyomavirus_VP1, HumanBocavirus1_NS1, HumanBocavirus1_VP1)
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
  data_points$TnonsynNC_C[count] = AllT_mean_value_nonsyn_nCpG/AllT_mean_value_nonsyn_CpG
 

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

#making points that are Nah or infint or 0 Not costly on the graph
for (i in row(data_points))
  for (j in 1:6){
    if(data_points[i,j] == "NaN")
      #no nonCpG or CpG mutations
      data_points[i,j] = 0.002
    if(data_points[i,j] == "Inf")
     # no CpG mutations
      data_points[i,j] = 700
    if(data_points[i,j] == "0")
      #0 no nonCpG mutations
      data_points[i,j] = .0005
}
#one of the error bars is too large at 2.47 e13 so we are placing it lower, but noting it high amount  
data_points[20,12] = data_points[20,14]/300 


print(data_points$TnonsynNC_LCLS/data_points$TnonsynC_LCLS)
# graphing 
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")
png("Costly_Graph_13.png", width = 15, height = 8, units = "in", res= 500)
par(mar=c(5,3,3,1), oma=c(10,4,1,1))
plot(data_points$Count-.3, data_points$AsynNC_C, main="How Costly is A CpG Mutation?", xlab=" ", yaxt = "n",
     ylab="Costly", pch=19, col= "red", log = 'y', xaxt = "n", ylim=c(0.0005, 700), xlim=c(1.5, length(my.list) +3.5), las= 1, cex.main=3)
# aty <- axTicks(2)
# labels <- sapply(aty,function(i)
#   as.expression(bquote(10^ .(i)))
# )
# axis(2,at=aty,labels=labels)
points(data_points$Count -.1, data_points$AnonsynNC_C, col= "green", pch=19)
points(data_points$Count+.1, data_points$TsynNC_C, col= "blue", pch=19)
points(data_points$Count + .3, data_points$TnonsynNC_C, col= "purple", pch=19)
# par(new=TRUE)

# hack: we draw arrows but with very special "arrowheads"
arrows(data_points$Count-.3, data_points$AsynNC_LCLS/data_points$AsynC_LCLS, data_points$Count -.3, data_points$AsynNC_UCLS/data_points$AsynC_UCLS, length=0.05, angle=90, code=3, col= "red")
arrows(data_points$Count-.1, data_points$AnonsynNC_LCLS/data_points$AnonsynC_LCLS, data_points$Count -.1, data_points$AnonsynNC_UCLS/data_points$AnonsynC_UCLS, length=0.05, angle=90, code=3, col= "green")
arrows(data_points$Count+.1, data_points$TsynNC_LCLS/data_points$TsynC_LCLS, data_points$Count+.1, data_points$TsynNC_UCLS/data_points$TsynC_UCLS, length=0.05, angle=90, code=3, col= "blue")
arrows(data_points$Count+.3, data_points$TnonsynNC_LCLS/data_points$TnonsynC_LCLS, data_points$Count +.3, data_points$TnonsynNC_UCLS/data_points$TnonsynC_UCLS, length=0.05, angle=90, code=3, col= "purple")

axis(2, at = c(.01, 0.5,1,2,5,10,20,50,100), labels = c(0.01, 0.5,1,2,5,10,20,50,100),  las=2)
axis.break(2, 0.007,breakcol="black",style="slash")
mtext('No nonCpG \n  or  CpG     \n    mutations  ', side=2, line=.005, at=0.002, las=1.1, cex = .7)
mtext('No nonCpG \n   mutations   ', side=2, line=.005, at=0.0005, las=1.1, cex = .7)

axis.break(2, 200,breakcol="black",style="slash")
mtext('2.47e+13', side=2, line=.005, at=300, las=1.1, cex = .9)

mtext('xTimes as Costly', side=2, line=5, at=.5, las=0, cex = 2)

mtext('No CpG \n mutations ', side=2, line=.005, at=700, las=1.1, cex = .7)

abline(h=c(0.002, 0.0005, .01,0.5,1,2,5,10,20,50,100, 300, 700), col="grey", lty=c(2,2))
abline(v=c(1.5,2.5,3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,12.5,13.5,14.5,16.5,17.5,18.5,19.5,15.5, 20.5, 21.5, 22.5, 23.5), col="grey", lty=c(1))


# xlab="Virus "

axis(1, at=1:length(my.list), labels=data_points$Virus, las= 2)
legend((length(my.list) + 1.4), 10, legend=c("A Syn", "A NonSyn", "T Syn", "T NonSyn"),
       col=c("red", "green", "blue", "purple"), lty=1, lwd= 3, cex=1)
dev.off()
# add horsontal gray lines




