
data<-read.csv("DengueVirus1.fasta_pruned.mu.trim05_DF.csv", header = T)



comparing_CpG_Syn_Nonsyn = function(data){
  library(ggplot2)
  library(dplyr)
  library(plotrix)
  library(scales)
  library(grid)
  library(gridExtra)
  # my solution for the error I was getting 
  require(extrafont)
  # # need only do this once!
  font_import(pattern="[A/a]rial", prompt=FALSE)
  require(ggplot2)

#subset into two groups yes makes cpg and no cpg
cpg.y<-subset(data, makesCpG==1)
cpg.n<-subset(data, makesCpG==0)
#subset further into letters nuclotideCpgforming or nucotideNonGpg
AC<-subset(cpg.y, wtnt=='a')
ANC<-subset(cpg.n, wtnt=='a')
GC<-subset(cpg.y, wtnt=='g')
GNC<-subset(cpg.n, wtnt=='g')
TC<-subset(cpg.y, wtnt=='t')
TNC<-subset(cpg.n, wtnt=='t')
CC<-subset(cpg.y, wtnt=='c')
CNC<-subset(cpg.n, wtnt=='c')

#Function to help create errorbars
sem<-function(x){
  return(sd(x,na.rm = FALSE)/sqrt(length(x)))
}

# addinf a small value to be able to graph on a log scale
# AC$MeanFreq = AC$MeanFreq +.00001
# ANC$MeanFreq = ANC$MeanFreq +.00001
# TC$MeanFreq = TC$MeanFreq +.00001
# TNC$MeanFreq = TNC$MeanFreq +.00001
# CNC$MeanFreq = CNC$MeanFreq +.00001
# GNC$MeanFreq = GNC$MeanFreq +.00001

#making the data frames with all information about a, t, c, g 
AllA = rbind(AC, ANC)
AllT = rbind(TC, TNC)
AllC = CNC
AllG = GNC

# added new colmuns for mean an error bars
AllA$mean_value <- .01
AllA$sem_vals<- 0
AllT$mean_value <- .01
AllT$sem_vals<- 0
AllC$mean_value <- .01
AllC$sem_vals<- 0
AllG$mean_value <- .01
AllG$sem_vals<- 0

  
# for loops to caculate mean ans errorbars and 1, 2, 3, 4 for position
for (i in 1:length(AllA$makesCpG)) {
  if (AllA$makesCpG[i] == 1 && AllA$TypeOfSite[i] == "syn") {
    AllA$graphit[i] <- "1"
    AllA$mean_value[i] <- mean(AllA$MeanFreq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") )])
    AllA$sem_vals[i]<-sem(AllA$MeanFreq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") )])
    }
  if (AllA$makesCpG[i] == 1 && AllA$TypeOfSite[i] == "nonsyn") {
    AllA$graphit[i] <- "3"
    AllA$mean_value[i] <- mean(AllA$MeanFreq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "nonsyn") )])
    AllA$sem_vals[i]<-sem(AllA$MeanFreq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "nonsyn") )])
    }
  if (AllA$makesCpG[i] == 0 && AllA$TypeOfSite[i] == "syn") {
    AllA$graphit[i] <- "2"
    AllA$mean_value[i] <- mean(AllA$MeanFreq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") )])
    AllA$sem_vals[i]<-sem(AllA$MeanFreq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") )])
  }
  if (AllA$makesCpG[i] == 0 && AllA$TypeOfSite[i] == "nonsyn") {
    AllA$graphit[i] <- "4"
    AllA$mean_value[i] <- mean(AllA$MeanFreq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "nonsyn") )])
    AllA$sem_vals[i]<-sem(AllA$MeanFreq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "nonsyn") )])
  }
}

for (i in 1:length(AllT$makesCpG)) {
  if (AllT$makesCpG[i] == 1 && AllT$TypeOfSite[i] == "syn") {
    AllT$graphit[i] <- "1"
    AllT$mean_value[i] <- mean(AllT$MeanFreq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )])
    AllT$sem_vals[i]<-sem(AllT$MeanFreq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )])
  }
  if (AllT$makesCpG[i] == 1 && AllT$TypeOfSite[i] == "nonsyn") {
    AllT$graphit[i] <- "3"
    AllT$mean_value[i] <- mean(AllT$MeanFreq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "nonsyn") )])
    AllT$sem_vals[i]<-sem(AllT$MeanFreq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "nonsyn") )])
  }
  if (AllT$makesCpG[i] == 0 && AllT$TypeOfSite[i] == "syn") {
    AllT$graphit[i] <- "2"
    AllT$mean_value[i] <- mean(AllT$MeanFreq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )])
    AllT$sem_vals[i]<-sem(AllT$MeanFreq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )])
  }
  if (AllT$makesCpG[i] == 0 && AllT$TypeOfSite[i] == "nonsyn") {
    AllT$graphit[i] <- "4"
    AllT$mean_value[i] <- mean(AllT$MeanFreq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "nonsyn") )])
    AllT$sem_vals[i]<-sem(AllT$MeanFreq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "nonsyn") )])
  }
}

for (i in 1:length(AllC$makesCpG)) {

  if (AllC$makesCpG[i] == 0 && AllC$TypeOfSite[i] == "syn") {
    AllC$graphit[i] <- "2"
    AllC$mean_value[i] <- mean(AllC$MeanFreq[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "syn") )])
    AllC$sem_vals[i]<-sem(AllC$MeanFreq[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "syn") )])
  }
  if (AllC$makesCpG[i] == 0 && AllC$TypeOfSite[i] != "syn") {
    AllC$graphit[i] <- "4"
    AllC$mean_value[i] <- mean(AllC$MeanFreq[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "nonsyn") )])
    AllC$sem_vals[i]<-sem(AllC$MeanFreq[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "nonsyn") )])
  }
}

for (i in 1:length(AllG$makesCpG)) {
  
  if (AllG$makesCpG[i] == 0 && AllG$TypeOfSite[i] == "syn") {
    AllG$graphit[i] <- "2"
    AllG$mean_value[i] <- mean(AllG$MeanFreq[(which(AllG$makesCpG == 0 & AllG$TypeOfSite == "syn") )])
    AllG$sem_vals[i]<-sem(AllG$MeanFreq[(which(AllG$makesCpG == 0 & AllG$TypeOfSite == "syn") )])
  }
  if (AllG$makesCpG[i] == 0 && AllG$TypeOfSite[i] != "syn") {
    AllG$graphit[i] <- "4"
    AllG$mean_value[i] <- mean(AllG$MeanFreq[(which(AllG$makesCpG == 0 & AllG$TypeOfSite == "nonsyn") )])
    AllG$sem_vals[i]<-sem(AllG$MeanFreq[(which(AllG$makesCpG == 0 & AllG$TypeOfSite == "nonsyn") )])
  }
}

# There are the upper and lower limits of the error bar
AllA$LCLS = AllA$mean_value - AllA$sem_vals
AllA$UCLS = AllA$mean_value + AllA$sem_vals

AllT$LCLS = AllT$mean_value - AllT$sem_vals
AllT$UCLS = AllT$mean_value + AllT$sem_vals

AllC$LCLS = AllC$mean_value - AllC$sem_vals
AllC$UCLS = AllC$mean_value + AllC$sem_vals

AllG$LCLS = AllG$mean_value - AllG$sem_vals
AllG$UCLS = AllG$mean_value + AllG$sem_vals

AllAT = rbind(AllA, AllT)
AllATCG = rbind(AllA, AllT, AllC, AllG)

# this loop adds a small value to the 0's so they may show up on the graph
for (i in 1:length(AllATCG$MeanFreq)){
  if (AllATCG$MeanFreq[i]==0){
    AllATCG$MeanFreq[i]= 0.000001
  }
}

# the plot
ggplot(aes(factor(graphit), MeanFreq, color=graphit), data = AllATCG)+
  #log scale to make the data eaisier to see
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #expand_limits(y = c(0.00001, 0.1)) +
  # the title
  labs(title = "Transition Mutations")+
  # jitter makes the points spred out width is how far apart the points are and alpha deals with opacicty
  geom_jitter(data= AllA,aes(x = factor(graphit)),position = position_jitter(width = .2), alpha = 0.3) +
  # point is the mean
  geom_point(data= AllA,mapping = aes (x = graphit, y =mean_value, colour =graphit ),size = 5.0, show.legend = FALSE) +
  # error bar is the uper and lower limits
  geom_errorbar(data = AllA,aes(x = graphit, ymin= LCLS, ymax= UCLS, color = graphit),width=.5) +
  #facet_wrap splits graph between a, t, c, g
  facet_wrap(~ wtnt)+
  geom_jitter(data= AllT,aes(x = factor(graphit)), position = position_jitter(width = .2), alpha = 0.3) +
  geom_point(data= AllT,mapping = aes (x = graphit, y =mean_value, colour =graphit ),size = 5.0, show.legend = FALSE) +
  geom_errorbar(data = AllT,aes(x = graphit, ymin= LCLS, ymax= UCLS, color = graphit),width=.5) +
  
  geom_jitter(data= AllC,aes(x = factor(graphit)),position = position_jitter(width = .2), alpha = 0.3) +
  geom_point(data= AllC,mapping = aes (x = graphit, y =mean_value, colour =graphit ),size = 5.0, show.legend = FALSE) +
  geom_errorbar(data = AllC,aes(x = graphit, ymin= LCLS, ymax= UCLS, color = graphit),width=.5) +
  
  geom_jitter(data= AllG,aes(x = factor(graphit)),position = position_jitter(width = .2), alpha = 0.3) +
  geom_point(data= AllG,mapping = aes (x = graphit, y =mean_value, colour =graphit ),size = 5.0, show.legend = FALSE) +
  geom_errorbar(data = AllG,aes(x = graphit, ymin= LCLS, ymax= UCLS, color = graphit),width=.5) +
 
   #give points new colors and lables the colors
  scale_color_manual(labels = c("CpG (syn)","nonCpG (syn)","CpG (nonsyn)", "nonCpG (nonsyn)"), values = c("firebrick", "royalblue3","goldenrod3", "darkolivegreen")) +
  #labels X and Y axis
  labs(x="Mutation Type", y="Mutation Frquency",col=" ")+
  #theme(text=element_text(family="Garamond", size=14))+
  annotation_logticks(sides="l") 


}


comparing_CpG_Syn_Nonsyn(data)
#make sure the packages are loaded before begining 

###################
#still testing

ggplot(aes(factor(graphit), MeanFreq, color=graphit), data = AllAT)+
  #log scale to make the data eaisier to see
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_discrete(labels=c("1" = "CpG ", "5" = "CpG ", "3" = "nonCpG ", "4"= "nonCpG "))+
  geom_boxplot(data= AllA,aes(x = factor(graphit))) +
  #facet_wrap splits graph between a and t
  facet_wrap(~ wtnt)+
  geom_boxplot(data= AllT,aes(x = factor(graphit))) +
  #give points new colors and lables the colors
  scale_color_manual(labels = c("CpG (syn)","Cpg (nonsyn)","nonCpG (syn)", "nonCpg (nonsyn)"), values = c("firebrick", "darkolivegreen","goldenrod3", "royalblue3")) +
  #labels X and Y axis
  labs(x="Mutation Type", y="Mutation Frquency",col=" ")+
  annotation_logticks(sides="l") 




plot(MeanFreq~TypeOfSite, data=AC, log='y',yaxt='n',xaxt='n',xlab="",ylab="",col=c("#601E0033","#FF000033","#0000FF33","#FF5B0033","#00FDFF33"), border=c("#601E0080","#FF000080","#0000FF80","#FF5B00","#00FDFF"))
par(new=TRUE)
plot(MeanFreq~TypeOfSite, data=ANC, log='y',yaxt='n',xaxt='n',xlab="",ylab="",col=c("#0000FF33","#FF5B0033","#00FDFF33"), border=c("#0000FF80","#FF5B00","#00FDFF"))

eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1)
axis.break(2,0.00012,style="slash") 
mtext('0', side=2, line=1.5, at=0.0001, las=1,cex=1.1)
#mtext("Mutation type", side=1, line=8, cex=1.5)
mtext("Mutation frequency", side=2, line=3, cex=1.1)

points(MF1~Type, data=stats_s,log='y',pch=16,cex=1)
