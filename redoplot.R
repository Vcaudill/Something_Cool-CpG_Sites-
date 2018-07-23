
#data<-read.csv("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/virus/DengueVirus1.fasta_pruned.mu.trim05_DF.csv", header = T)

#data<-read.csv("DengueVirus1.fasta_pruned.mu.trim05_DF.csv")

comparing_CpG_Syn_Nonsyn_new = function(data, truename){
  
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
      AllA$graphit[i] <- 2
      AllA$mean_value[i] <- mean(AllA$MeanFreq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") )])
      AllA$sem_vals[i]<-sem(AllA$MeanFreq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") )])
    }
    if (AllA$makesCpG[i] == 1 && AllA$TypeOfSite[i] == "nonsyn") {
      AllA$graphit[i] <- 4
      AllA$mean_value[i] <- mean(AllA$MeanFreq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "nonsyn") )])
      AllA$sem_vals[i]<-sem(AllA$MeanFreq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "nonsyn") )])
    }
    if (AllA$makesCpG[i] == 0 && AllA$TypeOfSite[i] == "syn") {
      AllA$graphit[i] <- 1
      AllA$mean_value[i] <- mean(AllA$MeanFreq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") )])
      AllA$sem_vals[i]<-sem(AllA$MeanFreq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") )])
    }
    if (AllA$makesCpG[i] == 0 && AllA$TypeOfSite[i] == "nonsyn") {
      AllA$graphit[i] <- 3
      AllA$mean_value[i] <- mean(AllA$MeanFreq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "nonsyn") )])
      AllA$sem_vals[i]<-sem(AllA$MeanFreq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "nonsyn") )])
    }
  }
  
  for (i in 1:length(AllT$makesCpG)) {
    if (AllT$makesCpG[i] == 1 && AllT$TypeOfSite[i] == "syn") {
      AllT$graphit[i] <- 2
      AllT$mean_value[i] <- mean(AllT$MeanFreq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )])
      AllT$sem_vals[i]<-sem(AllT$MeanFreq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )])
    }
    if (AllT$makesCpG[i] == 1 && AllT$TypeOfSite[i] == "nonsyn") {
      AllT$graphit[i] <- 4
      AllT$mean_value[i] <- mean(AllT$MeanFreq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "nonsyn") )])
      AllT$sem_vals[i]<-sem(AllT$MeanFreq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "nonsyn") )])
    }
    if (AllT$makesCpG[i] == 0 && AllT$TypeOfSite[i] == "syn") {
      AllT$graphit[i] <- 1
      AllT$mean_value[i] <- mean(AllT$MeanFreq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )])
      AllT$sem_vals[i]<-sem(AllT$MeanFreq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )])
    }
    if (AllT$makesCpG[i] == 0 && AllT$TypeOfSite[i] == "nonsyn") {
      AllT$graphit[i] <- 3
      AllT$mean_value[i] <- mean(AllT$MeanFreq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "nonsyn") )])
      AllT$sem_vals[i]<-sem(AllT$MeanFreq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "nonsyn") )])
    }
  }
  
  for (i in 1:length(AllC$makesCpG)) {
    
    if (AllC$makesCpG[i] == 0 && AllC$TypeOfSite[i] == "syn") {
      AllC$graphit[i] <- 1
      AllC$mean_value[i] <- mean(AllC$MeanFreq[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "syn") )])
      AllC$sem_vals[i]<-sem(AllC$MeanFreq[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "syn") )])
    }
    if (AllC$makesCpG[i] == 0 && AllC$TypeOfSite[i] != "syn") {
      AllC$graphit[i] <- 3
      AllC$mean_value[i] <- mean(AllC$MeanFreq[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "nonsyn") )])
      AllC$sem_vals[i]<-sem(AllC$MeanFreq[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "nonsyn") )])
    }
  }
  
  for (i in 1:length(AllG$makesCpG)) {
    
    if (AllG$makesCpG[i] == 0 && AllG$TypeOfSite[i] == "syn") {
      AllG$graphit[i] <- 1
      AllG$mean_value[i] <- mean(AllG$MeanFreq[(which(AllG$makesCpG == 0 & AllG$TypeOfSite == "syn") )])
      AllG$sem_vals[i]<-sem(AllG$MeanFreq[(which(AllG$makesCpG == 0 & AllG$TypeOfSite == "syn") )])
    }
    if (AllG$makesCpG[i] == 0 && AllG$TypeOfSite[i] != "syn") {
      AllG$graphit[i] <- 3
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
  
  
  
  
  ####################################################################################
  #layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
  # make 2 dataset 1 conting everything except 0's (meanfreq needs to be log)
  # datea set just for meanfreq = 0's
  # graph it with breaks 
  # add a key
  # make it more like jitter
  # try to make points transparent 30% via alpha or make error darker in color
  # make bars the same width
  # try to get good png size 
  
  library(scales)
  library(plotrix)
  library(sfsmisc)
  
  #png("Den1_jitter1.png", width = 405, height = 405, units = "px")
  
  #layout(matrix(c(1,2,3,4), nrow=2, byrow = TRUE))
  truenamepng = paste(truename,".png",sep="")
  png(truenamepng, width = 6.75, height = 6.75, units = "in", res= 300)
  par(mfrow=c(2,2)) 
  palette(alpha(c("#99FF99","#9999FF","#FF9900","#FF3300"),0.3))
  #graph_color = palette(c("#99FF99","#9999FF","#FF9900","#FF3300"))
  
  plot(jitter(AllA$graphit),AllA$MeanFreq + 0.0001,log='y',col=factor(AllA$graphit),pch=16, main="A",xlab = " ", ylab = "Mutation Frequency", yaxt="n", xaxt="n")
  #yaxt="n"  aty <- axTicks(8)
  # labels <- sapply(aty,function(i)
  #   as.expression(bquote(10^ .(i)))
  # )
  # axis(2,at=aty,labels=labels)
  points(AllA$graphit, AllA$mean_val, col= factor(AllA$graphit), pch=19, cex = 3)
  arrows(AllA$graphit, AllA$LCLS, AllA$graphit, AllA$UCLS, length=0.15,lwd=5, angle=90, code=3, col= "black")
  eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
  axis(1, at= c(1:4),labels = NA)
  axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
  mtext('0', side=2, line=1.5, at=0.0001, las=1.1)
  
  # plot(data_points$Count, data_points$AnonsynNC_C,
  #      ylim=range(c(data_points$AnonsynNC_LCLS/data_points$AnonsynC_LCLS, data_points$AnonsynNC_UCLS/data_points$AnonsynC_UCLS)),
  #      pch=19, col= "green", log = 'y'
  # )
  
  plot(jitter(AllT$graphit),AllT$MeanFreq+ 0.0001,log='y',col=factor(AllT$graphit),pch=16,main="T",xlab = " ", ylab = "Mutation Frequency", yaxt="n", xaxt = "n")
  points(AllT$graphit, AllT$mean_val, col= factor(AllT$graphit), pch=19, cex = 3)
  arrows(AllT$graphit, AllT$LCLS, AllT$graphit, AllT$UCLS, length=0.15, lwd = 5, angle=90, code=3, col= "black")
  eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
  axis(1, at= c(1:4),labels = NA)
  axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
  mtext('0', side=2, line=1.5, at=0.0001, las=1.1)
  
  
  palette(alpha(c("#99FF99","#FF9900"),0.3)) 
  plot(jitter(AllC$graphit, 0.6),AllC$MeanFreq+ 0.0001,log='y',col=factor(AllC$graphit),pch=16,main="C",xlab = "Mutation Type", ylab = "Mutation Frequency",  xlim=c(0.5,4.5),yaxt="n", xaxt = "n")
  
  points(AllC$graphit, AllC$mean_val, col= factor(AllC$graphit), pch=19, cex = 3)
  arrows(AllC$graphit, AllC$LCLS, AllC$graphit, AllC$UCLS, length=0.15,lwd = 5, angle=90, code=3, col= "black")
  eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
  axis(1, at= c(1:4),labels = NA)
  axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
  mtext('0', side=2, line=1.5, at=0.0001, las=1.1)
  
  
  plot(jitter(AllG$graphit, 0.6),AllG$MeanFreq+ 0.0001,log='y',col=factor(AllG$graphit),pch=16,main="G",xlab = "Mutation Type", ylab = "Mutation Frequency", xlim=c(0.5,4.5),yaxt="n", xaxt = "n")
  points(AllG$graphit, AllG$mean_val, col= factor(AllG$graphit), pch=19, cex = 3)
  arrows(AllG$graphit, AllG$LCLS, AllG$graphit, AllG$UCLS, length=0.15, lwd=5,angle=90, code=3, col= "black")
  eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
  axis(1, at= c(1:4),labels = NA)
  axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
  mtext('0', side=2, line=1.5, at=0.0001, las=1.1)
  
  dev.off()
  
  
  
  
  #gap.plot(G0$MeanFreq,gap=c(0.0001,2),add = TRUE)
  #plot(Gfreq$graphit,Gfreq$MeanFreq,log='y',col=factor(Gfreq$graphit),pch=16,main="G",xlab = "Mutation Type", ylab = "Mutation Frequency")
  #axis.break(2,0.01,style="gap")
  #axis.break(2,0.001*(1+0.02),breakcol="black",style="slash")
  
  #palette(c("#99FF99","#9999FF","#FF9900","#FF3300"))
  #gap.plot(A0$MeanFreq,gap=c(0.0001,2),col=factor(A0$graphit),main="A",  add = TRUE)
  #plot(Afreq$graphit,Afreq$MeanFreq,col=factor(Afreq$graphit))
  #axis.break(2,0.001,breakcol="snow",style="gap")
  #axis.break(2,0.0001*(1+0.02),breakcol="black",style="slash")
  
  #palette(c("#99FF99","#9999FF","#FF9900","#FF3300"))
  #plot(AllA$graphit, AllA$MeanFreq, log='y',col = factor(AllA$graphit),xlab="Mutation Type", ylab="MeanFreq",main="A", pch=16)
  #axis.break(2,0.001,style="slash") 
  #gap.plot(A0$graphit,A0$MeanFreq, gap=c(8,16))
  #par(new=TRUE)
  #plot(A0$graphit, A0$MeanFreq,pch=16)
  
  
  #gap.plot(A0$MeanFreq,gap=c(0.0001,2),col=factor(AllA$graphit),main="barplot with gap")
  #axis.break(2,0.0001,breakcol="snow",style="gap")
  #axis.break(2,0.0001*(1+0.02),breakcol="black",style="slash")
  #axis.break(4,0.0001*(1+0.02),breakcol="black",style="slash")
  #axis(2,at=0.0001)
  
  #axis(1, at=1:4, label)
  #xlim=c(0.5,4.5)
  
  
  # color by grpahit 
  plot(AllT$graphit, AllT$MeanFreq, log='y',col = factor(AllT$graphit),xlab="Mutation Type", ylab="MeanFreq",main="T")
  plot(AllG$graphit, AllG$MeanFreq, log='y',col = factor(AllG$graphit),xlab="Mutation Type", ylab="MeanFreq",main="G", xlim = c(-1,10))
  
  # color by grpahit 
  plot(AllC$graphit, AllC$MeanFreq, log='y',col = factor(AllC$graphit),xlab="Mutation Type", ylab="MeanFreq",main="C")
  
  
  
  
  
  col=c("red","blue","black","orange")[AllT$graphit]
  pch=c(16:18)[AllT$graphit]
  
  plot(AllA$makesCpG==1&AllA$TypeOfSite=="syn",AllA$MeanFreq+.000001,log='y', col = "red",xlab="Mutation Type", ylab="MeanFreq",main="A", xlim = c(-1,10))
  par(new=TRUE)
  plot(AllA$makesCpG==0&AllA$TypeOfSite=="syn",AllA$MeanFreq+.000001,log='y', col = "blue",xlab="Mutation Type", ylab="MeanFreq",main="A", xlim = c(-1,10))
  par(new=TRUE)
  plot(AllA$makesCpG==1&AllA$TypeOfSite!="syn",AllA$MeanFreq+.000001,log='y', col = "yellow",xlab="Mutation Type", ylab="MeanFreq",main="A")
  par(new=TRUE)
  plot(AllA$makesCpG==0&AllA$TypeOfSite!="syn",AllA$MeanFreq+.000001,log='y', col = "green",xlab="Mutation Type", ylab="MeanFreq",main="A")
  
  par(new=TRUE)
  plot(AllG$makesCpG==1&AllG$TypeOfSite=="syn",AllG$MeanFreq+.000001,log='y', col = "red",xlab="Mutation Type", ylab="MeanFreq",main="G")
  par(new=TRUE)
  plot(AllG$makesCpG==0&AllG$TypeOfSite=="syn",AllG$MeanFreq+.000001,log='y', col = "blue",xlab="Mutation Type", ylab="MeanFreq",main="G")
  par(new=TRUE)
  plot(AllG$makesCpG==1&AllG$TypeOfSite!="syn",AllG$MeanFreq+.000001,log='y', col = "yellow",xlab="Mutation Type", ylab="MeanFreq",main="G")
  par(new=TRUE)
  plot(AllG$makesCpG==0&AllG$TypeOfSite!="syn",AllG$MeanFreq+.000001,log='y', col = "green",xlab="Mutation Type", ylab="MeanFreq",main="G")
  
  
  ######################################################################################
  
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


#comparing_CpG_Syn_Nonsyn(DF)
#make sure the packages are loaded before begining 

###################
#still testing
# 
# ggplot(aes(factor(graphit), MeanFreq, color=graphit), data = AllAT)+
#   #log scale to make the data eaisier to see
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   scale_x_discrete(labels=c("1" = "CpG ", "5" = "CpG ", "3" = "nonCpG ", "4"= "nonCpG "))+
#   geom_boxplot(data= AllA,aes(x = factor(graphit))) +
#   #facet_wrap splits graph between a and t
#   facet_wrap(~ wtnt)+
#   geom_boxplot(data= AllT,aes(x = factor(graphit))) +
#   #give points new colors and lables the colors
#   scale_color_manual(labels = c("CpG (syn)","Cpg (nonsyn)","nonCpG (syn)", "nonCpg (nonsyn)"), values = c("firebrick", "darkolivegreen","goldenrod3", "royalblue3")) +
#   #labels X and Y axis
#   labs(x="Mutation Type", y="Mutation Frquency",col=" ")+
#   annotation_logticks(sides="l") 
# 
# 
# 
# 
# plot(MeanFreq~TypeOfSite, data=AC, log='y',yaxt='n',xaxt='n',xlab="",ylab="",col=c("#601E0033","#FF000033","#0000FF33","#FF5B0033","#00FDFF33"), border=c("#601E0080","#FF000080","#0000FF80","#FF5B00","#00FDFF"))
# par(new=TRUE)
# plot(MeanFreq~TypeOfSite, data=ANC, log='y',yaxt='n',xaxt='n',xlab="",ylab="",col=c("#0000FF33","#FF5B0033","#00FDFF33"), border=c("#0000FF80","#FF5B00","#00FDFF"))
# 
# eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1)
# axis.break(2,0.00012,style="slash") 
# mtext('0', side=2, line=1.5, at=0.0001, las=1,cex=1.1)
# #mtext("Mutation type", side=1, line=8, cex=1.5)
# mtext("Mutation frequency", side=2, line=3, cex=1.1)
# 
# points(MF1~Type, data=stats_s,log='y',pch=16,cex=1)
