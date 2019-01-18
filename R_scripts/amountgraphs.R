# graph for data vs significant 

library(scales)
library(plotrix)
library(sfsmisc)
library(gdata)
Virus_info<- read.csv("data/CpG_List.csv")

data_points_df<- read.csv("output/alldatapoints.csv")
#data_points_df=data_points_df[-(27), ] 
data_points_df <- drop.levels(data_points_df[-(27), ]) 
tally<- read.csv("data/recount2.csv")
data_points_df$SynCpg=0
data_points_df$NonSynCpg=0
data_points_df$SynNonSyn=0
for (i in 1:nrow(data_points_df)) {
 treat=data_points_df$Virus[i]
  
  print(treat)
 # treat=factor(data_points_df$Virus[i])
  for (j in 1:nrow(tally)){
    # splitnameAll<-unlist(strsplit(as.character(Virus_info$name[j]),".fasta"))
    # splitname<-unlist(strsplit(as.character(splitnameAll[1]),"_"))
   truename<-tally$Virus[j]
    print(truename)
   # truename=factor(tally$Virus[i])
    if(treat==truename){
      data_points_df$SynCpg[i]<-tally$A..G.Syn..CpG.v.NonCpG.[j]+tally$T..C.Syn..CpG.v.NonCpG.[j]
      data_points_df$NonSynCpg[i]<-tally$A..G.NonSyn..CpG.v.NonCpG[j]+tally$T..C.NonSyn..CpG.v.NonCpG[j]
      data_points_df$SynNonSyn[i]<-tally$A..G.Syn.v.NonSyn[j]+tally$T..C.Syn.v.NonSyn[j]
      
      
    }
  } 
}

#Syn CpG vs Non-CpG
  

  png("output/amount/Syn.png", width = 6.75, height = 6.75, units = "in", res= 300)
  #par(mfrow=c(2,2))#, bg = "darkseagreen1"
  palette(alpha(c("red","green","deepskyblue1")))
  plot(data_points_df$seq,data_points_df$nucl,log='xy',col=factor(data_points_df$SynCpg),pch=19,cex=1.5, main="Syn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")
  legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
         col= c("deepskyblue1","green","red"), horiz= FALSE, 
         #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
         cex=1, pch = c(16,16,16))
  
  png("output/amount/NonSyn.png", width = 6.75, height = 6.75, units = "in", res= 300)
  palette(alpha(c("red","green","deepskyblue1")))
  plot(data_points_df$seq,data_points_df$nucl,log='xy',col=factor(data_points_df$NonSynCpg),pch=19,cex=1.5, main="NonSyn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")
  legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
         col= c("deepskyblue1","green","red"), horiz= FALSE, 
         #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
         cex=1, pch = c(16,16,16))
  png("output/amount/SynNonSyn.png", width = 6.75, height = 6.75, units = "in", res= 300)
  palette(alpha(c("red","green","deepskyblue1")))
  plot(data_points_df$seq,data_points_df$nucl,log='xy',col=factor(data_points_df$SynNonSyn),pch=19,cex=1.5, main="Syn vs NonSyn",xlab = "# of Sequences", ylab = "# of Nucletides")
  legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
         col= c("deepskyblue1","green","red"), horiz= FALSE, 
         #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
         cex=1, pch = c(16,16,16))
dev.off()

png("output/amount/alllogy.png", width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(2,2))#, bg = "darkseagreen1"
palette(alpha(c("red","green","deepskyblue1")))
plot(data_points_df$seq,data_points_df$nucl,log='y',col=factor(data_points_df$SynCpg),pch=19,cex=1.5, main="Syn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(16,16,16))


palette(alpha(c("red","green","deepskyblue1")))
plot(data_points_df$seq,data_points_df$nucl,log='y',col=factor(data_points_df$NonSynCpg),pch=19,cex=1.5, main="NonSyn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(16,16,16))

palette(alpha(c("red","green","deepskyblue1")))
plot(data_points_df$seq,data_points_df$nucl,log='y',col=factor(data_points_df$SynNonSyn),pch=19,cex=1.5, main="Syn vs NonSyn",xlab = "# of Sequences", ylab = "# of Nucletides")
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(16,16,16))
dev.off()

png("output/amount/alllogxy.png", width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(2,2))#, bg = "darkseagreen1"
palette(alpha(c("red","green","deepskyblue1")))
plot(data_points_df$seq,data_points_df$nucl,log='xy',col=factor(data_points_df$SynCpg),pch=19,cex=1.5, main="Syn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")



palette(alpha(c("red","green","deepskyblue1")))
plot(data_points_df$seq,data_points_df$nucl,log='xy',col=factor(data_points_df$NonSynCpg),pch=19,cex=1.5, main="NonSyn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")

palette(alpha(c("red","green","deepskyblue1")))
plot(data_points_df$seq,data_points_df$nucl,log='xy',col=factor(data_points_df$SynNonSyn),pch=19,cex=1.5, main="Syn vs NonSyn",xlab = "# of Sequences", ylab = "# of Nucletides")
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(16,16,16))
dev.off()

png("output/amount/allnolog.png", width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(2,2))#, bg = "darkseagreen1"
palette(alpha(c("red","green","deepskyblue1")))
plot(data_points_df$seq,data_points_df$nucl,col=factor(data_points_df$SynCpg),pch=19,cex=1.5, main="Syn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(16,16,16))


palette(alpha(c("red","green","deepskyblue1")))
plot(data_points_df$seq,data_points_df$nucl,col=factor(data_points_df$NonSynCpg),pch=19,cex=1.5, main="NonSyn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(16,16,16))

palette(alpha(c("red","green","deepskyblue1")))
plot(data_points_df$seq,data_points_df$nucl,col=factor(data_points_df$SynNonSyn),pch=19,cex=1.5, main="Syn vs NonSyn",xlab = "# of Sequences", ylab = "# of Nucletides")
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(16,16,16))
dev.off()

