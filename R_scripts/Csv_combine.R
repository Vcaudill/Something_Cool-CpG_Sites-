library(seqinr)
library(ggplot2)
library(scales)
library(grid)
library(plotrix)
hyphy_virus<-read.csv("data/MEME_CpG_list1.csv")
i=1
for (i in 1: nrow(hyphy_virus)){
#these viruses dont have datamonkey files yet, except humanpap not using
if(hyphy_virus$name[i]== "DengueVirus1.fasta_pruned.mu.trim05"){
  next
}
if(hyphy_virus$name[i]== "DengueVirus2.fasta_pruned.mu.trim05"){
  next
}
if(hyphy_virus$name[i]== "DengueVirus3.fasta_pruned.mu.trim05"){
  next
}
if(hyphy_virus$name[i]== "InfluenzaAvirus_HA_H3N2.fasta.mu.trim05"){
    next
}
if(hyphy_virus$name[i]== "Humanpapillomavirus16.fasta_pruned.mu.trim05"){
    next
}
if(hyphy_virus$name[i]== "InfluenzaBvirus_HA.fasta.mu.trim05"){
    next
  }
  
  
  

name = as.character(hyphy_virus$name[i])
splitname<-unlist(strsplit(as.character(hyphy_virus$name[i]),".fasta"))
truename<-splitname[1]  
path_file<-paste("Hyphy/Consensus_Hyphy/",truename,".csv", sep="")
al1<-read.csv(path_file)
folder<-paste(hyphy_virus$Hyphy_folders[i],sep='')
path<-paste("Hyphy/",folder,"/datamonkey-table.csv",sep='')
datamonkey<-read.csv(path)
file<-al1$wtnt_consensus
nuc<-splitseq(file, frame = 0, word = 3)
makesCpG<-al1$makesCpG
CpG<-splitseq(makesCpG, frame = 0, word = 3)

datamonkey$consensus_codon<-nuc
datamonkey$consensus_AA<-translate(as.character(file))
datamonkey$potential_CpG<-"no"
datamonkey$makesCpG<-0
datamonkey$graphit<-0
datamonkey$graphit2<-0
datamonkey$graphit3<-0
for(j in 1:nrow(datamonkey)){
  letter<-datamonkey$consensus_AA[j]
  if(letter=="H"||letter=="E"||letter=="Y"||letter=="C"||letter=="F"||letter=="N"||letter=="K"||letter=="Q"||letter=="D"){
    datamonkey$potential_CpG[j]="yes"
  }
  if(sum(as.numeric(CpG[j])) > 0){
    datamonkey$makesCpG[j] = 1
  }
  if(datamonkey$makesCpG[j]==1 && datamonkey$potential_CpG[j] == "yes"){
    datamonkey$graphit[j]="1"
  }
  if(datamonkey$makesCpG[j]==1 && datamonkey$potential_CpG[j] == "yes"){
    datamonkey$graphit3[j]="1"
  }
  if(datamonkey$makesCpG[j]==0 && datamonkey$potential_CpG[j] == "yes"){
    datamonkey$graphit3[j]="2"
  }
  if(datamonkey$potential_CpG[j] == "yes"){
    datamonkey$graphit2[j]="1"
  }
  if(datamonkey$potential_CpG[j] == "no"){
    datamonkey$graphit2[j]="2"
  }
  if(datamonkey$makesCpG[j]==0){
    datamonkey$graphit[j]="2"
  }
  if(datamonkey$makesCpG[j]==1 && datamonkey$potential_CpG[j] == "no"){
    datamonkey$graphit[j]="3"
  }
}

write.csv(datamonkey,path)

}

#########
datamonkey<-read.csv("path.....")
#meme&fubar http://datamonkey.org/fubar/5cb3c3463994747a2e471e19 BK fubar
cpg.y<-subset(datamonkey, makesCpG==1 & potential_CpG == "yes")
#makes CpG in genetic code and is 2 fold 
cpg.n<-subset(datamonkey, makesCpG==0)
# does not make CpG in genetic code 
cpg.m<-subset(datamonkey, makesCpG==1 & potential_CpG == "no")
#does make CpG in genetic code for transition mutation or transversion, but may be 4 fold
ggplot(aes(factor(graphit),X.alpha.,color= graphit),data=datamonkey)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title = "BK Alpha Score Meme or Fubar")+
  geom_jitter(data=cpg.y,aes(x=factor(graphit)),position = position_jitter(width = .2), alpha = 0.8)+
  geom_jitter(data=cpg.n,aes(x=factor(graphit)),position = position_jitter(width = .2), alpha = 0.8)+
  geom_jitter(data=cpg.m,aes(x=factor(graphit)),position = position_jitter(width = .2), alpha = 0.8)+
  scale_color_manual(labels = c("makesCpG","nonCpG","maybeCpG"), values = c("darkgreen", "darkred","gold")) +
  #labels X and Y axis
  labs(x="CpG Type", y="Alpha Value",col=" ")+
  annotation_logticks(sides="l") 


cpg.py<-subset(datamonkey, makesCpG==1 & potential_CpG == "yes")
#makes CpG in genetic code and is 2 fold 
cpg.pn<-subset(datamonkey, potential_CpG == "no")
# does not make CpG in genetic code 
ggplot(aes(factor(graphit2),X.alpha., color=graphit2),data=datamonkey)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title = "BK Alpha Score Meme or Fubar")+
  geom_jitter(data=cpg.py,aes(x=factor(graphit2)),position = position_jitter(width = .2), alpha = 0.8)+
  geom_jitter(data=cpg.pn,aes(x=factor(graphit2)),position = position_jitter(width = .2), alpha = 0.8)+
  scale_color_manual(labels = c("Yes CpG","No CpG"), values = c("darkgreen", "darkred")) +
  #labels X and Y axis
  labs(x="CpG Type", y="Alpha Value",col=" ")+
  annotation_logticks(sides="l") 


cpg.yc<-subset(datamonkey, makesCpG==1 & potential_CpG == "yes")
#makes CpG in genetic code and is 2 fold 
cpg.ync<-subset(datamonkey, makesCpG==0 & potential_CpG == "yes")
# does not make CpG in genetic code 
ggplot(aes(factor(graphit3),X.alpha., color=graphit3),data=datamonkey)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title = "BK Alpha Score Fubar (from possible trasition mutations only)")+
  geom_jitter(data=cpg.yc,aes(x=factor(graphit3)),position = position_jitter(width = .2), alpha = 0.8)+
  geom_jitter(data=cpg.ync,aes(x=factor(graphit3)),position = position_jitter(width = .2), alpha = 0.8)+
  scale_color_manual(labels = c("Yes CpG","No CpG"), values = c("darkgreen", "darkred")) +
  #labels X and Y axis
  labs(x="CpG Type", y="Alpha Value",col=" ")+
  annotation_logticks(sides="l") 
print(wilcox.test(cpg.yc$X.alpha., cpg.ync$X.alpha., alternative='less'))


array1 = datamonkey$X.alpha.[datamonkey$potential_CpG == "yes" & datamonkey$makesCpG == 0]
array2 = datamonkey$X.alpha.[datamonkey$potential_CpG == "yes" & datamonkey$makesCpG == 1]

print(wilcox.test(array1, array2, alternative='greater'))

  
