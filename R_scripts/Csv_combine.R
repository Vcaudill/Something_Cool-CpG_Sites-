library(seqinr)
library(ggplot2)
library(scales)
library(grid)
library(plotrix)
#reads in MEME.CSV with all the viruses
hyphy_virus<-read.csv("data/MEME_CpG_list1.csv")
i=1
special_folder<-"meme"
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
  if(hyphy_virus$name[i]== "HepatitisB_precore.fasta_pruned.mu.trim05"){
    next
  }#Hep_precore no presence of nonCpG mutation 
  if(hyphy_virus$name[i]== "HepatitisB_pre_S.fasta_pruned.mu.trim05"){
    next
  }#Hep_precore no presence of nonCpG mutation for a last nuc

  name = as.character(hyphy_virus$name[i])
  splitname<-unlist(strsplit(as.character(hyphy_virus$name[i]),".fasta"))
  #making truename for shorter virus name
  truename<-splitname[1]  
  path_file<-paste("Hyphy/Consensus_Hyphy/",truename,".csv", sep="")
  al1<-read.csv(path_file)
  folder<-paste(hyphy_virus$Hyphy_folders[i],sep='')
  #accessing the meme folders 
  path<-paste("Hyphy/",folder,"/meme/datamonkey-table.csv",sep='')
  #path<-paste("Hyphy/",folder,"/fubar/datamonkey-table.csv",sep='')
  #reading in the meme csv for virus
  datamonkey<-read.csv(path)
  #making a new column for consensus codons
  file<-al1$wtnt_consensus
  nuc<-splitseq(file, frame = 0, word = 3)
  makesCpG<-al1$makesCpG
  CpG<-splitseq(makesCpG, frame = 0, word = 3)
  l_nuc<-splitseq(file, frame = 0, word = 1)
  last_nuc<-l_nuc[seq(3,length(l_nuc),3)]
  datamonkey$consensus_codon<-nuc
  datamonkey$last_nuc<-last_nuc
  datamonkey$consensus_AA<-translate(as.character(file))
  #new columns
  datamonkey$potential_CpG<-"no"
  datamonkey$makesCpG<-0
  datamonkey$graphit4<-0
  #need to calculate standard error
  sem<-function(x){
    return(sd(x,na.rm = FALSE)/sqrt(length(x)))
  }
  for(j in 1:nrow(datamonkey)){
    letter<-datamonkey$consensus_AA[j]
    #write another if if a or flag it if last nuc is a or t then flag as potential sites
    if(letter=="H"||letter=="E"||letter=="Y"||letter=="C"||letter=="F"||letter=="N"||letter=="K"||letter=="Q"||letter=="D"){
      if(datamonkey$last_nuc[j]== "a"||datamonkey$last_nuc[j]== "t"){
        datamonkey$potential_CpG[j]="yes"
      }
    }
    if(sum(as.numeric(CpG[j])) > 0){
      datamonkey$makesCpG[j] = 1
    }
    if(datamonkey$makesCpG[j]==1 && datamonkey$potential_CpG[j] == "yes"&& datamonkey$last_nuc[j]== 'a'){
      datamonkey$graphit4[j]="2"
    }
    if(datamonkey$makesCpG[j]==0 && datamonkey$potential_CpG[j] == "yes"&& datamonkey$last_nuc[j]== 'a'){
      datamonkey$graphit4[j]="1"
    } 
    if(datamonkey$makesCpG[j]==1 && datamonkey$potential_CpG[j] == "yes"&& datamonkey$last_nuc[j]== 't'){
      datamonkey$graphit4[j]="4"
    }
    if(datamonkey$makesCpG[j]==0 && datamonkey$potential_CpG[j] == "yes"&& datamonkey$last_nuc[j]== 't'){
      datamonkey$graphit4[j]="3"
    }
  
  }
  
  write.csv(datamonkey,path)
  


######### graph

nice_name <- as.character(hyphy_virus$nice_name[i])
#meme&fubar http://datamonkey.org/fubar/5cb3c3463994747a2e471e19 BK fubar


cpg.yca<-subset(datamonkey, makesCpG==1 & potential_CpG == "yes" & last_nuc == "a" )
#makes CpG in genetic code and is 2 fold 
cpg.ynca<-subset(datamonkey, makesCpG==0 & potential_CpG == "yes" & last_nuc == "a")
# does not make CpG in genetic code 
cpg.yct<-subset(datamonkey, makesCpG==1 & potential_CpG == "yes" & last_nuc == "t" )
#makes CpG in genetic code and is 2 fold 
cpg.ynct<-subset(datamonkey, makesCpG==0 & potential_CpG == "yes" & last_nuc == "t")

cpg.yca$mean_value <- .01
cpg.yca$sem_vals<- 0
cpg.ynca$mean_value <- .01
cpg.ynca$sem_vals<- 0
cpg.yct$mean_value <- .01
cpg.yct$sem_vals<- 0
cpg.ynct$mean_value <- .01
cpg.ynct$sem_vals<- 0

cpg.yca$mean_value <- mean(cpg.yca$X.alpha)
cpg.yca$sem_vals<-sem(cpg.yca$X.alpha)
cpg.ynca$mean_value <- mean(cpg.ynca$X.alpha)
cpg.ynca$sem_vals<-sem(cpg.ynca$X.alpha)
cpg.yct$mean_value <- mean(cpg.yct$X.alpha)
cpg.yct$sem_vals<-sem(cpg.yct$X.alpha)
cpg.ynct$mean_value <- mean(cpg.ynct$X.alpha)
cpg.ynct$sem_vals<-sem(cpg.ynct$X.alpha)

# There are the upper and lower limits of the error bar

cpg.yca$LCLS = cpg.yca$mean_value - cpg.yca$sem_vals
cpg.yca$UCLS = cpg.yca$mean_value + cpg.yca$sem_vals
cpg.ynca$LCLS = cpg.ynca$mean_value - cpg.ynca$sem_vals
cpg.ynca$UCLS = cpg.ynca$mean_value + cpg.ynca$sem_vals
cpg.yct$LCLS = cpg.yct$mean_value - cpg.yct$sem_vals
cpg.yct$UCLS = cpg.yct$mean_value + cpg.yct$sem_vals
cpg.ynct$LCLS = cpg.ynct$mean_value - cpg.ynct$sem_vals
cpg.ynct$UCLS = cpg.ynct$mean_value + cpg.ynct$sem_vals

#total contains dataframes of subsets 
total <- rbind(cpg.yca,cpg.ynca,cpg.yct,cpg.ynct)
#errorbar function
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sem(x)
  ymax <- m+sem(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
if(special_folder== "meme"){
  graph_title= paste(nice_name,"Alpha Score MEME")
}
if(special_folder== "fubar"){
  graph_title= paste(nice_name,"Alpha Score FUBAR")
}
truenamepng = paste("output/alpha_graphs/",al1$Virus,".png",sep="")
png(truenamepng, width = 6.75, height = 6.75, units = "in", res= 300)
ggplot(aes(factor(graphit4),X.alpha., color=factor(graphit4)),data=total)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title = graph_title)+
  geom_jitter(data=cpg.ynca,aes(x=factor(graphit4)),position = position_jitter(width = .2), alpha = 0.8)+
  geom_jitter(data=cpg.yca,aes(x=factor(graphit4)),position = position_jitter(width = .2), alpha = 0.8)+
  geom_jitter(data=cpg.ynct,aes(x=factor(graphit4)),position = position_jitter(width = .2), alpha = 0.8)+
  geom_jitter(data=cpg.yct,aes(x=factor(graphit4)),position = position_jitter(width = .2), alpha = 0.8)+
  scale_x_discrete(labels = c("a \n No CpG ","a \n CpG ","t \n No CpG ","t \n CpG "))+
  scale_color_manual(labels = c("No CpG a","CpG a","No CpG t", "CpG t"), values = c("darkmagenta", "forestgreen","darkred", "dodgerblue4")) +
  #labels X and Y axis
  labs(x="CpG Type", y="Alpha Value",col=" ")+
  annotation_logticks(sides="l") +
  stat_summary(fun.data = data_summary, geom = "errorbar",width=.4, size = 1.2)+
  stat_summary(fun.data= data_summary, size = 1)
dev.off()

}


print(wilcox.test(cpg.yc$X.alpha., cpg.ync$X.alpha., alternative='less'))


array1 = datamonkey$X.alpha.[datamonkey$potential_CpG == "yes" & datamonkey$makesCpG == 0]
array2 = datamonkey$X.alpha.[datamonkey$potential_CpG == "yes" & datamonkey$makesCpG == 1]

print(wilcox.test(array1, array2, alternative='greater'))


