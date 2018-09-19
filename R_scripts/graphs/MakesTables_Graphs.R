Virus_info =read.csv("data/CpG_List.csv")
for(i in 1:nrow(Virus_info)){
  
  viruplace = paste('data/fasta/', Virus_info$name[i], sep="")
  name = as.character(Virus_info$name[i])
  splitname<-unlist(strsplit(as.character(Virus_info$name[i]),".fasta"))
  truename<-splitname[1]
  DF=Tables(truename)
  Pvalues=Wilcox_test(DF, truename)
  makeTable(Pvalues, truename)
  source("R_scripts/graphs/redoplot.R")
  comparing_CpG_Syn_Nonsyn_new(truename)
  #bring back the image and set up so put two graphs in one page
  truenamepng = paste("output/GT/",truename,".png",sep="")
  png(truenamepng, width = 6.75, height = 6.75, units = "in", res= 300)
  
  }
DF=Tables(truename)
Pvalues=Wilcox_test(DF, truename)
makeTable(Pvalues, truename)