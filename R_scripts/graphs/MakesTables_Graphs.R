Virus_info =read.csv("data/CpG_List.csv")
for(i in 1:nrow(Virus_info)){
  
  viruplace = paste('data/fasta/', Virus_info$name[i], sep="")
  name = as.character(Virus_info$name[i])
  splitname<-unlist(strsplit(as.character(Virus_info$name[i]),".fasta"))
  truename<-splitname[1]
  #source("R_scripts/Tables/HowToMakeWilcoxTables.R")
  DF=Tables(truename)
  Pvalues=Wilcox_test(DF, truename)#get error x must be numeric
  makeTable(Pvalues, truename)
  source("R_scripts/graphs/redoplot.R")
  redo<-comparing_CpG_Syn_Nonsyn_new(truename)
  
  #bring back the image and set up so put two graphs in one page
  truenamepng = paste("output/GT/",truename,".png",sep="")
  png(truenamepng, width = 6.75, height = 6.75, units = "in", res= 300)
  par(mfrow=c(1,2),mar=c(3.9, 4.1, 2.1, 0.8))
  return(redo)
  dev.off()
  

  
  }
DF=Tables(truename)
Pvalues=Wilcox_test(DF, truename)
makeTable(Pvalues, truename)