Virus_info =read.csv("data/CpG_List.csv")
for(i in 1:nrow(Virus_info)){
  
  viruplace = paste('data/fasta/', Virus_info$name[i], sep="")
  name = as.character(Virus_info$name[i])
  splitname<-unlist(strsplit(as.character(Virus_info$name[i]),".fasta"))
  truename<-splitname[1]
  DF=Tables(truename)
  Pvalues=Wilcox_test(DF, truename)
  makeTable(Pvalues, truename)
  
}
DF=Tables(truename)
Pvalues=Wilcox_test(DF, truename)
makeTable(Pvalues, truename)