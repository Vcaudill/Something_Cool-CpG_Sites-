# your floder where all pages are saved
Tables = function(truename){  
  setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-")
  
  #How to load data
  #virusname = 'HumanBocavirus1_NS1.fasta_pruned.mu.trim05.txt'
  #splitname<-unlist(strsplit(virusname,".fasta"))
  #truename<-splitname[1]
  
  truenameRda= paste(truename, ".Rda", sep="")
  setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-/Rda_Files")
  load(truenameRda)
  
  
  setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-")
  # Wilcox test
  source("RyanWilcox.R")
  Wilcox_test(DF, truename)
}
namelist=c("DengueVirus1", "DengueVirus2", "DengueVirus3", "DengueVirus4", "HumanBocavirus1_NS1", "HumanBocavirus1_VP1", "humanparainfluenzavirus1_F", "InfluenzaAvirus_HA_H1N1","InfluenzaAvirus_HA_H3N2", "EnterovirusA_VP1", "EnterovirusA_VP2","EnterovirusB_VP1", "EnterovirusB_VP2","EnterovirusC_VP1","EnterovirusD_VP1", "humanparainfluenzavirus1_HN", "InfluenzaAvirus_NA_H1N1", "InfluenzaBvirus_NA", "BKpolyomavirus_VP1")
for(i in namelist){
  Tables(i)
}
########
########


