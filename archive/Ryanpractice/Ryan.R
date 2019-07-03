# your floder where all pages are saved


namelist=c("DengueVirus1", "DengueVirus2", "DengueVirus3", "DengueVirus4", "humanparainfluenzavirus1_F", "humanparainfluenzavirus1_HN", "humanparainfluenzavirus3_HN", "InfluenzaAvirus_HA_H1N1","InfluenzaAvirus_HA_H3N2", "InfluenzaAvirus_NA_H1N1", "InfluenzaAvirus_NA_H3N2","InfluenzaBvirus_HA", "InfluenzaBvirus_NA", "EnterovirusA_VP1", "EnterovirusA_VP2","EnterovirusB_VP1", "EnterovirusB_VP2","EnterovirusC_VP1","EnterovirusC_VP2","EnterovirusD_VP1", "BK_polyomavirus_VP1", "HumanBocavirus1_NS1", "HumanBocavirus1_VP1")
  Tables = function(truename){  
  
    setwd("~/Desktop/Something_Cool-CpG_Sites-")
    truenameRda= paste(truename, ".Rda", sep="")
    setwd("~/Desktop/Something_Cool-CpG_Sites-/data/Rda")
    print(truenameRda)
    load(truenameRda)
    # Wilcox test
    setwd("~/Desktop/Something_Cool-CpG_Sites-")
    source("RyanWilcox.R")
    Wilcox_test(DF, truename)
  }
  
  
for(truename in namelist){
  
  Tables(truename)}
  
Tables("DengueVirus1")
truename = "DengueVirus1"
########
########


  



