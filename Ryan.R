# your floder where all pages are saved
sort = function(pair){
  if(pair == "AA"){AA = c(pair)}
  if(pair == "GG"){GG = c(pair)}
  if(pair == "TT"){TT = c(pair)}
  if(pair == "CC"){CC = c(pair)}
  if(pair == "AT"|| pair == "TA"){AT = c(pair)}
  if(pair == "AG"|| pair == "GA"){AG = c(pair)}
  if(pair == "AC"|| pair == "CA"){AC = c(pair)}
  if(pair == "TG"|| pair == "GT"){TG = c(pair)}
  if(pair == "TC"|| pair == "CT"){TG = c(pair)}
  if(pair == "GC"|| pair == "CG"){TG = c(pair)}
}
processFile = function(filepath){
  con = file(filepath, "r")
  purine = c( 'A', 'G')
  pyrimidine = c('T', 'C')
  skip = c('<', '>', '-', ' ')
  comb = c("AA", "AG", "AT", "AC", "GA", "GG", "GT", "GC", "TA", "TG", "TT", "TC", "CA", "CG", "CT", "CC" )
  nuc = readChar(con, n = 1)
  x = 0
    while (x %% 2 == 0){
      if (nuc %in% purine){
        first = nuc
        if(n == 1){second = ""}
        }
      if (nuc %in% pyrimidine){
        first = nuc
        if(n == 1){second = ""}}
      if(nuc %in% skip){next}
    }
  
  while (x %% 2 != 1){
      if (nuc %in% purine){
        second = nuc
        pur = TRUE
        }
      if (nuc %in% pyrimidine){
        second = nuc
        pur = FALSE
      }
      if(nuc %in% skip){
        first = ""
        second = ""
        next
      }
  }
  
  if(pur == TRUE){
    if(first == "A"){
      Acount = Acount + 1
      }else{
        Gcount =Gcount + 1}
    }else{
      if(first== "T"){
        Tcount = Tcount +1
        }else{Ccount = Ccount +1}
    }
  pair = paste(first, second, sep = "")
  sort(pair)


  close(con)
  x = x+1
  print(Acount, Gcount, Tcount, Ccount)
  }
nucfreq = function(){
  setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-/allVirus/")
  allVir = c()
  allVir = c(list.files(path = "C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-/allVirus"))
  #print(allVir)
  processFile(allVir[1])
  for(filepath in allVir){
    #setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-")
    #print(filepath)
    
    #processFile(filepath)
  }
}




Tables = function(truename){  
  setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-")
  
  truenameRda= paste(truename, ".Rda", sep="")
  setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-/Rda_Files")
  load(truenameRda)
  
  
  setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-")
  # Wilcox test
  source("RyanWilcox.R")
  Wilcox_test(DF, truename)
}

  

namelist=c(DengueVirus1, DengueVirus2, DengueVirus3, DengueVirus4, humanparainfluenzavirus1_F, humanparainfluenzavirus1_HN, humanparainfluenzavirus3_HN, InfluenzaAvirus_HA_H1N1,InfluenzaAvirus_HA_H3N2, InfluenzaAvirus_NA_H1N1, InfluenzaAvirus_NA_H3N2,InfluenzaBvirus_HA, InfluenzaBvirus_NA, EnterovirusA_VP1, EnterovirusA_VP2,EnterovirusB_VP1, EnterovirusB_VP2,EnterovirusC_VP1,EnterovirusD_VP1, BK_polyomavirus_VP1, HumanBocavirus1_NS1, HumanBocavirus1_VP1)
for(i in namelist){
  Tables(i)}
  

Acount = 0
Gcount = 0
Tcount = 0
Ccount = 0
nucfreq()
########
########


