library("Biostrings")
library(ggplot2)
library(scales)
library(plotrix)
library(sfsmisc)
nucfreq = function(){
  setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-/allVirus/")
  allVir = c()
  allVir = c(list.files(path = "C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-/allVirus"))
  #print(allVir)
  #processFile(allVir[1])
  #for(filepath in allVir){
  #setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-")
  #print(filepath)
  
  #processFile(filepath)
}
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
  
  #if(pur == TRUE){
  #   if(first == "A"){
  #    Acount = Acount + 1
  #   }else{
  #    Gcount =Gcount + 1}
  #}else{
  # if(first== "T"){
  #  Tcount = Tcount +1
  # }else{Ccount = Ccount +1}
  #}
}
processFile = function(filepath){
  con = file(filepath, "r")
  print(filepath)
  pur = FALSE
  stop = FALSE
  purine = c( 'A', 'G')
  pyrimidine = c('T', 'C')
  skip = c('<', '>', '-', ' ')
  comb = c("AA", "AG", "AT", "AC", "GA", "GG", "GT", "GC", "TA", "TG", "TT", "TC", "CA", "CG", "CT", "CC" )
  nuc = readChar(con, n = 1)
  x = 0
  print(nuc)
  if(stop == TRUE){
    x = x+1
    if(nuc == "p"){stop = FALSE}
    next
  }else{
    if (x %% 2 == 0){
      x = x+1
      if (nuc %in% purine){
        first = nuc
        if(x == 1){second = ""}
      }
      if (nuc %in% pyrimidine){
        first = nuc
        if(x == 1){second = ""}
      }
      if(nuc == ">"){stop = TRUE}else{next}
    }
    
    if (x %% 2 != 1){
      if(first != ""){next
      }else{
        x= x+1
        break}
      x = x+1
      print("ln 36")
      if (nuc %in% purine){
        second = nuc
        pur = TRUE
      }
      if (nuc %in% pyrimidine){
        second = nuc
        pur = FALSE
      }
      if(nuc == '>'){
        stop = TRUE
        first = ""
        second = ""
        if(nuc == "p"){stop = FALSE}
        next
      }
      pair = paste(first, second, sep = "")
      sort(pair)
    }
    
    
    close(con)
    print("bye")
    #print(Acount, Gcount, Tcount, Ccount)
  }
}



allVir = c(list.files(path = "C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-/allVirus/"))
for (file in allVir){
  setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-/allVirus/")
  fastaFile <- readDNAStringSet(file)
  name0 = split(name, ".txt", sep="")
  name = paste(file, ".csv", sep="")
  print(name)
  sequence = paste(fastaFile)
  df <- data.frame(sequence)
  #sort.df<- with(df[order(description) , ])
  setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-/csvVir")
  write.csv(df, file= name)  
}



setwd("C:/Users/ryanw/Desktop/codeLab/CpG/Something_Cool-CpG_Sites-")
df = read.csv("Counts.csv", header = TRUE)
ggplot(df, aes(x=i, y= )) + geom_point()
plot(df$i, df$mean_AA, col ='red', ylim = c(1, 2000), pch=19, log = 'y' )
lines(c(df$mean_AA), col = "red")
points(df$i, df$mean_AT, col ='red4',pch=19)
points(df$i, df$mean_AG, col ='firebrick',pch=19)
points(df$i, df$mean_AC, col ='tomato',pch=19)
points(df$i, df$mean_TA, col ='orange',pch=19)
points(df$i, df$mean_TT, col ='orange3',pch=19)
points(df$i, df$mean_TG, col ='yellow',pch=19)
points(df$i, df$mean_TC, col ='gold',pch=19)
points(df$i, df$mean_GA, col ='dodgerblue',pch=19)
points(df$i, df$mean_GT, col ='skyblue',pch=19)
points(df$i, df$mean_GG, col ='deepskyblue4',pch=19)
points(df$i, df$mean_GC, col ='blue1',pch=19)
points(df$i, df$mean_CA, col ='forestgreen',pch=19)
points(df$i, df$mean_CT, col ='darkgreen',pch=19)
points(df$i, df$mean_CG, col ='green',pch=19)
points(df$i, df$mean_CC, col ='limegreen',pch=19)
