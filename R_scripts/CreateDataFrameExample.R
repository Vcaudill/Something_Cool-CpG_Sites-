# freq calc v2 GL

library(seqinr)
library(stringi)

bk <- read.fasta("bk.txt")
#bk<-read.fasta("InfluenzaAvirus_HA_H1N1.txt")
#bk<-read.fasta("DengueVirus1.txt")
nightcrewBK= function(data) {
  # dataframe columns
  #In bk data 1089 is the amount of observations it will be different for ever virus
  num <- c(1:1089)
  WTnt <- c()
  MeanFreq <- c()
  #ours 1089
  # for MeanFreq calculation later
  absfreq <- c(rep(0, 1089))
  totalcount <- c(rep(0, 1089))
  
  # average WT calculation
  # counts number of each nucleotide in each position
  acount <- c(rep(0, 1089))
  gcount <- c(rep(0, 1089))
  ccount <- c(rep(0, 1089))
  tcount <- c(rep(0, 1089))
  nuc <- c()
  
  # same as line 20 comment
  for (i in 1:length(bk)) {
    sequence <- bk[[i]]
    for (j in 1:length(sequence)) {
      if (sequence[j] == 'a') {
        acount[j] = acount[j] + 1
      }
      if (sequence[j] == 'g') {
        gcount[j] = gcount[j] + 1
      }
      if (sequence[j] == 'c') {
        ccount[j] = ccount[j] + 1
      }
      if (sequence[j] == 't') {
        tcount[j] = tcount[j] + 1
      }
    }
  }
  
  # assigns wtnt based on most frequent nucleotide across all sequences per position
  for (j in 1:length(sequence)) {
    nuc[j] <- max(c(acount[j], gcount[j], ccount[j], tcount[j]))
    if (max(nuc[j]) == acount[j]) {
      WTnt[j] <- 'a'
    }
    if (max(nuc[j]) == gcount[j]) {
      WTnt[j] <- 'g'
    }
    if (max(nuc[j]) == ccount[j]) {
      WTnt[j] <- 'c'
    }
    if (max(nuc[j]) == tcount[j]) {
      WTnt[j] <- 't'
    }
    nuc <- c()
  }
  
  # gives absolute totals to be used for frequency calculation
  for (i in 1:length(bk)) {
    sequence <- bk[[i]]
    for (j in 1:length(sequence)) {
      if (WTnt[j] == 'a') {
        if (sequence[j] == 'g') {
          absfreq[j] <- absfreq[j] + 1
          totalcount[j] <- totalcount[j] + 1
        }
        else if (sequence[j] == 'a') {
          totalcount[j] <- totalcount[j] + 1
        }
      }
      if (WTnt[j] == 'g') {
        if (sequence[j] == 'a') {
          absfreq[j] <- absfreq[j] + 1
          totalcount[j] <- totalcount[j] + 1
        }
        else if (sequence[j] == 'g') {
          totalcount[j] <- totalcount[j] + 1
        }
      }
      if (WTnt[j] == 'c') {
        if (sequence[j] == 't') {
          absfreq[j] <- absfreq[j] + 1
          totalcount[j] <- totalcount[j] + 1
        }
        else if (sequence[j] == 'c') {
          totalcount[j] <- totalcount[j] + 1
        }
      }
      if (WTnt[j] == 't') {
        if (sequence[j] == 'c') {
          absfreq[j] <- absfreq[j] + 1
          totalcount[j] <- totalcount[j] + 1
        }
        else if (sequence[j] == 't') {
          totalcount[j] <- totalcount[j] + 1
        }
      }
    }
  }
  
  # calculates frequency as percentage
  for (i in 1:length(absfreq)) {
    MeanFreq[i] <- absfreq[i] / totalcount[i]
  }
  # translation and comparison setup
  TypeOfSite <- c()
  MUTAA <- c()
  WTAAp <- seqinr::translate(WTnt)
  
  #fig out true WTAA
  WTAAs <- stri_dup(WTAAp, 3)
  WTAA <- unlist(strsplit(WTAAs, ""))
  
  WTAA <- as.character(WTAA)
  # creates dataframe containing all data
  bk_dataY <- data.frame(num, WTnt, MeanFreq, WTAA)
  
}

bk_data<-nightcrewBK(bk)
MUTAA= function(bk_data){
  x=1
  for (x in 1:nrow(bk_data)) {
    if (bk_data$WTnt[x] == "a") {
      bk_data$A[x] <- "g"
    }
    if (bk_data$WTnt[x] == "g") {
      bk_data$A[x] <- "a"
    }
    if (bk_data$WTnt[x] == "c") {
      bk_data$A[x] <- "t"
    }
    if (bk_data$WTnt[x] == "t") {
      bk_data$A[x] <- "c"
    }
    if (bk_data$WTnt[x] == "a") {
      bk_data$A[x] <- "g"
    }
    if (bk_data$WTnt[x] == "g") {
      bk_data$A[x] <- "a"
    }
    if (bk_data$WTnt[x] == "c") {
      bk_data$A[x] <- "t"
    }
    if (bk_data$WTnt[x] == "t") {
      bk_data$A[x] <- "c"
    } 
  }
  
  x=1
  y=1
  count<-1
  bk_data$MA<-c(0)
  bk_data$MUTAA<-c(0)
  for (x in 1:(nrow(bk_data)/3)) {
    for(y in 1:3){
      if(y==1){
        bk_data$MUTAA[count]<-translate(r<-c(bk_data$A[count],as.character(bk_data$WTnt[count+1]),as.character(bk_data$WTnt[count+2])))
      }
      if(y==2){
        bk_data$MUTAA[count]<-seqinr::translate(a<-c(as.character(bk_data$WTnt[count-1]),as.character(bk_data$A[count]),as.character(bk_data$WTnt[count+1])))
      }
      if(y==3){
        bk_data$MUTAA[count]<-seqinr::translate(p<-c(as.character(bk_data$WTnt[count-2]),as.character(bk_data$WTnt[count-1]),as.character(bk_data$A[count])))
      }
      count<-count+1
    }
  }
  bk_data <-subset(bk_data, select = -c(MA, A))
  return(bk_data)
  
}
bk_data<-MUTAA(bk_data)

write.csv(bk_data, "bk_data.csv")
save(bk_data, file = "bk_data.Rda")
#load("bk_data.Rda")
