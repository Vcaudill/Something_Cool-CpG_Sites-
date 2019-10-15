slimdir <- "/Users/ryanwinstead/Something_Cool-CpG_Sites-/data/Slim"
setwd(slimdir)
options(warn=-1)
options(scipen=999)
simdata <-list.files(slimdir)
POSITIONS <- 1000
TRIALS <- 100


getSlimData <- function(simdata){
  #reads data into a string of 0's and 1's and creates DF 
  #Each DF colomn = 1 file. 
  col <- 1
  for (file in simdata){
    bigstring <-readChar(file, file.info(file)$size)
    bigstrings <- gsub("\n", "", bigstring)
    
    positions<-c()
    for (i in bigstrings){
      positions<- c(positions,strsplit(i, "")[[1]])
    }
    seqvec<-c()
    #1000 mutations in vector seqvec of 0's and 1's
    for (i in 0:1000){
      seqvec<- c(seqvec, positions[i])
    }
    print(length(seqvec))
    #append data to DF **Not working**
    DF[,col]<- seqvec
    
    col <- col + 1
  }
  
  #return(DF)
}
#making new column of frequencies
interpret <- function(DF){
  freqs<-c()
  col<-0
  row<-0
  for (row in 1:nrow(DF)){
    ones<-0
    for (col in 1:ncol(DF)){
      print(DF[row, col])
      print(row)
      print(col)
      print("______")
      if (col=="1"){
        ones <- ones + 1
      }
      freqs<- c(freqs, (ones/TRIALS))
    }
    break
  }
  DF$Freq <- freqs
  return(DF)
}


#Defined Dataframe Shape. 100 trials. 1000 positions
DF = data.frame(matrix(NA, nrow = POSITIONS, ncol = TRIALS))
#d[,i+3] = seqvector

getSlimData(simdata)
interpret(DF)



'''


for (file in simdata){
  bigstring <-readChar(file, file.info(file)$size)
bigstrings <- gsub("\n", "", bigstring)

positions<-c()
for (i in bigstrings){
positions<- c(positions,strsplit(i, "")[[1]])

}
print(length(positions))
}



for (file in simdata){ 
  classes <- c()
  bigstring <-readChar(file, file.info(file)$size)
  bigstrings <- gsub("\n", " ", bigstring)
  for (i in strsplit(bigstrings," ")){classes<-c(classes, i)}
  positions <- c()
  for (i in classes){
    if (as.numeric(i) > 1 || is.na(as.numeric(i))){next()}
    positions<- c(positions, as.numeric(i))
  }
  ughhhh <- c()
  for (i in positions){
  i<- i * 1000
  i <- as.integer(i)
  ughhhh<-c(ughhhh, i)
  }
  CpGseq <- ""
  for (i in 0:999){
  if (i %in% ughhhh){
  CpGseq <- paste(CpGseq, "1", sep=",", collapse = ",")
  }
  else{CpGseq <- paste(CpGseq, "0",sep=",", collapse = ",")}
  }
  
  seqvector <- c(unlist(strsplit(classes,split = '')))
  seqvector <-seqvector[seqvector != ""]
  #a vector of comma sep 0s and 1s
  print(col)
  DF[,col]<- seqvector
  col <- col+1
}



'''
 #bigstring <-readChar(file, file.info(file)$size)
#bigstrings <- gsub("\n", " ", bigstring)
#positions <- as.list(strsplit(bigstrings, " "))
#slimlist <- as.list(strsplit(gsub("\n", " ", readChar(file, file.info(file)$size)), " "))