#slimdir <- "/Users/ryanwinstead/Something_Cool-CpG_Sites-"
simdata <-list.files(paste (getwd(),"/data/Slim", sep = ""))
options(warn=-1)
options(scipen=999)
simdata <-list.files(paste (getwd(),"/data/Slim", sep = ""))
POSITIONS <- 1000
TRIALS <- 100


getSlimData <- function(simdata, TRIALS, POSITIONS){
  DF = data.frame(matrix(NA, nrow = POSITIONS, ncol = TRIALS))
  #reads data into a string of 0's and 1's and creates DF 
  #Each DF colomn = 1 file. 
  col <- 1
  for (file in simdata){
    file<- paste (getwd(),"/data/Slim/", file, sep = "")
    bigstring <-readChar(file, file.info(file)$size)
    bigstrings <- gsub("\n", "", bigstring)
    
    positions<-c()
    for (i in bigstrings){
      positions<- c(positions,strsplit(i, "")[[1]])
    }
    seqvec<-c()
    #1000 mutations in vector seqvec of 0's and 1's
    for (i in 0:POSITIONS){
      seqvec<- c(seqvec, positions[i])
    }
    print(file)
    #append data to DF **Not working**
    DF[,col]<- seqvec
    
    col <- col + 1
  }
  
  return(DF)
}
#making new column of frequencies
interpret <- function(DF, TRIALS){
  row<-0
  DF$Freq<-0
  for (row in 1:nrow(DF)){
    DF$Freq[row]<-sum(as.numeric(DF[row,]))/TRIALS
    #DF$Freq[row]<-DF$Freq[row]/ncol(DF)
    }
  return(DF)
}


#Defined Dataframe Shape. 100 trials. 1000 positions

#d[,i+3] = seqvector

DF = getSlimData(simdata, TRIALS, POSITIONS)
interpretedDF = interpret(DF, TRIALS)



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