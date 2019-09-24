slimdir <- "/Users/ryanwinstead/Something_Cool-CpG_Sites-/data/Slim"
setwd(slimdir)
options(warn=-1)
options(scipen=999)


simdata <-list.files(slimdir)

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
  
}



#bigstring <-readChar(file, file.info(file)$size)
#bigstrings <- gsub("\n", " ", bigstring)
#positions <- as.list(strsplit(bigstrings, " "))
#slimlist <- as.list(strsplit(gsub("\n", " ", readChar(file, file.info(file)$size)), " "))