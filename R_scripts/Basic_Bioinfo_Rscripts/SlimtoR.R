slimdir <- "/Users/ryanwinstead/Something_Cool-CpG_Sites-/data/Slim"
setwd(slimdir)

simdata <-list.files(slimdir)

for (file in simdata){ 
  classes <- c("positions:")
  positions <- as.list(strsplit(gsub("\n", " ", readChar(file, file.info(file)$size)), " "))
  for (entry in positions){
    if (grepl("^[[:digit:]]+L",entry)){
      positions[[entry]] <- NULL
      print(entry)
    }
  }
  #positions[["//"]] <- NULL
  #DF <-read.table(file, sep =" ", colClasses = classes, header=TRUE, fill = TRUE)
  DF <- positions
  print(positions)
  break
}



#bigstring <-readChar(file, file.info(file)$size)
#bigstrings <- gsub("\n", " ", bigstring)
#positions <- as.list(strsplit(bigstrings, " "))