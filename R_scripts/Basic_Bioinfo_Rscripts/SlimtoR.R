slimdir <- "~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/data/Slim"
setwd(slimdir)
simdata <-list.files(slimdir)
for (file in simdata){ 
  DataSet <-read.delim(file, sep =" ")
  break
}