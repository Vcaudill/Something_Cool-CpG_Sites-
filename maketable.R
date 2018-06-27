
makeTable <- function(Pvalues){
  pdf("rectangles.pdf", width = 7, height= 5)
  col1 <- c("A-G", "T-C")
  col2 <- c("Syn v CpG", "SynCpG v Non CpG", "NonSynCpG v NonCpG")
  ycoor <- c(4*100/5 , 3*100/5 + 4, 3*100/5 - 12, 2*100/5 -7, 1*100/5-1, 100/5-14- 4, 100)
  df = data.frame(col1, col2, Pvalues)
  
  #layout(matrix(1:1, nrow = 1))
  
  
  par(xpd=F)
  plot(1, 2, xlim=c(0,100),ylim=c(0,100), col=0, xaxt="n", yaxt="n", xlab="d", ylab="bf")
  
  abline(v = 100/5)
  abline(v = 2*100/3)
  abline(h = 100-100/7 + 3)
  abline(h = 100- 2*100/7+2)
  abline(h = 100 - 3*100/7+1)
  abline(h= 100 - 4*100/7 -1)
  abline(h= 100 -5*100/7 - 2)
  abline(h= 100 -6*100/7 - 3)
  
  text(x=100/7- 6, y= 5*100/5-3, "Mutation Type")
  text(x=3*100/7, y = 5*100/5-3, "Comparison")
  text(x= 6*100/7, y=5*100/5-3, "p-value")
  text(x= 100/7, y= 4*100/5+1, "A-G")
  text(x= 100/7, y = 2*100/5 - 6, "T-C")
  text(x= 3*100/7, y = ycoor[1], labels= col2[1])
  text(x= 3*100/7, y = ycoor[2], labels= col2[2])
  text(x= 3*100/7, y = ycoor[3], labels= col2[3])
  text(x= 3*100/7, y = ycoor[4], labels= col2[1])
  text(x= 3*100/7, y = ycoor[5], labels= col2[2])
  text(x= 3*100/7, y = ycoor[6], labels= col2[3])
  
  num <- 1
  print(Pvalues)
  for (i in Pvalues){
    print(i)
    
    
  library(scales)
  if (i < 0.05){
    a = 0.4 
  }
  else if(i <0.5){
    a = 0.3
  }
  else if(i >0.5){
    a = 0.1
  }
  
  rect(xleft = 2*100/3, xright = 200, ybottom = ycoor[num]-7, ytop = ycoor[num]+8, col = alpha("firebrick", a))
  text(x= 6*100/7, y =ycoor[num], labels = i)
  num = num + 1
  
  }
  print("end")
  
  
}
