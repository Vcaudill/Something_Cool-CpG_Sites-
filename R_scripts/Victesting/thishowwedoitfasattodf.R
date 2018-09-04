library("Biostrings")

fastaFile <- readDNAStringSet("~/Desktop/env_all.fasta.txt")
description = (names(fastaFile))
id = names(fastaFile)
gene = "env"
sequence = paste(fastaFile)
df <- data.frame(description, id, gene, sequence)
#strsplit('A1.dgstdfdg', split("."))
sort.df <- with(df,  df[order(description) , ])
write.csv(sort.df, file = "Env_All_ali_try2.csv")

library(cowsay)
say("Hello world!")
someone_say_my_fortune()

require(MASS)
library(animation)

mkmovie = TRUE  #guarantees some form of output

param <- c(50-30i,18+8i,12-10i,-14-60i,1+20i)
parar <- param * exp(1i*pi/2)  #rotate by 90 degrees
pinky <- function() {
  Cx <- as.complex(rep(0,length(param)))
  Cy <- as.complex(rep(0,length(param))) 
  tv <- seq(0,2*pi,length=1000)
  
  for (i in 1:2) { #movie frames
    Cx[1] <- parar[1] + Im(param[1])
    Cx[2] <- parar[2] + Im(param[2])
    Cx[3] <- Re(param[3])
    Cx[4] <- Re(param[5]) - (i-1)
    Cx[5] <- Re(param[4])
    
    Cy[1] <- param[1] - Re(param[1]) + Im(param[4])
    Cy[2] <- param[2] - Re(param[2])
    Cy[3] <- param[3] - Re(param[3])
    
    x <- c(fourier(tv, Cx))
    y <- c(fourier(tv, Cy))
    
    plot(y, -x, type="l", col='red', lwd=10, axes=FALSE, ylab='', xlab='')
    lines(y, -x, type="l", col='pink', lwd=4)
    if (i > 1) points(Im(param[5]), Im(param[5]), col='black', pch=126, cex=2)
    else points(Im(param[5]), Im(param[5]), col='black', pch=20, cex=2)
  }
}


fourier <- function(tt,cc) {
  wt <- rep(0, length(tt)) 
  fsum <- function(n) {
    if (n > 0) wt <- wt + fsum(n-1) + Re(cc[n]) * cos(n*tt) + Im(cc[n]) * sin(n*tt)
    return(wt)
  }
  fsum(length(cc))
}


if (mkmovie) {
  aopt = ani.options(interval = 0, nmax = 2)
  saveMovie(pinky(), interval = 0.25, outdir = getwd(), width = 400, height = 400)
  ani.options(aopt)
} else pinky()


library(fortunes)
fortune()
fortune(204)
fortune("memory")

someone_say_hello <- function() {
  animal <- sample(names(animals), 1)
  say(paste("Hello, I'm a ", animal, ".", collapse = ""), by = animal)
}

someone_say_my_fortune <- function(x) {
  animal <- animal <- sample(names(animals), 1)
  say(paste(fortune(), collapse = "\n"), by = animal)
}
someone_say_hello()
someone_say_my_fortune()
