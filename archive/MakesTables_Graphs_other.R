Virus_info =read.csv("data/CpG_List_RyanVersion.csv")
for(i in 1:nrow(Virus_info)){
  name = as.character(Virus_info[i,1])
  print(name)
}
for(i in 50:nrow(Virus_info)){
  i = 51
  viruplace = paste('data/fasta/', Virus_info[i,1], sep="")
  name = as.character(Virus_info[i,1])
  splitname<-unlist(strsplit(as.character(Virus_info[i,1]),".fasta"))
  truename<-splitname[1]
  source("R_scripts/Tables/HowToMakeWilcoxTables.R")
  DF=Tables(truename)
  Pvalues=Wilcox_test(DF, truename)#get error x must be numeric
  makeTable(Pvalues, truename, Virus_info$nice_name[i])
  source("R_scripts/graphs/redoplot.R")
  comparing_CpG_Syn_Nonsyn_new(truename)
  dev.off()
  library(png)
  library(grid)
  library(gridExtra)
  img1 <-  rasterGrob(as.raster(readPNG("output/Redoplot/BkpolyomaVirus_VP1.png")), interpolate = FALSE)
  img2 <-  rasterGrob(as.raster(readPNG("output/tables_blue/BKpolyomaVirusVP1.pdf")), interpolate = FALSE)
  grid.arrange(img1, img2, ncol = 2)
  
}

library(gridExtra)
library(png)
img<-readPNG("output/redeploy/BKpolyomavirus_VP1.png")
h<-dim(img)[1]
w<-dim(img)[2]

#open new file for output
png("out.png", width=w, height=h)
par(mar=c(0,0,0,0), xpd=NA, mgp=c(0,0,0), oma=c(0,0,0,0), ann=F)
plot.new()
plot.window(0:1, 0:1)

#fill plot with image
usr<-par("usr")    
rasterImage(img, usr[1], usr[3], usr[2], usr[4])

#add text
text(.5,.5, "hello", cex=5, col=rgb(.2,.2,.2,.7))

#close image
dev.off()


#bring back the image and set up so put two graphs in one page
truenamepng = paste("output/GT/",truename,".png",sep="")
png(truenamepng, width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(1,2),mar=c(3.9, 4.1, 2.1, 0.8))
return(redo)
dev.off()

DF=Tables(truename)
Pvalues=Wilcox_test(DF, truename)
makeTable(Pvalues, truename)