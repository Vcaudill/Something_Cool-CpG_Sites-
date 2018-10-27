Virus_info =read.csv("data/CpG_List.csv")
library(png)
for(i in 1:nrow(Virus_info)){
  name = as.character(Virus_info[i,1])
  print(name)
}
for(i in 1:nrow(Virus_info)){
  if (Virus_info[i,1] == "Humanherpesvirus2_gD.fasta_pruned.mu.trim05"){
    next
  }
  viruplace = paste('data/fasta/', Virus_info[i,1], sep="")
  name = as.character(Virus_info[i,1])
  splitname<-unlist(strsplit(as.character(Virus_info[i,1]),".fasta"))
  truename<-splitname[1]
  # source("R_scripts/Tables/HowToMakeWilcoxTables.R")
  # DF=Tables(truename)
  # Pvalues=Wilcox_test(DF, truename)#get error x must be numeric
  # makeTable(Pvalues, truename)
  # source("R_scripts/graphs/redoplot.R")
  # redo<-comparing_CpG_Syn_Nonsyn_new(truename)
  truenameboth= paste(truename, 'both', ".png", sep="")
  sarplot= paste('output/redeploy/', truename, ".png", sep="")
  ryaplot= paste('output/redeploy/', truename,"tables", ".png", sep="")
  isg<-readPNG(sarplot)
  isng<-readPNG(ryaplot)
  h<-dim(isg)[1]
  w<-dim(isg)[2]
  
  png(truenameboth, width=w, height=h)
  par(mfrow=c(1,2))
  # par(mar=c(0,0,0,0), xpd=NA, mgp=c(0,0,0), oma=c(0,0,0,0), ann=F)
  plot.new()
  plot.window(0:1, 0:1)
  
  #fill plot with image
  usr<-par("usr")    
  rasterImage(isg, usr[1], usr[3], usr[2], usr[4])
  
  plot.new()
  plot.window(0:1, 0:1)
  
  #fill plot with image
  usr<-par("usr")    
  rasterImage(isng, usr[1], usr[3], usr[2], usr[4])
  #add text
  
  #close image
  dev.off()
  
}

library(png)
library(magick)
img<-image_read("output/redeploy/BKpolyomavirus_VP1.png","output/redeploy/BKpolyomavirus_VP1tables.png")
isg<-readPNG("output/redeploy/BKpolyomavirus_VP1.png")
isng<-readPNG("output/redeploy/BKpolyomavirus_VP1tables.png")
ing <- pdf("output/redeploy/BKpolyomavirus_VP1.pdf")
print(img)
image_append(image_scale(img, "100"), stack = TRUE)
layout(matrix(1:2, ncol=1, byrow=TRUE))
combine(img, isg)

h<-dim(img)[1]
w<-dim(img)[2]

#open new file for output
png("out.png", width=w, height=h)
par(mfrow=c(1,2))
# par(mar=c(0,0,0,0), xpd=NA, mgp=c(0,0,0), oma=c(0,0,0,0), ann=F)
plot.new()
plot.window(0:1, 0:1)

#fill plot with image
usr<-par("usr")    
rasterImage(isg, usr[1], usr[3], usr[2], usr[4])

plot.new()
plot.window(0:1, 0:1)

#fill plot with image
usr<-par("usr")    
rasterImage(isng, usr[1], usr[3], usr[2], usr[4])
#add text

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

