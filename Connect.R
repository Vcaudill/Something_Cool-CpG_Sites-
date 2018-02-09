# your floder where all pages are saved
setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")
source("MeaFreq.R")


DF<-meanFreq("DataSet")
DF$wtnt<-as.character(DF$wtnt)

source("WTAA.R")
DF<-getWTAA(DF)

source("MUTAA.R")
DF<-getMUTAA(DF)

source("Drastic_AA_Change.R")
DF<-big_aa_change(DF)

source("SynNonSyn.R")
DF<-functionSynNonSyn(DF)

source("CPG_Function.R")
DF<-CPG_site(DF)


#########
########
# Plots Not on Github Yet
# source("Plot1.R")
# LvsF_CpG_Printer(DF)
# 
# source("PLot2.R")
# plotsyn(DF)
# 
# source("Plot3.R")
# Fig3(DF)

