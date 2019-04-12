library(seqinr)
Virus_info <- read.csv("data/MEME_CpG_list.csv")
k=2
i=1

for(k in 1:nrow(Virus_info)){
name = as.character(Virus_info$name[k])
splitname<-unlist(strsplit(as.character(Virus_info$name[k]),".fasta"))
truename<-splitname[1]  
path_file<-paste("Hyphy/Consensus_Hyphy/",truename,".csv", sep="")
al<-read.csv(path_file)

for (i in 1: length(Virus_info)){
folder<-Virus_info$Hyphy_folders[i]
path<-paste("Hyphy/",folder,"/datamonkey-table.csv")
datamonkey<-read.csv()
file<-datamonkey$wtnt_consensus
nuc<-splitseq(file, frame = 0, word = 3)
datamonkey$wtnt_codon<-nuc[i]

}
for(j in 1:length(Virus_info)){

#list.files(path = ".", pattern = NULL, all.files = FALSE,
#           full.names = FALSE, recursive = FALSE,
#           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

#dir(path = ".", pattern = NULL, all.files = FALSE,
#    full.names = FALSE, recursive = FALSE,
#    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

#list.dirs(path = ".", full.names = TRUE, recursive = TRUE
foldername<-list.files("Hyphy/", pattern = "stuff", all.files = TRUE)
Virus_info$Hyphy_folders<-c(foldername)
new_pathfile<-paste("Hyphy/Consensus_Hyphy/",truename,".csv", sep="")
}  #folder EnteroCVP1, InfluenzaBvirus_Ha
  }
