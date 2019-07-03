source("R_scripts/graphs/redoplot.R")
source("R_scripts/Tables/HowToMakeWilcoxTables.R")
truename <- list('DengueVirus1', 'DengueVirus2', 'DengueVirus3', 'DengueVirus4', 'humanparainfluenzavirus1_F', 'humanparainfluenzavirus1_HN', 'humanparainfluenzavirus3_HN', 'InfluenzaAvirus_HA_H1N1','InfluenzaAvirus_HA_H3N2', 'InfluenzaAvirus_NA_H1N1', 'InfluenzaAvirus_NA_H3N2','InfluenzaBvirus_HA', 'InfluenzaBvirus_NA', 'EnterovirusA_VP1', 'EnterovirusA_VP2','EnterovirusB_VP1', 'EnterovirusB_VP2','EnterovirusC_VP1','EnterovirusC_VP2','EnterovirusD_VP1', 'HumanBocavirus1_NS1', 'HumanBocavirus1_VP1', 'HepatitisB', 'Humanherpesvirus2_glycoprotein_G', 'Humanpapillomavirus16', 'Humanpapillomavirus16_L1', 'Humanrespiratorysyncytialvirus', 'Humanrespiratorysyncytialvirus_G', 'JCpolyomavirus_VP1', 'Measles','Measles_hemagglutinin_OR_haemagglutinin','ParvovirusB19_NS1', 'ParvovirusB19_VP1', 'RhinovirusB', 'RhinovirusB_polyprotein', 'RhinovirusC', 'RotavirusA_VP6', 'humanparainfluenzavirus1', 'humanparainfluenzavirus3','BKpolyomavirus_VP1','HepatitisB_core','HepatitisB_s','HepatitisB_precore','HepatitisB_pre_S','HepatitisB_polymerase_truncated_precore','HepatitisB_polymerase')
Virus_info<- read.csv("data/CpG_List.csv")

for (i in truename){
  
  print(i)
  comparing_CpG_Syn_Nonsyn_new(i)
}

for(i in 1:nrow(Virus_info)){
  #print(Virus_info$name[i])
  name = as.character(Virus_info$name[i])
  splitname<-unlist(strsplit(as.character(Virus_info$name[i]),".fasta"))
  truename<-splitname[1]
  print(truename)
  nice_name = as.character(Virus_info$nice_name[i])
  if (truename == "Humanherpesvirus2_gD") {
    next
  }
#redo plots with titles
  comparing_CpG_Syn_Nonsyn_new(truename, nice_name)
# Ryans tables now with the nice name column
  DF=Tables(truename)
  Pvalues=Wilcox_test(DF, truename)
  makeTable(Pvalues, truename, nice_name)
}
