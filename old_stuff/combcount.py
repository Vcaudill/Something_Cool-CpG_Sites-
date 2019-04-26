import math
import csv
from os import walk
import pandas as pd
from numpy import mean
fastaFolder = "csvVir"


def comcount(fasta, name, i):
    lisAA = []
    lisAT = []
    lisAG = []
    lisAC = []
    lisTA = []
    lisTT = []
    lisTG = []
    lisTC = []
    lisGA = []
    lisGT = []
    lisGG = []
    lisGC = []
    lisCA = []
    lisCT = []
    lisCG = []
    lisCC = []

    j = 0
    while j < len(fasta):
        seq = fasta.loc[j, 'sequence']
        lisAA.append(seq.count("AA"))
        lisAT.append(seq.count("AT"))
        lisAG.append(seq.count("AG"))
        lisAC.append(seq.count("AC"))
        lisTA.append(seq.count("TA"))
        lisTT.append(seq.count("TT"))
        lisTG.append(seq.count("TG"))
        lisTC.append(seq.count("TC"))
        lisGA.append(seq.count("GA"))
        lisGT.append(seq.count("GT"))
        lisGG.append(seq.count("GG"))
        lisGC.append(seq.count("GC"))
        lisCA.append(seq.count("CA"))
        lisCT.append(seq.count("CT"))
        lisCG.append(seq.count("CG"))
        lisCC.append(seq.count("CC"))
        j += 1
    mean_AA = mean(lisAA)
    mean_AT = mean(lisAT)
    mean_AG = mean(lisAG)
    mean_AC = mean(lisAC)
    mean_TA = mean(lisTA)
    mean_TT = mean(lisTT)
    mean_TG = mean(lisTG)
    mean_TC = mean(lisTC)
    mean_GA = mean(lisGA)
    mean_GT = mean(lisGT)
    mean_GG = mean(lisGG)
    mean_GC = mean(lisGC)
    mean_CA = mean(lisCA)
    mean_CT = mean(lisCT)
    mean_CG = mean(lisCG)
    mean_CC = mean(lisCC)
    meanlist = [name, i, mean_AA, mean_AT, mean_AG, mean_AC, mean_TA, mean_TT, mean_TG,
                mean_TC, mean_GA, mean_GT, mean_GG, mean_GC, mean_CA, mean_CT, mean_CG, mean_CC]
    return meanlist


myfiles = ["DengueVirus1.csv", "DengueVirus2.csv", 'DengueVirus3.csv', 'DengueVirus4.csv', 'humanparainfluenzavirus1_F.csv', "humanparainfluenzavirus1_HN.csv", 'humanparainfluenzavirus3_HN.csv', 'InfluenzaAvirus_HA_H1N1.csv', 'InfluenzaAvirus_HA_H3N2.csv', 'InfluenzaAvirus_NA_H1N1.csv',
           'InfluenzaAvirus_NA_H3N2.csv', 'InfluenzaBvirus_HA.csv', 'InfluenzaBvirus_NA.csv', 'EnterovirusA_VP1.csv', 'EnterovirusA_VP2.csv', 'EnterovirusB_VP1.csv', 'EnterovirusB_VP2.csv', 'EnterovirusC_VP1.csv', 'EnterovirusD_VP1.csv', 'HumanBocavirus1_NS1.csv', 'HumanBocavirus1_VP1.csv']

for (dirpath, dirnames, filenames) in walk(fastaFolder):
    myfiles.extend(filenames)
    break
headerlist = ['name', 'i', 'mean_AA', 'mean_AT', 'mean_AG', 'mean_AC', 'mean_TA', 'mean_TT', 'mean_TG',
              'mean_TC', 'mean_GA', 'mean_GT', 'mean_GG', 'mean_GC', 'mean_CA', 'mean_CT', 'mean_CG', 'mean_CC']
thenum = []
print(headerlist)
# myfiles = ['myfile1.csv', "",]
i = 0
with open('/Users/victoria/Desktop/Counts.csv', 'w') as f1:

    headerlist = ['name', 'i', 'mean_AA', 'mean_AT', 'mean_AG', 'mean_AC', 'mean_TA', 'mean_TT', 'mean_TG',
                  'mean_TC', 'mean_GA', 'mean_GT', 'mean_GG', 'mean_GC', 'mean_CA', 'mean_CT', 'mean_CG', 'mean_CC']
    writer = csv.DictWriter(f1, fieldnames=headerlist)
    writer.writeheader()

    for file in myfiles:

        print(file)
        fasta = pd.read_csv(file)
        fasta.drop(fasta.columns[[0]], axis=1)

        thenums = comcount(fasta, myfiles[i], i)
        csvWri1 = csv.writer(f1, delimiter=',')
        csvWri1.writerow(thenums)
        # print(the_seq)
        i += 1
        print(i)


# namelist = [DengueVirus1.csv, DengueVirus2.csv, DengueVirus3, DengueVirus4, humanparainfluenzavirus1_F, humanparainfluenzavirus1_HN, humanparainfluenzavirus3_HN, InfluenzaAvirus_HA_H1N1, InfluenzaAvirus_HA_H3N2, InfluenzaAvirus_NA_H1N1,
#             InfluenzaAvirus_NA_H3N2, InfluenzaBvirus_HA, InfluenzaBvirus_NA, EnterovirusA_VP1, EnterovirusA_VP2, EnterovirusB_VP1, EnterovirusB_VP2, EnterovirusC_VP1, EnterovirusD_VP1, BK_polyomavirus_VP1, HumanBocavirus1_NS1, HumanBocavirus1_VP1]
