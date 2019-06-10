#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 13:10:28 2019

@author: yolandatiao
"""

############################## Find the most likely complex for input genes ##############################
# Author: Huitian (Yolanda) Diao
# June 10th, 2019

############################## Import ##############################
import os
from astropy.io import ascii
import csv

############################## Main ##############################
def findComplex(inList):
    #inList = ["BRPF3", "BRPF1", "KAT6A", "KAT6B", "PHF17", "ING4","BRD1", "KAT7", "ING5"]
    
    refFile = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/2_GO_terms/CRM_complexes_count.csv"
    temp1 = "temp.csv"
    with open(refFile, "r") as fin:
        with open(temp1, "w") as fout:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            header = next(rfin)
            wfout.writerow(header)
            for row in rfin:
                if row[0].upper() in inList:
                    wfout.writerow([row[0].upper()] + row[1:])
    tempTab = ascii.read(temp1)
    for i in inList:
        if i not in list(tempTab["gene_name"]):
            newrow = [i] + ["None" for x in range(1, len(tempTab.colnames))]
            tempTab.add_row(newrow)
    topComplex = []
    topNum = 0

    for i in tempTab.colnames:
        if "Yes" not in list(tempTab[i]) and i != "gene_name":
            del tempTab[i]
        else:
            yes_n = list(tempTab[i]).count("Yes")
            if yes_n > topNum:
                topComplex = []
                topComplex.append(i)
                topNum = yes_n
            elif yes_n == topNum:
                topComplex.append(i)
    for i in tempTab.colnames:
        if i not in topComplex and i != "gene_name":
            del tempTab[i]
    return(tempTab)

def ListComplexes(inGeneLists, outname):
    with open(outname, "w") as fout:
        wfout = csv.writer(fout, delimiter=",")
        for listx in inGeneLists:
            tabx = findComplex(listx)
            complexes = tabx.colnames[1:]
            wfout.writerow(["; ".join(complexes)])
            wfout.writerow(["", "pctg"] + list(tabx["gene_name"]))
            for complex_x in complexes:
                complex_count = list(tabx[complex_x]).count("Yes")
                complex_pctg = float(complex_count)/len(list(tabx[complex_x]))*100
                wfout.writerow([complex_x] + [complex_pctg] + list(tabx[complex_x]))
                print ("----------------")
                print(listx)
                print("%s percentage: %s"%(complex_x, complex_pctg))
            wfout.writerow([""])
        
        
        

wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/3_CRM_protein_HGSScore/2_findComplexes"
os.chdir(wk_dir) 

g1 = ["BRPF3", "BRPF1", "KAT6A", "KAT6B", "PHF17", "ING4","BRD1", "KAT7", "ING5"]
g2 = ["SUX12", "EZH1", "PHF1", "PHF19", "JARID2", "MTF2", "EZH2","EED2"]
g3 = ["CBX2", "CBX6", "SCML2", "SCMH1", "PHC1", "CBX8", "PCGF2", "RNF2",
      "BMI1", "PCGF1", "PHC2", "RING1", "PCGF5", "KDM2B"]
g4 = ["PRMT5", "SETD7", "WHSC1", "PRDM9", "SIRT2"]
g5 = ["L3MBTL4", "L3MBTL3", "L3MBTL1"]
g6 = ["HDAC1", "HDAC2"]

geneLists = [g1, g2, g3, g4, g5, g6]
ListComplexes(geneLists, "findComplexes.csv")

    
    