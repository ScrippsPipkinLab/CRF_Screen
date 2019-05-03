#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 19:48:57 2019

@author: yolandatiao
"""

########## shRNA library analysis ##########
# Author: Huitian (Yolanda) Diao
# Apr. 2nd, 2019

# Input files:
# - csv format
# - 3 columns: type, shRNA, count
#    type: target / control
#    shRNA: genename.1
#    count: integer

########## Import ##########
import glob
import os
from astropy.io import ascii
import csv



########## Self-defined functions ##########
def avg_ctrl(inFile):
    #inFile = "Q4minusQ1_byGene.csv"
    outName = inFile.replace(".csv", "_ctrl_avg.csv")
    inTab = ascii.read(inFile)
    del inTab['Pool']
    
    # Slice table into 'duplicate' and 'unique'
    gene_list = list(inTab['Gene'])
    gene_uniq = list(set(gene_list))
    gene_dup = [x for x in gene_uniq if gene_list.count(x) > 1]
    gene_dup_anno = []
    for i in gene_list:
        if i in gene_dup:
            gene_dup_anno.append("dup")
        else:
            gene_dup_anno.append("uniq")
    inTab["dup"] = gene_dup_anno
    inTab_by_dup = inTab.group_by("dup")
    
    if inTab_by_dup.groups[0]["dup"][0] == "dup":
        dupTab = inTab_by_dup.groups[0]
        uniqTab = inTab_by_dup.groups[1]
    else:
        dupTab = inTab_by_dup.groups[1]
        uniqTab = inTab_by_dup.groups[0]
    
    # Calculate average percentage for 'duplicate'
    dupTab_by_gene = dupTab.group_by("Gene")
    dup_gene_names = [str(x).replace("\n","").replace("Gene","").replace(" ","").replace("-","") for x in dupTab_by_gene.groups.keys]
    dup_n = len(dup_gene_names)
    dup_gene_names
    for x in range(0, dup_n):
        dup_key = dup_gene_names[x]
        dup_vals = list(dupTab_by_gene.groups[x]['nbPctgShift'])
        dup_vals = [float(x) for x in dup_vals]
        dup_mean = sum(dup_vals)/len(dup_vals)
        uniqTab.add_row([dup_key, dup_mean, "avg_dup"])
    
    del uniqTab["dup"]
    ascii.write(uniqTab, outName, format="csv", overwrite=True)
    


########## Main ##########
# Calculate mean of duplicated controls for each file
os.chdir("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190403_Exp35Exp56_nbPctgToAll/GateComparisons")
for file in glob.glob("*byGene.csv"):
    avg_ctrl(file)



















