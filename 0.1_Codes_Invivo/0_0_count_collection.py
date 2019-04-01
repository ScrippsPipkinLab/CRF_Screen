#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 21:41:41 2019

@author: yolandatiao
"""

########## Count data collection ##########
# Author: Huitian (Yolanda) Diao
# Jan 22nd, 2019
# Input files:
# *.count

########## Import ##########
import os
import csv
import glob
from astropy.io import ascii

########## Self-defined functions ##########
def count_shRNA(inFile, refList):
    #refList = shRNA_list
    #inFile = "/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/out/BC_1.R_2014_12_18_12_20_30_user_SN2-12-Exp56_45pM.out"
    count_dict = dict((ele, 0) for ele in refList)
    with open(inFile, "r") as fin:
        rfin = csv.reader(fin, delimiter = "\t")
        lastID = ""
        for row in rfin:
            currID = row[0]
            rowshRNA = row[1]
            if currID != lastID:
                if rowshRNA in shRNA_list:
                    count_dict[rowshRNA] = count_dict[rowshRNA] + 1
                else:
                    print("Error!!! %s not found in reference!!!"%rowshRNA)
            lastID = currID
    out_name = inFile.replace(".out", "_count.csv")
    with open(out_name, "w") as fout:
        wfout = csv.writer(fout, delimiter = ",")
        wfout.writerow(["shRNA", "count"])
        for key, value in count_dict.items():
            wfout.writerow([key, value])
            
    
    

########## Main ##########
ref_file = "/Volumes/Yolanda/CRF_Screen/InVivo/shRNA_ref.csv"
ref_tab = ascii.read(ref_file)
shRNA_list = list(ref_tab.columns[1])
len(shRNA_list)

##### Count collection
'''
### Exp35 batch 1
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/out/firstrun"
os.chdir(wk_dir)
for file in glob.glob("*.out"):
    count_shRNA(file, shRNA_list)

### Exp35 batch 2
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/out/secondrun_2014-12_replicates"
os.chdir(wk_dir)
for file in glob.glob("*.out"):
    count_shRNA(file, shRNA_list)
'''

### Exp 56
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/0_out/Exp56"
os.chdir(wk_dir)
for file in glob.glob("*.out"):
    count_shRNA(file, shRNA_list)


