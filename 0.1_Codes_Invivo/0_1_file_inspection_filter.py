#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 01:55:51 2019

@author: yolandatiao
"""

########## Count data collection ##########
# Author: Huitian (Yolanda) Diao
# Jan 22nd, 2019
# Input files:
# - *count.csv
# Dependencies:
# - 0_0_count_collection.py

########## Import ##########
import os
import csv
import glob
from astropy.io import ascii
from astropy.table import Table, join

########## Self-defined functions ##########
def annoCount(inFile, expList, ctrlList):
    #inFile = "/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/count/bc_1_count.csv"
    #expList = p1_shRNA
    #ctrlList = ctrl_shRNA
    
    inFile_name_noformat = inFile.split("/")[-1].replace("_count.csv", "")
    outname = inFile.replace(".csv", "_flt.csv")
    outname_else = inFile.replace(".csv", "_contaminants.csv")
    outname_sum = inFile.replace(".csv", "_summary.csv")
    
    inTab = ascii.read(inFile)
    shRNA_list = list(inTab['shRNA'])
    tgt_names = []
    ctrl_names = []
    else_names = []
    for index, shRNA in enumerate(shRNA_list):
        if shRNA in ctrlList:
            ctrl_names.append(index)
        elif shRNA in expList:
            tgt_names.append(index)
        else:
            else_names.append(index)
    tgt_tab = inTab[tgt_names]
    ctrl_tab = inTab[ctrl_names]
    else_tab = inTab[else_names]
    tgt_tab.sort(inTab.colnames[1])
    ctrl_tab.sort(inTab.colnames[1])
    else_tab.sort(inTab.colnames[1])
    tgt_tab.reverse()
    ctrl_tab.reverse()
    else_tab.reverse()
    
    tgt_sum = sum(list(tgt_tab.columns[1]))
    ctrl_sum = sum(list(ctrl_tab.columns[1]))
    else_sum = sum(list(else_tab.columns[1]))
    all_sum = tgt_sum + ctrl_sum + else_sum

    with open(outname, "w") as fout:
        with open(outname_else, "w") as fout_else:
            with open(outname_sum, "w") as fout_sum:
                wfout = csv.writer(fout)
                wfout_else = csv.writer(fout_else)
                wfout_sum = csv.writer(fout_sum)
                wfout.writerow(["type", "shRNA", "count"])
                wfout_else.writerow(["type", "shRNA", "count"])
                wfout_sum.writerow(["x", inFile_name_noformat])
                # Write outfile
                for x in range(0, len(tgt_tab)):
                    wfout.writerow(["target"] + list(tgt_tab[x]))
                for x in range(0, len(ctrl_tab)):
                    wfout.writerow(["control"] + list(ctrl_tab[x]))
                # Write else file
                for x in range(0, len(else_tab)):
                    wfout_else.writerow(["contaminants"] + list(else_tab[x]))
                # Write summary
                wfout_sum.writerow(["Total reads"] + [all_sum])
                wfout_sum.writerow(["Target percentage"] + [float(tgt_sum)/all_sum*100])
                wfout_sum.writerow(["Control percentage"] + [float(ctrl_sum)/all_sum*100])
                wfout_sum.writerow(["Contaminant percentage"] + [float(else_sum)/all_sum*100])
                wfout_sum.writerow(["No.1 Target"] + ["%s: %s"%(tgt_tab[0][0], tgt_tab[0][1])])
                wfout_sum.writerow(["No.2 Target"] + ["%s: %s"%(tgt_tab[1][0], tgt_tab[1][1])])
                wfout_sum.writerow(["No.3 Target"] + ["%s: %s"%(tgt_tab[2][0], tgt_tab[2][1])])
                wfout_sum.writerow(["Top 3 target %"] + 
                                   [float(tgt_tab[0][1] + tgt_tab[1][1] + tgt_tab[2][1])/tgt_sum*100])  
                wfout_sum.writerow(["No.1 Control"] + ["%s: %s"%(ctrl_tab[0][0], ctrl_tab[0][1])])
                wfout_sum.writerow(["No.2 Control"] + ["%s: %s"%(ctrl_tab[1][0], ctrl_tab[1][1])])
                wfout_sum.writerow(["No.3 Control"] + ["%s: %s"%(ctrl_tab[2][0], ctrl_tab[2][1])])
                wfout_sum.writerow(["Top 3 control %"] + 
                                   [float(ctrl_tab[0][1] + ctrl_tab[1][1] + ctrl_tab[2][1])/ctrl_sum*100]) 
                wfout_sum.writerow(["No.1 Contaminant"] + ["%s: %s"%(else_tab[0][0], else_tab[0][1])])
                wfout_sum.writerow(["No.2 Contaminant"] + ["%s: %s"%(else_tab[1][0], else_tab[1][1])])
                wfout_sum.writerow(["No.3 Contaminant"] + ["%s: %s"%(else_tab[2][0], else_tab[2][1])])
                wfout_sum.writerow(["Top 3 contaminant %"] + 
                                   [float(else_tab[0][1] + else_tab[1][1] + else_tab[2][1])/else_sum*100])                

        
########## Main ##########

### Read reference file
ref_file = "/Volumes/Yolanda/CRF_Screen/InVivo/shRNA_ref.csv"
ref_tab = ascii.read(ref_file)
ref_tab_pool = list(ref_tab.columns[0])
ref_tab_shRNA = list(ref_tab.columns[1])
p1_list = [index for (index, x) in enumerate(ref_tab_pool) if x == "pool1-7|bc1_3_17_20"]
p2_list = [index for (index, x) in enumerate(ref_tab_pool) if x == "pool8-14|bc6-10"]
p3_list = [index for (index, x) in enumerate(ref_tab_pool) if x == "p15-20|bc11-15"]
p1_shRNA = [ref_tab_shRNA[index] for index in p1_list]
p2_shRNA = [ref_tab_shRNA[index] for index in p2_list]
p3_shRNA = [ref_tab_shRNA[index] for index in p3_list]

ctrl_file = "/Volumes/Yolanda/CRF_Screen/InVivo/shRNA_control.csv"
ctrl_tab = ascii.read(ctrl_file)
ctrl_shRNA = list(ctrl_tab.columns[0])

##### Exp35

### Process all files
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/1_count/Exp35"
os.chdir(wk_dir)
for file in glob.glob("*.csv"):
    file_bc_n = int(file.replace("_firstrun", "").replace("_count.csv", "").replace("bc_",""))
    if ((file_bc_n >= 1) and (file_bc_n <= 5)) or ((file_bc_n == 17) or (file_bc_n == 20)):
        annoCount(file, p1_shRNA, ctrl_shRNA)
    elif (file_bc_n >= 6) and (file_bc_n <= 10):
        annoCount(file, p2_shRNA, ctrl_shRNA)
    elif (file_bc_n >= 11) and (file_bc_n <= 15):
        annoCount(file, p3_shRNA, ctrl_shRNA)


### Moved all output files to seperate folders

### Merge all summary files

wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/2_summary"
os.chdir(wk_dir)

file1 = "bc_1_count_summary.csv"
tab1 = ascii.read(file1)

for file in glob.glob("*.csv"):
    if file != file1:
        newtab = ascii.read(file)
        tab1 = join(tab1, newtab)

ascii.write(tab1, "Exp35_reads_summary.csv", format="csv", overwrite=True)

### Merge two sequencing rounds
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/2_flt"
os.chdir(wk_dir)

ref_file = "/Volumes/Yolanda/CRF_Screen/InVivo/sample_barcode_exp35.csv"
ref_tab = ascii.read(ref_file)
ref_tab_bc = list(ref_tab.columns[0])
ref_tab_group = list(ref_tab.columns[1])

file_list = []
for file in glob.glob("*.csv"):
    file_list.append(file)
file_list

file_list_bc = [int(x.replace("_count_flt.csv","").replace("_firstrun","").replace("bc_","")) for x in file_list]

ref_dict = {}
for idx, (refg, refbc) in enumerate(zip(ref_tab_group, ref_tab_bc)):
    ref_dict[refg] = [int(x) for x in refbc.split(",")]
ref_dict
    
for key, value in ref_dict.items():
    #key = 'P1-7_Input'
    #value = [1]
    
    use_files = []
    for idx, (bc, file) in enumerate(zip(file_list_bc, file_list)):
        if bc in value:
            use_files.append(file)
    tab1 = ascii.read(use_files[0])
    
    if len(use_files) == 2:
        tab1['count'].name = 'count1'
        tabx = ascii.read(use_files[1])
        del tabx['type']
        tabx['count'].name = 'count2'
        tab1 = join(tab1, tabx, keys = 'shRNA', join_type = "outer")
        count_all = [x+y for index, (x, y) in enumerate(zip(list(tab1['count1']), list(tab1['count2'])))]
        tab1['count'] = count_all
        del tab1['count1']
        del tab1['count2']
    elif len(use_files) >2:
        print("Error!! More than 3 seq rounds?")    
        
    outname = "%s_combined.csv"%key
    ascii.write(tab1, outname, format="csv", overwrite=True)
    


##### Exp56

### Process all files
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/1_count/Exp56"
os.chdir(wk_dir)
for file in glob.glob("*.csv"):
    file_bc_n = int(file.replace("_firstrun", "").replace("_count.csv", "").replace("BC_",""))
    if ((file_bc_n >= 1) and (file_bc_n <= 5)) or ((file_bc_n == 17) or (file_bc_n == 18) or (file_bc_n == 20)):
        annoCount(file, p1_shRNA, ctrl_shRNA)
    elif (file_bc_n >= 6) and (file_bc_n <= 10):
        annoCount(file, p2_shRNA, ctrl_shRNA)
    elif (file_bc_n >= 11) and (file_bc_n <= 15):
        annoCount(file, p3_shRNA, ctrl_shRNA)

### Merge all summary files
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/2_summary/Exp56"
os.chdir(wk_dir)

file1 = "BC_1_count_summary.csv"
tab1 = ascii.read(file1)

for file in glob.glob("*.csv"):
    if file != file1:
        newtab = ascii.read(file)
        tab1 = join(tab1, newtab)

ascii.write(tab1, "Exp56_reads_summary.csv", format="csv", overwrite=True)

