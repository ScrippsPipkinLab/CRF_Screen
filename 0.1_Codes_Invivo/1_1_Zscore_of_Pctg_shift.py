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
from astropy.table import Table, Column, join   # For using astropy table functions
import statistics as st



########## Self-defined functions ##########
def UniformCase(strx):
    return(strx[0].upper() + strx[1:].lower())    
    
    
def avg_ctrl(inFile):
    #inFile = "Q4minusQ1_byGene.csv"
    outName = inFile.replace(".csv", "_ctrl_avg.csv")
    inTab = ascii.read(inFile)
    inTab['Gene'] = [UniformCase(x) for x in list(inTab['Gene'])]
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
    
def ZScore(inFile):
    #inFile = "Q4minusQ1_byGene_ctrl_avg.csv"
    outFile = inFile.replace(".csv", "_Z.csv")
    intab = ascii.read(inFile)
    shift = list(intab['nbPctgShift'])
    shift_dev = st.stdev(shift)
    shift_mean = sum(shift)/len(shift)
    z_list = [(x-shift_mean)/shift_dev for x in shift]
    intab['ZScore'] = z_list
    ascii.write(intab, outFile, format="csv", overwrite=True)

########## Main ##########
# Calculate mean of duplicated controls for each file
os.chdir("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/3_Gate_comparisons_compiled")
for file in glob.glob("*byGene.csv"):
    avg_ctrl(file)

os.chdir("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/3_Gate_comparisons_compiled")
for file in glob.glob("*avg.csv"):
    ZScore(file)

#####---------- Compile percentile shift z-score files
os.chdir("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/4_gate_comparisons_combined")

i_a = ascii.read("InputMinusAvg_byGene_ctrl_avg_Z.csv")
q3_o = ascii.read("Q3minusOther_byGene_ctrl_avg_Z.csv")
q4_q1 = ascii.read("Q4minusQ1_byGene_ctrl_avg_Z.csv")

del i_a['nbPctgShift']
del q3_o['nbPctgShift']
del q4_q1['nbPctgShift']
i_a.columns[1].name = "InputMinusAvg"
q3_o.columns[1].name = "Q3minusOther"
q4_q1.columns[1].name = "Q4minusQ1"

cb_tab = join(i_a, q3_o, keys = "Gene")
cb_tab = join(cb_tab, q4_q1, keys = "Gene")
cb_tab.columns[0].name = "gene_name"
cb_tab.colnames

ascii.write(cb_tab, "all_z-score.csv", format="csv", overwrite=True)

#####---------- Adjust z score table to make p > 0.9 genes into zscore=0
os.chdir("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/4_t-test_by_gene")

i_a = ascii.read("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/3_Gate_comparisons_compiled/InputMinusAvg_byGene_ctrl_avg_Z.csv")
q3_o = ascii.read("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/3_Gate_comparisons_compiled/Q3minusOther_byGene_ctrl_avg_Z.csv")
q4_q1 = ascii.read("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/3_Gate_comparisons_compiled/Q4minusQ1_byGene_ctrl_avg_Z.csv")
i_a_p = ascii.read("InputMinusAvg_t-test.by.gene.csv")
q3_o_p = ascii.read("Q3minusOther_t-test.by.gene.csv")
q4_q1_p = ascii.read("Q4minusQ1_t-test.by.gene.csv")


del i_a['nbPctgShift']
del q3_o['nbPctgShift']
del q4_q1['nbPctgShift']
del i_a_p['avg_percentile']
del q3_o_p['avg_percentile']
del q4_q1_p['avg_percentile']
i_a_p.columns[0].name = "Gene"
q3_o_p.columns[0].name = "Gene"
q4_q1_p.columns[0].name = "Gene"

i_a = join(i_a, i_a_p, keys = "Gene")
q3_o = join(q3_o, q3_o_p, keys = "Gene")
q4_q1 = join(q4_q1, q4_q1_p, keys = "Gene")


def p_cutoff(intab, cutoff):
    for x in range(0, len(intab)):
        if intab[x][2] > cutoff:
            intab[x][1] = 0
    return(intab)

i_a = p_cutoff(i_a, 0.6)
q3_o = p_cutoff(q3_o, 0.6)
q4_q1 = p_cutoff(q4_q1, 0.6)



i_a.columns[1].name = "InputMinusAvg_z"
q3_o.columns[1].name = "Q3minusOther_z"
q4_q1.columns[1].name = "Q4minusQ1_z"
i_a.columns[2].name = "InputMinusAvg_p"
q3_o.columns[2].name = "Q3minusOther_p"
q4_q1.columns[2].name = "Q4minusQ1_p"

cb_tab = join(i_a, q3_o, keys = "Gene")
cb_tab = join(cb_tab, q4_q1, keys = "Gene")
cb_tab.columns[0].name = "gene_name"
cb_tab.colnames

ascii.write(cb_tab, "all_z-score_p.csv", format="csv", overwrite=True)

















