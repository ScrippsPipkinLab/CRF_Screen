#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 11:54:28 2019

@author: yolandatiao
"""

########## Networkx pagerank ##########
# Author: Huitian (Yolanda) Diao
# May 17th, 2019

########## Import ##########
import networkx as nx
import os
import csv
from astropy.io import ascii
from astropy.table import Table, join, vstack
import math

########## Self-defined functions ##########
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/PageRank"
os.chdir(wk_dir)

#####---------- Remove duplicate interactions in interaction table
biogrid_int = "/Volumes/Yolanda/CRF_Screen/Ref/BIOGRID.use_simp.txt"
int_tab = ascii.read(biogrid_int)
# Convert all genename to uppercase
int_tab["gene_name_upper"] = [x.upper() for x in list(int_tab.columns[0])]
int_tab["gene_name_upper_2"] = [x.upper() for x in list(int_tab.columns[1])]
# Remove duplicated interactions
gene1 = list(int_tab["gene_name_upper"])
gene2 = list(int_tab["gene_name_upper_2"])
int_all = [[x,y] for index, (x, y) in enumerate(zip(gene1, gene2))]
for x in range(0, len(int_all)):
    int_x = int_all[x]
    int_x.sort()
    int_all[x] = int_x
int_all = ["%s::%s"%(x[0], x[1]) for x in int_all]
int_all = list(set(int_all))
int_all = [x.split("::") for x in int_all]

out_file = "/Volumes/Yolanda/CRF_Screen/Ref/BIOGRID.use_simp.dupR.csv"
with open(out_file, "w") as fout:
    wfout = csv.writer(fout, delimiter = ",")
    wfout.writerow(["gene1", "gene2"])
    for x in int_all:
        wfout.writerow(x)
        
#####---------- Remove duplicates in protein abundance table
wt_protein_abundance = "/Volumes/Yolanda/RNAseq_Compilation/CRM_analysis_2/9_proteomics/0_input/WT_whole.csv"
p_tab = ascii.read(wt_protein_abundance)
p_colnames = p_tab.colnames
p_tab_by_name = p_tab.group_by(p_colnames[0])

out_file = "/Volumes/Yolanda/RNAseq_Compilation/CRM_analysis_2/9_proteomics/0_input/WT_whole_dupr-max.csv"
with open(out_file, "w") as fout:
    wfout = csv.writer(fout, delimiter=",")
    wfout.writerow(['gene_name','0h.a', '0h.b', '2h.a', '2h.b', '8h.a', '8h.b', '16h.a', '16h.b'])
    for i in range(0, len(p_tab_by_name.groups.keys)):
        new_row = []
        gene_tab = p_tab_by_name.groups[i]
        if len(gene_tab) > 1:
            print(gene_tab)
            new_row.append(gene_tab[0][0])
            for x in range(1, len(p_colnames)):
                new_row.append(max(list(gene_tab.columns[x])))
        else:
            new_row=gene_tab[0]
        wfout.writerow(new_row)

#####---------- Calculate average protein abundance in duplicate removed table
in_file = "/Volumes/Yolanda/RNAseq_Compilation/CRM_analysis_2/9_proteomics/0_input/WT_whole_dupr-max.csv"
out_file = "/Volumes/Yolanda/RNAseq_Compilation/CRM_analysis_2/9_proteomics/0_input/WT_whole_dupr-max_avg.csv"
in_tab = ascii.read(in_file)
in_tab.colnames
in_tab["0h"] = [(x+y)/2 for index, (x,y) in enumerate(zip(list(in_tab['0h.a']), list(in_tab['0h.b'])))]
in_tab["2h"] = [(x+y)/2 for index, (x,y) in enumerate(zip(list(in_tab['2h.a']), list(in_tab['2h.b'])))]
in_tab["8h"] = [(x+y)/2 for index, (x,y) in enumerate(zip(list(in_tab['8h.a']), list(in_tab['8h.b'])))]
in_tab["16h"] = [(x+y)/2 for index, (x,y) in enumerate(zip(list(in_tab['16h.a']), list(in_tab['16h.b'])))]
del in_tab['0h.a']
del in_tab['0h.b']
del in_tab['2h.a']
del in_tab['2h.b']
del in_tab['8h.a']
del in_tab['8h.b']
del in_tab['16h.a']
del in_tab['16h.b']
ascii.write(in_tab, out_file, format="csv", overwrite=True)

        
        
#####---------- Read protein abundance & interaction files and merge
wt_protein_abundance = "/Volumes/Yolanda/RNAseq_Compilation/CRM_analysis_2/9_proteomics/0_input/WT_whole_dupr-max_avg.csv"
dupr_int = "/Volumes/Yolanda/CRF_Screen/Ref/BIOGRID.use_simp.dupR.csv"
p_tab = ascii.read(wt_protein_abundance)
dupr_int_tab = ascii.read(dupr_int)
dupr_int_tab.colnames

# Convert all genename to uppercase
gnames = p_tab.columns[0]
del p_tab[p_tab.colnames[0]]
p_tab["gene1"] = [x.upper() for x in list(gnames)]

# Calculate log10 tpm for all columns in protein expression table
"""
tps = p_tab.colnames[:-1]
for tp in tps:
    tp_col_list = list(p_tab[tp])
    tp_col_new = [math.log10(x)/10 for x in tp_col_list]
    del p_tab[tp]
    p_tab[tp] = tp_col_new
"""

# Merge table
dupr_int_tab = join(dupr_int_tab, p_tab, keys="gene1", join_type='inner')
for x in range(0, len(p_tab.colnames)):
    p_tab.columns[x].name = p_tab.colnames[x].replace("a", "a_2").replace("b", "b_2").replace("gene1", "gene2")
dupr_int_tab = join(dupr_int_tab, p_tab, keys="gene2", join_type='inner')
dupr_int_tab.colnames

ascii.write(dupr_int_tab, "Biogrid-interaction_CRF-ProteinAbundance.csv", format="csv", overwrite=True)


# Filter self-interactions
in_file = "Biogrid-interaction_CRF-ProteinAbundance.csv"
out_file = "Biogrid-interaction_CRF-ProteinAbundance_fltself.csv"
with open(in_file, "r") as fin:
    with open(out_file, "w") as fout:
        rfin = csv.reader(fin, delimiter=",")
        wfout = csv.writer(fout, delimiter=",")
        header = next(rfin)
        wfout.writerow(header)
        for row in rfin:
            if row[0] != row[1]:
                wfout.writerow(row)

