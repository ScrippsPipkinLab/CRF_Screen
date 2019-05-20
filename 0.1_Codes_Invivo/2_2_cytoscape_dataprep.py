#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 18:56:14 2019

@author: yolandatiao
"""

########################################  Cytoscape dataprep ######################################## 
# Author: Huitian (Yolanda) Diao
# May 17th, 2019

########################################  Import ######################################## 
import os
import csv
from astropy.io import ascii
from astropy.table import Table, join, vstack
import math

########################################  Self-defined functions ########################################
def prune_tree(inTab, pruneN):
    print("Table row number: %s" %len(inTab))
    source_list = list(inTab["source"])
    tgt_list = list(inTab["target"])
    total_list = source_list + tgt_list
    source_list_app = [total_list.count(i) for i in source_list]
    tgt_list_app = [total_list.count(i) for i in tgt_list]
    
    inTab['source_n'] = source_list_app
    inTab['target_n'] = tgt_list_app
    del_rows = [index for index, x in enumerate(source_list_app) if x <= pruneN]
    del_rows = del_rows + [index for index, x in enumerate(tgt_list_app) if x <= pruneN]
    inTab.remove_rows(del_rows)
    print("Pruned table row number: %s" %len(inTab))
    return(inTab)
    
def select_int_and(inFile, selectList, out_apdx):
    outFile = inFile.replace(".csv", "_%s.csv"%out_apdx)
    with open(inFile, "r") as fin:
        with open(outFile, "w") as fout:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            header = next(rfin)
            wfout.writerow(header)
            for row in rfin:
                if row[0] in selectList and row[1] in selectList:
                    wfout.writerow(row)

def select_int_or(inFile, selectList, out_apdx):
    outFile = inFile.replace(".csv", "_%s.csv"%out_apdx)
    with open(inFile, "r") as fin:
        with open(outFile, "w") as fout:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            header = next(rfin)
            wfout.writerow(header)
            for row in rfin:
                if row[0] in selectList or row[1] in selectList:
                    wfout.writerow(row)
    
def convertGN(listX, inNames, outNames):
    for i in range(0, len(listX)):
        if listX[i] in inNames:
            idx_i = inNames.index(listX[i])
            listX[i] = outNames[idx_i]
            print(listX[i])
    return(listX)
######################################## Main ########################################

##########-------------------- Create cytoscape source file with biogrid data
'''
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein"
os.chdir(wk_dir)

in_file = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/CRF_target_biogrid.csv"
in_tab = ascii.read(in_file)
'''

###----- Delete self interactions
'''
del_rows = []
for x in range(0, len(in_tab)):
    if in_tab[x][0] == in_tab[x][1]:
        del_rows.append(x)
in_tab.remove_rows(del_rows)
ascii.write(in_tab, "CRF_target_biogrid_rmself.csv", format="csv", overwrite=True)
'''

###----- Calculate appearance
'''
in_tab_all_genes = list(in_tab["gene_name"]) + list(in_tab["gene_name2"])
in_tab['source_appearance'] = [in_tab_all_genes.count(i) for i in list(in_tab["gene_name"])]
in_tab['source_appearance_sqrt'] = [math.sqrt(i) for i in list(in_tab['source_appearance'])]
ascii.write(in_tab, "CRF_target_biogrid_rmself_count.csv", format="csv", overwrite=True)
'''

###----- Count numbers of appearance
'''
in_tab = prune_tree(in_tab, 0)
in_tab = prune_tree(in_tab, 5)
in_tab = prune_tree(in_tab, 5)
in_tab = prune_tree(in_tab, 5)

len(set(list(in_tab['source']) + list(in_tab['target'])))
'''

###---- Recount
''' 
source_list = list(in_tab["source"])
tgt_list = list(in_tab["target"])
total_list = source_list + tgt_list
source_list_app = [total_list.count(i) for i in source_list]
tgt_list_app = [total_list.count(i) for i in tgt_list]
in_tab['source_n'] = source_list_app
in_tab['target_n'] = tgt_list_app
in_tab['sourcesize'] = [math.sqrt(x) for x in source_list_app]
in_tab['targetsize'] = [math.sqrt(x) for x in tgt_list_app]

ascii.write(in_tab, "CRF_target_biogrid_source_pruned.csv", format="csv", overwrite=True)
'''


##########-------------------- Select genes for plotting based on adjusted z-score
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein"
os.chdir(wk_dir)

###----- Find top and bottom quarter protein interaction for CRFs (in vivo screen)
# Read file
adj_z_file = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190516/5_zscore_div_sqrt_pval/all_z-score_div_sqrt-p_sqrt.csv"
adj_z_tab = ascii.read(adj_z_file)
# Convert gene names
alt_gn_file = "/Volumes/Yolanda/CRF_Screen/Ref/CRF_alternative_gn.csv"
alt_gn_tab = ascii.read(alt_gn_file)
adj_z_tab['gene_name'] = convertGN(list(adj_z_tab['gene_name']), list(alt_gn_tab.columns[1]), list(alt_gn_tab.columns[0]))
adj_z_tab['gene_name'] = [x.upper() for x in list(adj_z_tab['gene_name'] )]
# Different comparisons
adj_z_tab_len = len(adj_z_tab)
qt = math.floor(adj_z_tab_len / 4)
adj_z_tab.sort('Q4minusQ1')
q4_q1_dn = list(adj_z_tab['gene_name'])[0:qt]
q4_q1_up = list(adj_z_tab['gene_name'])[adj_z_tab_len - qt:]
adj_z_tab.sort('Q3minusOther')
q3_o_dn = list(adj_z_tab['gene_name'])[0:qt]
q3_o_up = list(adj_z_tab['gene_name'])[adj_z_tab_len - qt:]
adj_z_tab.sort('InputMinusAvg')
in_a_dn = list(adj_z_tab['gene_name'])[0:qt]
in_a_up = list(adj_z_tab['gene_name'])[adj_z_tab_len - qt:]

###----- Slice data
int_file = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/CRF_target_biogrid_rmself_count.csv"

select_int_and(int_file, q4_q1_dn, "q4_q1_dn")
select_int_and(int_file, q4_q1_up, "q4_q1_up")
select_int_and(int_file, q3_o_dn, "q3_o_dn")
select_int_and(int_file, q3_o_up, "q3_o_up")
select_int_and(int_file, in_a_dn, "in_a_dn")
select_int_and(int_file, in_a_up, "in_a_up")


select_int_and(int_file, q4_q1_dn, "q4_q1_dn_or")
select_int_and(int_file, q4_q1_up, "q4_q1_up_or")
select_int_and(int_file, q3_o_dn, "q3_o_dn_or")
select_int_and(int_file, q3_o_up, "q3_o_up_or")
select_int_and(int_file, in_a_dn, "in_a_dn_or")
select_int_and(int_file, in_a_up, "in_a_up_or")











