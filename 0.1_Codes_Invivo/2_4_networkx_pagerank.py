#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 11:54:28 2019

@author: yolandatiao
"""

#################### Networkx pagerank ####################
# Author: Huitian (Yolanda) Diao
# May 17th, 2019

#################### Import ####################
import networkx as nx
import os
import csv
from astropy.io import ascii
from astropy.table import Table, join, vstack
import math

#################### Self-defined functions ####################
def prune_tree(inTab, pruneN):
    print("Table row number: %s" %len(inTab))
    source_list = list(inTab["gene1"])
    tgt_list = list(inTab["gene2"])
    total_list = source_list + tgt_list
    source_list_app = [total_list.count(i) for i in source_list]
    tgt_list_app = [total_list.count(i) for i in tgt_list]
    
    del_rows = [index for index, x in enumerate(source_list_app) if x <= pruneN]
    del_rows = del_rows + [index for index, x in enumerate(tgt_list_app) if x <= pruneN]
    inTab.remove_rows(del_rows)
    print("Pruned table row number: %s" %len(inTab))
    return(inTab)

def pageRankProtein(inFile, outFile, sCol, tCol, swCol, twCol):
    #inFile = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/PageRank/Biogrid-interaction_CRF-ProteinAbundance_fltself.csv"
    #outFile = "protein_pr_0h.csv"
    #sCol = "gene1"
    #tCol = "gene2"
    #swCol = "0h_1"
    #twCol = "0h_2"
    
    in_tab = ascii.read(inFile)
    source_list = list(in_tab[sCol])
    target_list = list(in_tab[tCol])
    source_wt_list = list(in_tab[swCol])
    target_wt_list = list(in_tab[twCol])
    
    source_wt_list = [(x**(1/2)) for x in source_wt_list]
    target_wt_list = [(x**(1/2)) for x in target_wt_list]
    source_and_wt = ["%s::%s"%(x, str(y)) for index, (x,y) in enumerate(zip(source_list, source_wt_list))]
    target_and_wt = ["%s::%s"%(x, str(y)) for index, (x,y) in enumerate(zip(target_list, target_wt_list))]
    name_wt_list = list(set(source_and_wt + target_and_wt))
    p_dict = {}
    for i in name_wt_list:
        p_dict[i.split("::")[0]] = float(i.split("::")[1])    
    
    DG = nx.DiGraph()
    for x in range(0, len(in_tab)):
        s_x = source_list[x]
        t_x = target_list[x]
        #sw_x = source_wt_list[x]
        #tw_x = target_wt_list[x]
        DG.add_edges_from([(s_x, t_x), (t_x, s_x)])
    
    pr = nx.pagerank(DG, personalization=p_dict, max_iter=10000)
    with open(outFile, "w") as fout:
        wfout = csv.writer(fout, delimiter=",")
        wfout.writerow(["gene_name", "PageRankScore"])
        for key, val in pr.items():
            wfout.writerow([key,val])
    

#################### Main ####################
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/PageRank"
os.chdir(wk_dir)


###----- PageRank
in_file = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/PageRank/Biogrid-interaction_CRF-ProteinAbundance_fltself.csv"

out_file = "protein_pr_0h.csv"
source_col = "gene1"
target_col = "gene2"
source_weight = "0h_1"
target_weight = "0h_2"
pageRankProtein(in_file, out_file, source_col, target_col, source_weight, target_weight)

out_file = "protein_pr_16h.csv"
source_col = "gene1"
target_col = "gene2"
source_weight = "16h_1"
target_weight = "16h_2"
pageRankProtein(in_file, out_file, source_col, target_col, source_weight, target_weight)
