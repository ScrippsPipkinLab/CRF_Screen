#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 18:56:14 2019

@author: yolandatiao
"""

########## Cytoscape dataprep ##########
# Author: Huitian (Yolanda) Diao
# May 17th, 2019

########## Import ##########
import os
import csv
from astropy.io import ascii
from astropy.table import Table, join, vstack
import math

########## Self-defined functions ##########
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
    


########## Main ##########

wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein"
os.chdir(wk_dir)

in_file = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/CRF_target_biogrid_source.csv"
in_tab = ascii.read(in_file)

###----- Delete self interactions
del_rows = []
for x in range(0, len(in_tab)):
    if in_tab[x][0] == in_tab[x][1]:
        del_rows.append(x)
in_tab.remove_rows(del_rows)

###----- Count numbers of appearance
in_tab = prune_tree(in_tab, 5)
in_tab = prune_tree(in_tab, 5)
in_tab = prune_tree(in_tab, 5)
in_tab = prune_tree(in_tab, 5)

len(set(list(in_tab['source']) + list(in_tab['target'])))

###---- Recount
  
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

