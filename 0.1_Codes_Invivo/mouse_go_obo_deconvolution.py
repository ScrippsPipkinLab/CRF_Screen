#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 16:45:37 2019

@author: yolandatiao
"""

########################################  Cytoscape dataprep ######################################## 
# Author: Huitian (Yolanda) Diao
# May 20th, 2019


########################################  Import ######################################## 
import os
import csv

########################################  Main ########################################
wk_dir = "/Volumes/Yolanda/CRF_Screen/Ref"
os.chdir(wk_dir)

in_file = "mouse_go_obo_20190509.txt"
out_file = "mouse_go_obo_20190509.csv"

in_file = "mouse_go_basic_obo_20190509.txt"
out_file = "mouse_go_basic_obo_20190509.csv"

with open(in_file, "r") as fin:
    with open(out_file, "w") as fout:
        rfin = csv.reader(fin, delimiter = "\t")
        wfout = csv.writer(fout, delimiter = ",")
        header = ["id", "namespace", "name", "def"]
        wfout.writerow(header)
        newrow = ["", "", "", ""]
        alt_id = []
        for row in rfin:
            if len(row) >= 1: 
                if row[0] == "[Term]":
                    if newrow != ["", "", "", ""]:
                        if newrow[0] != "":
                            wfout.writerow(newrow)
                            if len(alt_id) > 0:
                                for idx in alt_id:
                                    wfout.writerow([idx] + newrow[1:])
                    newrow = ["", "", "", ""]
                    alt_id = []
                elif "id: " in row[0]:
                    if "alt_id: " in row[0]:
                        alt_id.append(row[0].replace("alt_id: ", ""))
                    else:
                        newrow[0] = row[0].replace("id: ", "")                    
                elif "name: " in row[0]:
                    newrow[2] = row[0].replace("name: ", "")
                elif "namespace: " in row[0]:
                    newrow[1] = row[0].replace("namespace: ", "")
                elif "def: " in row[0]:
                    newrow[3] = row[0].replace("def: ", "")
            