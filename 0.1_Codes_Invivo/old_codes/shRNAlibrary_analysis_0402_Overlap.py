#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 19:53:55 2019

@author: yolandatiao
"""

########## shRNA library analysis ##########
# Author: Huitian (Yolanda) Diao
# Apr. 2nd, 2019

########## Import ##########
import csv
import glob
import os
from astropy.io import ascii
from astropy.table import Table, join, vstack

########## Main ##########
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190402_Overlap"
os.chdir(wk_dir)

###----- Dataset: filtered Exp35&Exp56 
###----- Overlap top 30 hits in two different batches of experiments
Exp35_Input_avg = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190402_Exp35_rmOutlier/GateComparisons/InputMinusAvg_byGene_Exp35_flt.csv"
Exp35_Q3_other = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190402_Exp35_rmOutlier/GateComparisons/Q3minusOther_byGene_Exp35_flt.csv"
Exp35_Q4_Q1 = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190402_Exp35_rmOutlier/GateComparisons/Q4minusQ1_byGene_Exp35_flt.csv"

Exp56_Input_avg = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190402_Exp56_rmOutlier/GateComparisons/InputMinusAvg_byGene_Exp56_flt.csv"
Exp56_Q3_other = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190402_Exp56_rmOutlier/GateComparisons/Q3minusOther_byGene_Exp56_flt.csv"
Exp56_Q4_Q1 = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190402_Exp56_rmOutlier/GateComparisons/Q4minusQ1_byGene_Exp56_flt.csv"


Exp35_Input_avg_tab = ascii.read(Exp35_Input_avg)
Exp35_Q3_other_tab = ascii.read(Exp35_Q3_other)
Exp35_Q4_Q1_tab = ascii.read(Exp35_Q4_Q1)

Exp56_Input_avg_tab = ascii.read(Exp56_Input_avg)
Exp56_Q3_other_tab = ascii.read(Exp56_Q3_other)
Exp56_Q4_Q1_tab = ascii.read(Exp56_Q4_Q1)

Exp35_Input_avg_tab.sort('nbPctgShift')
Exp35_Q3_other_tab.sort('nbPctgShift')
Exp35_Q4_Q1_tab.sort('nbPctgShift')

Exp56_Input_avg_tab.sort('nbPctgShift')
Exp56_Q3_other_tab.sort('nbPctgShift')
Exp56_Q4_Q1_tab.sort('nbPctgShift')

Exp56_Input_avg_tab[0]

Input_avg_set = list(set(list(Exp35_Input_avg_tab['Gene'][0:50]) +
                     list(Exp56_Input_avg_tab['Gene'][0:50]) ))
Input_avg_list = list(set(list(Exp35_Input_avg_tab['Gene'][0:50]))) + list(set(list(Exp56_Input_avg_tab['Gene'][0:50]))) 
len(Input_avg_set)
len(Input_avg_list)


Q3_other_set = list(set(list(Exp35_Q3_other_tab['Gene'][0:50]) +
                     list(Exp56_Q3_other_tab['Gene'][0:50]) ))
Q3_other_list = list(set(list(Exp35_Q3_other_tab['Gene'][0:50]))) + list(set(list(Exp56_Q3_other_tab['Gene'][0:50]))) 
len(Q3_other_set)
len(Q3_other_list)

