#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 2 11:22:56 2019

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
from __future__ import print_function
import csv
import glob
import os
from astropy.io import ascii
from astropy.table import Table, join, vstack
from scipy.special import gammaln
from scipy.special import psi
from scipy.misc import factorial
from scipy.optimize import fmin_l_bfgs_b as optim
import sys
from sklearn.preprocessing import quantile_transform
import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt

########## Self defined functions ##########
def fit_nbinom(X, initial_params=None):
    # Copyright (C) 2014 Gokcen Eraslan
    # https://github.com/gokceneraslan/fit_nbinom/blob/master/fit_nbinom.py
    # X is a numpy array representing the data
    # initial params is a numpy array representing the initial values of
    # size and prob parameters
    infinitesimal = np.finfo(np.float).eps

    def log_likelihood(params, *args):
        r, p = params
        X = args[0]
        N = X.size

        #MLE estimate based on the formula on Wikipedia:
        # http://en.wikipedia.org/wiki/Negative_binomial_distribution#Maximum_likelihood_estimation
        result = np.sum(gammaln(X + r)) \
            - np.sum(np.log(factorial(X))) \
            - N*(gammaln(r)) \
            + N*r*np.log(p) \
            + np.sum(X*np.log(1-(p if p < 1 else 1-infinitesimal)))

        return -result

    def log_likelihood_deriv(params, *args):
        r, p = params
        X = args[0]
        N = X.size

        pderiv = (N*r)/p - np.sum(X)/(1-(p if p < 1 else 1-infinitesimal))
        rderiv = np.sum(psi(X + r)) \
            - N*psi(r) \
            + N*np.log(p)

        return np.array([-rderiv, -pderiv])

    if initial_params is None:
        #reasonable initial values (from fitdistr function in R)
        m = np.mean(X)
        v = np.var(X)
        size = (m**2)/(v-m) if v > m else 10

        #convert mu/size parameterization to prob/size
        p0 = size / ((size+m) if size+m != 0 else 1)
        r0 = size
        initial_params = np.array([r0, p0])

    bounds = [(infinitesimal, None), (infinitesimal, 1)]
    optimres = optim(log_likelihood,
                     x0=initial_params,
                     #fprime=log_likelihood_deriv,
                     args=(X,),
                     approx_grad=1,
                     bounds=bounds)

    params = optimres[0]
    return {'size': params[0], 'prob': params[1]}

def pctgTotal(inFile):
    #inFile = "/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_flt/P1-7_Input_flt.csv"
    #print(inFile)
    inFileName = inFile.split("/")[-1].replace(".csv", "")
    outFileName = inFileName + "_pctg.csv"
    inTab = ascii.read(inFile)
    countList = inTab.columns[len(inTab.colnames)-1]
    countSum = sum(countList)
    pctgList = [float(x)/countSum*100 for x in countList]
    inTab['pctg'] = pctgList
    ascii.write(inTab, outFileName, format="csv", overwrite="True")


def ZScore(inFile):
    #inFile = "/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_flt/P1-7_Input_flt.csv"
    #print(inFile)
    inFileName = inFile.split("/")[-1].replace(".csv", "")
    outFileName = inFileName + "_ZScore.csv"
    inTab = ascii.read(inFile)
    countList = inTab.columns[2]
    countSum = sum(countList)
    countAvg = countSum/len(countList)
    countStd = np.std(np.array(countList))
    ZList = [((float(x) - countAvg)/countStd) for x in countList]
    inTab['ZScore'] = ZList
    ascii.write(inTab, outFileName, format="csv", overwrite="True")

def fltOutlier(inFile):
    inFileName = inFile.split("/")[-1].replace("_ZScore.csv", "")
    outFileName = inFileName + "_fltOutlier.csv"
    with open(inFile, "r") as fin:
        with open(outFileName, "w") as fout:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            header = next(rfin)
            newheader = header[0:-1]
            wfout.writerow(newheader)
            for row in rfin:
                if ((float(row[3]) <= 2.5) and (float(row[3]) >= -2.5)):
                    wfout.writerow(row[0:-1])
        
def nbPctg(inFile):
    #inFile = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190401/P1-7_Input_flt_pctg.csv"
    print(inFile)
    inFileName = inFile.split("/")[-1].replace(".csv", "")
    outFileName = inFileName + "_nbPctg.csv"
    inTab = ascii.read(inFile)
    pctgList = inTab.columns[3]
    pctgList_million = [x*1000000 for x in pctgList]
    allNb = fit_nbinom(np.array(pctgList_million))
    nbPctgList = [st.nbinom.cdf(x, allNb['size'], allNb['prob'])*100 for x in pctgList_million]    
    inTab['nbPctg'] = nbPctgList
    ascii.write(inTab, outFileName, format="csv", overwrite="True")  

def nbPctgTotal(inFile, allNb):
    #inFile = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190401/P1-7_Input_flt_pctg.csv"
    print(inFile)
    inFileName = inFile.split("/")[-1].replace(".csv", "")
    outFileName = inFileName + "_nbPctgTotal.csv"
    inTab = ascii.read(inFile)
    pctgList = inTab.columns[len(inTab.colnames)-1]
    pctgList_million = [x*1000000 for x in pctgList]
    nbPctgList = [st.nbinom.cdf(x, allNb['size'], allNb['prob'])*100 for x in pctgList_million]    
    inTab['nbPctg'] = nbPctgList
    ascii.write(inTab, outFileName, format="csv", overwrite="True")  
    
def nbRawCount(inFile):
    #inFile = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190401/P1-7_Input_flt_pctg.csv"
    print(inFile)
    inFileName = inFile.split("/")[-1].replace(".csv", "")
    outFileName = inFileName + "_nbRawCount.csv"
    inTab = ascii.read(inFile)
    countList = inTab.columns[2]
    allNb = fit_nbinom(np.array(countList))
    nbRawList = [st.nbinom.cdf(x, allNb['size'], allNb['prob'])*100 for x in countList]    
    inTab['nbRawCount'] = nbRawList
    ascii.write(inTab, outFileName, format="csv", overwrite="True")     

########## Main ##########
#----- Directory
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190513_Exp35Exp56_nbPctl-All"
os.chdir(wk_dir)

#----- Calculate Z-Score, filter outliers
for file in glob.glob("/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/2_flt/Exp56/*.csv"):
    ZScore(file)
for file in glob.glob("/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/2_flt/Exp35/combined/*.csv"):
    ZScore(file)    
for file in glob.glob("*ZScore.csv"):
    fltOutlier(file)

#----- Add up raw reads from Exp35 and Exp56
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190513_Exp35Exp56_nbPctl-All/0_fltOutlier"
os.chdir(wk_dir)

filelist = []
for file in glob.glob("*fltOutlier.csv"):
    filelist.append(file)
typelist = [x.replace("_combined", "").replace("_flt", "").replace("Outlier.csv","") for x in filelist]
typelist = list(set(typelist))

type_dict = {}
for typex in typelist:
    filex = [x for x in filelist if typex in x]
    type_dict[typex] = filex

for typex, filex in type_dict.items():
    #typex = 'P1-7_Q4'
    #filex = ['P1-7_Q4_flt_fltOutlier.csv', 'P1-7_Q4_combined_fltOutlier.csv']
    outname = typex + "_Exp35Exp56.csv"
    tab1 = ascii.read(filex[0])
    tab2 = ascii.read(filex[1])
    del tab1["type"]
    del tab2["type"]
    tab1["count"].name = "count1"
    tab2["count"].name = "count2"
    newtab = join(tab1, tab2, keys="shRNA", join_type="inner")
    total_count = [x+y for index,(x,y) in enumerate(zip(list(newtab["count1"]), list(newtab["count2"])))]
    newtab["count"] = total_count
    del newtab["count1"]
    del newtab["count2"]
    ascii.write(newtab, outname, format="csv", overwrite=True)

#----- Filter out bench contaminants
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190513_Exp35Exp56_nbPctl-All/1_Exp35Exp56_combined"
os.chdir(wk_dir)
# Contaminant: CDK9, Ccnt1
# Off target shRNA: Cd19
cont_list = ["CDK9", "Ccnt1", "Cd19"]
def flt_ct(inFile):
    outName = inFile.replace(".csv", "flt-ct.csv")
    with open(inFile, 'r') as fin:
        with open(outName, 'w') as fout:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            for row in rfin:
                if row[0].split(".")[0] not in cont_list:
                    wfout.writerow(row)
for file in glob.glob("*_Exp35Exp56.csv"):
    flt_ct(file)


#----- Calculate percentages of counts in each groups
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190513_Exp35Exp56_nbPctl-All/2_flt_comtaminants"
os.chdir(wk_dir)
for file in glob.glob("*Exp35Exp56flt-ct.csv"):
    pctgTotal(file)

#----- Calculate percentiles of count percentage in total distribution
pctgmillion_list = []
for file in glob.glob("*ct_pctg.csv"):
    tab = ascii.read(file)
    pctgmillion_list += [x*1000000 for x in list(tab["pctg"])]
pctgmillion_list.sort()
plt.plot(pctgmillion_list)
plt.show
all_nb = fit_nbinom(np.array(pctgmillion_list))

for file in glob.glob("*Exp35Exp56flt-ct_pctg.csv"):
    nbPctgTotal(file, all_nb)
    
   
###--- Gate comparisons
#----- Q1vQ4
def Q4minusQ1(q4File, q1File):
    #q4File = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190401/P1-7_Q4_Exp35Exp56_pctg_nbPctgTotal.csv"
    #q1File = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190401/P1-7_Q1_Exp35Exp56_pctg_nbPctgTotal.csv"
    group = q4File.split("/")[-1].split("_")[0]
    outName = group + "_Q4minusQ1.csv"
    q4Tab = ascii.read(q4File)
    q1Tab = ascii.read(q1File)
    del q4Tab['count']
    del q4Tab['pctg']
    del q1Tab['count']
    del q1Tab['pctg']
    q4Tab["nbPctg"].name = "nbPctg_Q4"
    q1Tab["nbPctg"].name = "nbPctg_Q1"
    allTab = join(q4Tab, q1Tab, join_type="inner", keys="shRNA")
    q4minusq1 = [x-y for index, (x,y) in enumerate(zip(list(allTab["nbPctg_Q4"]), list(allTab["nbPctg_Q1"])))]
    allTab["q4minusq1_nbPctg"] = q4minusq1
    ascii.write(allTab, outName, format="csv", overwrite=True)
    

P17Q4 = "P1-7_Q4_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P17Q1 = "P1-7_Q1_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
Q4minusQ1(P17Q4, P17Q1)

P814Q4 = "P8-14_Q4_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P814Q1 = "P8-14_Q1_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
Q4minusQ1(P814Q4, P814Q1)

P1521Q4 = "P15-21_Q4_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P1521Q1 = "P15-21_Q1_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
Q4minusQ1(P1521Q4, P1521Q1)



def Q3minusOther(q4File, q3File, q2File, q1File):
    #q4File = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190401/P1-7_Q4_Exp35Exp56_pctg_nbPctgTotal.csv"
    #q1File = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190401/P1-7_Q1_Exp35Exp56_pctg_nbPctgTotal.csv"
    group = q4File.split("/")[-1].split("_")[0]
    outName = group + "_Q3minusOther.csv"
    q4Tab = ascii.read(q4File)
    q3Tab = ascii.read(q3File)
    q2Tab = ascii.read(q2File)
    q1Tab = ascii.read(q1File)
    del q4Tab['count']
    del q3Tab['count']
    del q2Tab['count']
    del q1Tab['count']
    del q4Tab['pctg']
    del q3Tab['pctg']
    del q2Tab['pctg']
    del q1Tab['pctg']
    q4Tab["nbPctg"].name = "nbPctg_Q4"
    q3Tab["nbPctg"].name = "nbPctg_Q3"
    q2Tab["nbPctg"].name = "nbPctg_Q2"
    q1Tab["nbPctg"].name = "nbPctg_Q1"
    allTab = join(q4Tab, q3Tab, join_type="inner", keys="shRNA")
    allTab = join(allTab, q2Tab, join_type="inner", keys="shRNA")
    allTab = join(allTab, q1Tab, join_type="inner", keys="shRNA")
    avg124 = [(x+y+z)/3 for index, (x,y,z) in enumerate(zip(list(allTab["nbPctg_Q1"]),list(allTab["nbPctg_Q2"]),list(allTab["nbPctg_Q4"])))]
    q3minusOther = [x-y for index, (x,y) in enumerate(zip(list(allTab["nbPctg_Q3"]), avg124))]
    allTab["q3minusOther_nbPctg"] = q3minusOther
    ascii.write(allTab, outName, format="csv", overwrite=True)

P17Q4 = "P1-7_Q4_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P17Q3 = "P1-7_Q3_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P17Q2 = "P1-7_Q2_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P17Q1 = "P1-7_Q1_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
Q3minusOther(P17Q4, P17Q3, P17Q2, P17Q1)

P814Q4 = "P8-14_Q4_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P814Q3 = "P8-14_Q3_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P814Q2 = "P8-14_Q2_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P814Q1 = "P8-14_Q1_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
Q3minusOther(P814Q4, P814Q3, P814Q2, P814Q1)

P1521Q4 = "P15-21_Q4_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P1521Q3 = "P15-21_Q3_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P1521Q2 = "P15-21_Q2_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P1521Q1 = "P15-21_Q1_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
Q3minusOther(P1521Q4, P1521Q3, P1521Q2, P1521Q1)


def InputVsRest(q4File, q3File, q2File, q1File,inputFile):
    #q4File = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190401/P1-7_Q4_Exp35Exp56_pctg_nbPctgTotal.csv"
    #q1File = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190401/P1-7_Q1_Exp35Exp56_pctg_nbPctgTotal.csv"
    group = q4File.split("/")[-1].split("_")[0]
    outName = group + "_InputMinusAvg.csv"
    q4Tab = ascii.read(q4File)
    q3Tab = ascii.read(q3File)
    q2Tab = ascii.read(q2File)
    q1Tab = ascii.read(q1File)
    inputTab = ascii.read(inputFile)
    del q4Tab['count']
    del q3Tab['count']
    del q2Tab['count']
    del q1Tab['count']
    del inputTab['count']
    del q4Tab['pctg']
    del q3Tab['pctg']
    del q2Tab['pctg']
    del q1Tab['pctg']
    del inputTab['pctg']
    q4Tab["nbPctg"].name = "nbPctg_Q4"
    q3Tab["nbPctg"].name = "nbPctg_Q3"
    q2Tab["nbPctg"].name = "nbPctg_Q2"
    q1Tab["nbPctg"].name = "nbPctg_Q1"
    inputTab["nbPctg"].name = "nbPctg_input"
    allTab = join(q4Tab, q3Tab, join_type="inner", keys="shRNA")
    allTab = join(allTab, q2Tab, join_type="inner", keys="shRNA")
    allTab = join(allTab, q1Tab, join_type="inner", keys="shRNA")
    allTab = join(allTab, inputTab, join_type="inner", keys="shRNA")
    avgAll = [(a+b+c+d)/4 for index, (a,b,c,d) in enumerate(zip(list(allTab["nbPctg_Q1"]),list(allTab["nbPctg_Q2"]),list(allTab["nbPctg_Q3"]),list(allTab["nbPctg_Q4"])))]
    InminusAvg = [x-y for index, (x,y) in enumerate(zip(list(allTab["nbPctg_input"]), avgAll))]
    allTab["inputMinusAvg_nbPctg"] = InminusAvg
    ascii.write(allTab, outName, format="csv", overwrite=True)

P17Q4 = "P1-7_Q4_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P17Q3 = "P1-7_Q3_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P17Q2 = "P1-7_Q2_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P17Q1 = "P1-7_Q1_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P17input = "P1-7_Input_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
InputVsRest(P17Q4, P17Q3, P17Q2, P17Q1, P17input)

P814Q4 = "P8-14_Q4_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P814Q3 = "P8-14_Q3_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P814Q2 = "P8-14_Q2_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P814Q1 = "P8-14_Q1_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P814input = "P8-14_Input_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
InputVsRest(P814Q4, P814Q3, P814Q2, P814Q1, P814input)

P1521Q4 = "P15-21_Q4_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P1521Q3 = "P15-21_Q3_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P1521Q2 = "P15-21_Q2_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P1521Q1 = "P15-21_Q1_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
P1521input = "P15-21_Input_Exp35Exp56flt-ct_pctg_nbPctgTotal.csv" 
InputVsRest(P1521Q4, P1521Q3, P1521Q2, P1521Q1, P1521input)


###--- Average effect by pool (calculate only the last column)
def avgByGene(inFile):
    #inFile = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190401/GateComparisons/P1-7_InputMinusAvg.csv"
    outFile = inFile.replace(".csv", "_byGene.csv")
    inTab = ascii.read(inFile)
    inTabCols = inTab.colnames
    for x in range(1, len(inTabCols)-1):
        del inTab[inTabCols[x]]
    inTabshRNA = list(inTab["shRNA"])
    inTabGene = [x.split(".")[0] for x in inTabshRNA]
    inTab["geneName"] = inTabGene
    inTab_byGeneName = inTab.group_by("geneName")
    
    with open(outFile, "w") as fout:
        wfout = csv.writer(fout, delimiter=",")
        wfout.writerow(["GeneName", inTab.colnames[1]])
        for groupX in inTab_byGeneName.groups:
            nameX = groupX["geneName"][0]
            groupXAvg = sum(list(groupX.columns[1]))/len(list(groupX.columns[1]))
            newRow = [nameX, groupXAvg]
            wfout.writerow(newRow)

os.chdir("/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190513_Exp35Exp56_nbPctl-All/3_gate_comparisons")
for file in glob.glob("*.csv"):
    avgByGene(file)

###--- Integrate data
def mergeTables(fileList):
    #fileList = ["P15-21_Q4minusQ1_byGene.csv"]
    outFileName = fileList[0].split("_")[-2] + "_" + fileList[0].split("_")[-1]
    with open(outFileName, "w") as fout:
        wfout = csv.writer(fout, delimiter=",")
        wfout.writerow(["Pool", "Gene", "nbPctgShift"])
        for file in fileList:
            with open(file, "r") as fin:
                pool = file.split("_")[0]
                rfin = csv.reader(fin, delimiter=",")
                next(rfin)
                for row in rfin:
                    wfout.writerow([pool] + row)
                    
mergeTables(["P1-7_InputMinusAvg_byGene.csv",
            "P8-14_InputMinusAvg_byGene.csv",
            "P15-21_InputMinusAvg_byGene.csv"])

mergeTables(["P1-7_Q3minusOther_byGene.csv",
            "P8-14_Q3minusOther_byGene.csv",
            "P15-21_Q3minusOther_byGene.csv"])

mergeTables(["P1-7_Q4minusQ1_byGene.csv",
            "P8-14_Q4minusQ1_byGene.csv",
            "P15-21_Q4minusQ1_byGene.csv"])




