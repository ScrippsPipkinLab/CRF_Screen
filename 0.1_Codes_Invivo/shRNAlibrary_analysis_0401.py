#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 11:22:56 2019

@author: yolandatiao
"""

########## shRNA library analysis ##########
# Author: Huitian (Yolanda) Diao
# Apr. 1st, 2019

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
    #inFile = "/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_combined/P1-7_Input_combined.csv"
    #print(inFile)
    inFileName = inFile.split("/")[-1].replace(".csv", "")
    outFileName = inFileName + "_pctg.csv"
    inTab = ascii.read(inFile)
    countList = inTab.columns[2]
    countSum = sum(countList)
    pctgList = [float(x)/countSum*100 for x in countList]
    inTab['pctg'] = pctgList
    ascii.write(inTab, outFileName, format="csv", overwrite="True")

def nbPctg(inFile):
    #inFile = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190401/P1-7_Input_combined_pctg.csv"
    print(inFile)
    inFileName = inFile.split("/")[-1].replace(".csv", "")
    outFileName = inFileName + "_nbPctg.csv"
    inTab = ascii.read(inFile)
    pctgList = inTab.columns[3]
    pctgList_10000 = [x*10000 for x in pctgList]
    allNb = fit_nbinom(np.array(pctgList_10000))
    nbPctgList = [st.nbinom.cdf(x, allNb['size'], allNb['prob'])*100 for x in pctgList_10000]    
    inTab['nbPctg'] = nbPctgList
    ascii.write(inTab, outFileName, format="csv", overwrite="True")  
    
    

########## Main ##########
#----- Directory
ref_file = "/Volumes/Yolanda/CRF_Screen/InVivo/sample_reference.csv"
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190401"
os.chdir(wk_dir)

#----- Calculate percentage of total
for file in glob.glob("/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_combined/*.csv"):
    pctgTotal(file)

#----- Calculate negative binomial distribution percentiles
for file in glob.glob("%s/*.csv"%wk_dir):
    nbPctg(file)


###--- Gate comparisons
#----- Q1vQ4
def Q4minusQ1(q4File, q1File):
    q4File = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190401/P1-7_Q4_combined_pctg_nbPctg.csv"
    q1File = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190401/P1-7_Q1_combined_pctg_nbPctg.csv"
    group = q4File.split("/")[-1].split("_")[0]
    outName = group + "_Q4minusQ1.csv"
    q4Tab = ascii.read(q4File)
    q1Tab = ascii.read(q1File)
    del q4Tab['count']
    del q4Tab['pctg']
    del q4Tab['type']
    del q1Tab['count']
    del q1Tab['pctg']
    del q1Tab['type']
    q4Tab["nbPctg"].name = "nbPctg_Q4"
    q1Tab["nbPctg"].name = "nbPctg_Q1"
    allTab = join(q4Tab, q1Tab, join_type="inner", keys="shRNA")
    q4minusq1 = [x-y for index, (x,y) in enumerate(zip(list(allTab["nbPctg_Q4"]), list(allTab["nbPctg_Q1"])))]
    allTab["q4minusq1_nbPctg"] = q4minusq1
    ascii.write(allTab, outName, format="csv", overwrite=True)
    





















