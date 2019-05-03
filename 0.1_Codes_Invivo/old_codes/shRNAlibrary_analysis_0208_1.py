#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 13:25:19 2019

@author: yolandatiao
"""

########## shRNA library analysis ##########
# Author: Huitian (Yolanda) Diao
# Jan 23nd, 2019

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
import scipy.stats as st
import numpy as np

from scipy.special import gammaln
from scipy.special import psi
from scipy.misc import factorial
from scipy.optimize import fmin_l_bfgs_b as optim
import sys

from sklearn.preprocessing import quantile_transform


########## Self-defined functions ##########
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

def NormZscore(valuex, listx):
    listx_std = np.std(listx)
    listx_mean = float(sum(listx))/len(listx)
    return((valuex - listx_mean)/listx_std)

def NormPercentile(valuex, listx):
    zScore = NormZscore(valuex, listx)
    return(st.norm.cdf(zScore)*100)

def pctgCtrl(tabX):
    #tabX = ascii.read("/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_combined/P1-7_Input_combined.csv")
    tabX_types = list(tabX['type'])
    tabX_count = list(tabX['count'])
    ctrl_idx = [idx for idx, x in enumerate(tabX_types) if x == "control"]
    tabCtrl_count = [tabX_count[x] for x in ctrl_idx]
    percentile = [NormPercentile(x, tabCtrl_count) for x in tabX_count]
    tabX['NormCtrlPercentile'] = percentile
    return(tabX)

def sltCtrl(tabX):
    #tabX = ascii.read("/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_combined/P1-7_Input_combined.csv")
    tabX_shRNA = list(tabX["shRNA"])
    tabX_gene = [x.split(".")[0].lower() for x in tabX_shRNA]
    slt_ctrl_idx = [idx for idx, x in enumerate(tabX_gene) if ((x=='cd19') or x=='cd4')]
    slt_ctrl_tab = tabX[slt_ctrl_idx]
    ctrl_avg = float(sum(list(slt_ctrl_tab['count'])))/len(slt_ctrl_tab)
    return([ctrl_avg, list(slt_ctrl_tab['count']), list(slt_ctrl_tab['shRNA'])])

def sltCtrlTabs(fileList):
    #fileList = ['/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_combined/P1-7_Input_combined.csv',
    #            '/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_combined/P1-7_Q1_combined.csv',
    #            '/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_combined/P1-7_Q2_combined.csv',
    #            '/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_combined/P1-7_Q3_combined.csv',
    #            '/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_combined/P1-7_Q4_combined.csv']
    
    fileNameNoFormatList = [x.split("/")[-1].replace(".csv","") for x in fileList]
    tabList = []
    for file in fileList:
        tabList.append(ascii.read(file))
    ctrlTabList = []
    for index,(tabx, name) in enumerate(zip(tabList, fileNameNoFormatList)):
        ctrl = sltCtrl(tabx)
        ctrlTabx = Table()
        ctrlTabx['shRNA'] = ctrl[2]
        ctrlTabx[name] = ctrl[1]
        ctrlTabList.append(ctrlTabx)
    allCtrlTab = ctrlTabList[0]
    for x in range(1, len(ctrlTabList)):
        allCtrlTab = join(allCtrlTab, ctrlTabList[x], join_type = 'outer', keys='shRNA')
    return(allCtrlTab)
 
def normInput(fileList):
    #fileList = ['/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_combined/P1-7_Input_combined.csv',
    #            '/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_combined/P1-7_Q1_combined.csv',
    #            '/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_combined/P1-7_Q2_combined.csv',
    #            '/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_combined/P1-7_Q3_combined.csv',
    #            '/Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/3_combined/P1-7_Q4_combined.csv']
    
    fileNameNoFormatList = [x.split("/")[-1].replace(".csv","") for x in fileList]  
    tabList = []
    for file in fileList:
        tabList.append(ascii.read(file))
    
    # Convert to percentages
    # Quantile transformation
    for idx, tab in enumerate(tabList):
        tab_counts = list(tab['count'])
        tab_counts_sum = sum(tab_counts)
        tab_counts_pctg = [float(x)/tab_counts_sum*100 for x in tab_counts]
        tab_counts_pctg_array = np.array(tab_counts_pctg).reshape(-1,1)
        tab_counts_quantile = quantile_transform(tab_counts_pctg_array, n_quantiles= len(tab_counts_pctg_array), random_state=0)
        tabList[idx]['quantile'] = [j for sub in tab_counts_quantile for j in sub]

    # Z-Score
    newTab=Table()
    newTab['shRNA'] = tabList[0]['shRNA']
    sample_dict = {}
    for x in fileNameNoFormatList:
        sample_dict[x] = []
    for x in range(0, len(newTab)):
        list_x = []
        for tab in tabList:
            list_x.append(tab['quantile'][x])
        zscore_x = st.zscore(np.array(list_x))
        for idx, key in enumerate(sample_dict.keys()):
            sample_dict[key] = sample_dict[key] + [zscore_x[idx]]
    
    for idx, key in enumerate(sample_dict.keys()):
        newTab[key] = sample_dict[key]
    
    newTab[0]
    del newTab[newTab.colnames[1]]
    return(newTab)
    


########## Main ##########
ref_file = "/Volumes/Yolanda/CRF_Screen/InVivo/sample_reference.csv"
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/1_1_Norm/20190208_1"
os.chdir(wk_dir)

#####----- Calculate percentage shift versus control --> summarize data into a table
### Build a reference dictionary with reference table
ref_tab = ascii.read(ref_file)
group_replicate_set = ["%s_rep%s"%(g, r) for index, (g,r) in enumerate(zip(list(ref_tab['screen_group']), list(ref_tab['replicate_number'])))]
group_replicate_set = list(set(group_replicate_set))
group_replicate_set
conditions_set = list(set(list(ref_tab['condition'])))
ref_dict = {}
for gr in group_replicate_set:
    gr_list = gr.replace("_rep",",").split(",")
    tablex = Table()
    tablex['conditions'] = conditions_set
    sample_names = []
    sample_files = []
    input_control_files = []
    for cond in conditions_set:
        found_cond = []
        for x in range(0, len(ref_tab)):
            if (ref_tab['screen_group'][x] == gr_list[0]) and (ref_tab['replicate_number'][x] == int(gr_list[1])) and (ref_tab['condition'][x] == cond):
                found_cond.append(x)
        if len(found_cond) == 1:
            sample_names.append(ref_tab.columns[0][found_cond[0]])
            sample_files.append(ref_tab['sample_file'][found_cond[0]])
            input_control_files.append(ref_tab['input_file'][found_cond[0]])
        else:
            sample_names.append('NA')
            sample_files.append('NA')
            input_control_files.append('NA')
            if len(found_cond) > 1:
                print("Error! Multiple matches found for %s: %s"%(gr, cond))
            else:
                print("Error! No match found for %s: %s"%(gr, cond))
    tablex['sample_names'] = sample_names
    tablex['sample_files'] = sample_files
    tablex['input_control_files'] = input_control_files
    ref_dict[gr] = tablex

### Loop through each of the groups
table_list = []
for key, tab in ref_dict.items():
    #key = 'P8-14_rep1'
    #tab = ref_dict[key]
    
    group_table_list = []
    sample_files = list(tab['sample_files'])
    input_control_file = list(tab['input_control_files'])[0]
    all_files_list = [input_control_file] + sample_files
    newtab = normInput(all_files_list) # Input file normalization
    for x in range(1,len(newtab.colnames)):
        newtab[newtab.colnames[x]].name = tab['conditions'][(x-1)]
    newtab['exp'] = [key for x in range(0, len(newtab))]
    table_list.append(newtab)
whole_tab = vstack(table_list)
ascii.write(whole_tab, "PctgQuantileZscore.csv", format="csv", overwrite=True)



#####----- Statistic tests for significant difference between mean v.s. shRNA
# Comparisons: 
# - shX.x in Q1 v.s. Q1
# - shX.x in Q2 v.s. Q2
# - shX.x in Q3 v.s. Q3
# - shX.x in Q4 v.s. Q4
# - shX.x in Q1 v.s. Q2
# - shX.x in Q1 v.s. Q3
# - shX.x in Q1 v.s. Q4
# - shX.x in Q2 v.s. Q3
# - shX.x in Q2 v.s. Q4
# - shX.x in Q3 v.s. Q4
'''
shRNA_list = list(whole_tab['shRNA_names'])
gene_names = [x.split("_")[0].split(".")[0] for x in shRNA_names]
gene_names_set = list(set(gene_names))
shRNA_names = [x.split("_")[0] for x in shRNA_names]
shRNA_names_set = list(set(shRNA_names))

gene_dict = {}
for genex in gene_names_set:
    genex_list = [idx for idx, name in enumerate(gene_names) if name == genex]
    gene_dict[genex] = genex_list

whole_tab['Q1-Q2'] = [(x-y) for index, (x, y) in enumerate(zip(list(whole_tab['Q1']), list(whole_tab['Q2'])))]
whole_tab['Q1-Q3'] = [(x-y) for index, (x, y) in enumerate(zip(list(whole_tab['Q1']), list(whole_tab['Q3'])))]
whole_tab['Q1-Q4'] = [(x-y) for index, (x, y) in enumerate(zip(list(whole_tab['Q1']), list(whole_tab['Q4'])))]
whole_tab['Q2-Q3'] = [(x-y) for index, (x, y) in enumerate(zip(list(whole_tab['Q2']), list(whole_tab['Q3'])))]
whole_tab['Q2-Q4'] = [(x-y) for index, (x, y) in enumerate(zip(list(whole_tab['Q2']), list(whole_tab['Q4'])))]
whole_tab['Q3-Q4'] = [(x-y) for index, (x, y) in enumerate(zip(list(whole_tab['Q3']), list(whole_tab['Q4'])))]


p_val_file = "byGene_p_val_NormPctgShift.csv"
avg_file = "byGene_avg_NormPctgShift.csv"
with open(p_val_file, "w") as foutp:
    with open(avg_file, "w") as fouta:
        wfoutp = csv.writer(foutp, delimiter=",")
        wfouta = csv.writer(fouta, delimiter=",")
        wfoutp.writerow(['gene_name'] + whole_tab.colnames[1:])
        wfouta.writerow(['gene_name'] + whole_tab.colnames[1:])
        for gene, gene_list in gene_dict.items():
            gene_tab = whole_tab[gene_list]
            p_val_list = []
            avg_list = []
            for x in range(1, len(whole_tab.colnames)):
                p_val = st.ttest_ind(list(gene_tab.columns[x]), list(whole_tab.columns[x])).pvalue
                avg = sum(list(gene_tab.columns[x])) / len(list(gene_tab.columns[x]))
                p_val_list.append(str(p_val))
                avg_list.append(str(avg))
            wfoutp.writerow([gene] + p_val_list)
            wfouta.writerow([gene] + avg_list)
            
'''        




























    