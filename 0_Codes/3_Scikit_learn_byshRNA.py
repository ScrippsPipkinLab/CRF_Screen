#!/usr/bin/env python3
# -*- coding: utf-8 -*-
########## In vitro screen CRF - Scikit_learn ##########
# Author: Huitian (Yolanda) Diao
# Dec 3, 2018

########## Import ##########
import os
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import decomposition
from sklearn import datasets
import pandas as pd
from sklearn.decomposition import PCA
from astropy.io import ascii # For using ascii table to open csv
from astropy.table import Table, Column, join   # For using astropy table functions
from sklearn import manifold, datasets
import matplotlib.patches as mpatches
import csv

########## Main ##########
wkdir = "/Users/yolandatiao/Desktop/CRF_Screen/1_2_normtocontrol_ZP/Compiled_Amt"
os.chdir(wkdir)

#--- Read tpm file
zscore_file = "Amt_normbycontrolZP_z-score.csv"
zscore_df = ascii.read(zscore_file)
del zscore_df["Position"]
zscore_names = zscore_df.colnames
shRNA_names = list(zscore_df["shRNA"])

#--- Convert to numpy array
zscore_arr = np.array([list(zscore_df[0])[1:]])
zscore_arr.shape
for x in range(1, len(zscore_df)):
    #print(len(list(tpm_df.columns[x])))
    zscore_arr = np.append(zscore_arr, [list(zscore_df[x])[1:]], axis = 0)
zscore_arr.shape

#--- PCA analysis
pca = PCA(n_components=15)
zscore_arr_new = pca.fit_transform(zscore_arr)
zscore_arr_new.shape # Check new array shape
pca.components_.shape # Principal axes in feature space
pca.explained_variance_ratio_
np.sum(pca.explained_variance_ratio_) # Percentage explained by PCs
plt.scatter(zscore_arr_new[:, 0], zscore_arr_new[:, 1], c="red")

#--- Tsne
def tsne_conversion(zscore_arr, per_n): 
    n_components = 2
    tsne = manifold.TSNE(n_components=n_components, 
                         perplexity = per_n, init='pca', random_state=0)
    tpm_arr_tsne = tsne.fit_transform(zscore_arr)
    tpm_arr_tsne.shape
    plt.scatter(tpm_arr_tsne[:, 0], tpm_arr_tsne[:, 1], c="red")
    
    out_name = "AAmt_normbycontrolZP_z-score_tsne_per%s.csv"%per_n
    with open(out_name, "w") as fout:
        wfout = csv.writer(fout, delimiter = ",")
        wfout.writerow(["gene_name", "t1", "t2"])
        for x in range(0,len(tpm_arr_tsne)):
            wfout.writerow([shRNA_names[x]] + list(tpm_arr_tsne[x]))
       
tsne_conversion(zscore_arr, 5)
tsne_conversion(zscore_arr, 10)
tsne_conversion(zscore_arr, 15)