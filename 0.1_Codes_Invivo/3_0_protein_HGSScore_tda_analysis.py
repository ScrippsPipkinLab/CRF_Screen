# -*- coding: utf-8 -*-

############################## Protein TDA network analysis ##############################
# Author: Huitian (Yolanda) Diao
# June 7th, 2019

##### References #####
# HGSScore: Marcon, Edyta, et al. "Human-chromatin-related protein interactions identify a demethylase complex required for chromosome segregation." Cell reports 8.1 (2014): 297-310.
# UMAP: https://umap-learn.readthedocs.io/en/latest/how_umap_works.html


############################## Import ##############################
import umap
import pandas as pd
import os
import matplotlib.pyplot as plt
from astropy.io import ascii
import numpy as np
#import seaborn as sns

os.chdir("/Volumes/Yolanda/CRF_Screen/0.1_Codes_Invivo")
import umap_umap

wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/3_CRM_protein_HGSScore"
os.chdir(wk_dir)


####################  Test conditions ####################
dist_file = "HGSCore_only-CRF_sq_minus_dist.csv"
dist_df = pd.read_csv(dist_file, sep=",", header=0, index_col=0)
dist_df.index.values # al row names

###----- Color selection
complex_ref = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/2_GO_terms/CRM_complexes_count.csv"
complex_tab = ascii.read(complex_ref)
complex_tab_genes = list(complex_tab["gene_name"])
swisnf_idx = list(complex_tab["SWI/SNF complex"])
nurd_idx = list(complex_tab["NuRD complex"])
mll_idx = list(complex_tab["MLL1 complex"])
swisnf_gn = [x for index, (x, y) in enumerate(zip(complex_tab_genes, swisnf_idx)) if y == "Yes"]
nurd_gn = [x for index, (x, y) in enumerate(zip(complex_tab_genes, nurd_idx)) if y == "Yes"]
mll_gn = [x for index, (x, y) in enumerate(zip(complex_tab_genes, mll_idx)) if y == "Yes"]
swisnf_gn = [x.upper() for x in swisnf_gn]
nurd_gn = [x.upper() for x in nurd_gn]
mll_gn = [x.upper() for x in mll_gn]
color_names = []
for i in dist_df.index.values:
    if i.split("|")[1] in swisnf_gn:
        color_names.append("tomato")
    elif i.split("|")[1] in nurd_gn:
        color_names.append("deepskyblue")
    elif i.split("|")[1] in mll_gn:
        color_names.append("limegreen")
    else:
        color_names.append("grey")
color_names

###----- Plot
z = 0
fig = plt.figure()
fig.set_size_inches(45, 60)
for index_i, nn in enumerate([2,4,5,6,8,10,12,15,20]):
    for index_j, mdist in enumerate([0, 0.1, 0.3, 0.4, 0.5, 0.75, 1]):
        z += 1
        embedding = umap.UMAP(n_neighbors=nn, min_dist=mdist, metric='precomputed', 
                              random_state=1234, n_epochs=500,).fit_transform(dist_df.values)
        ax = plt.subplot(9,7, z)
        ax.set_title("n_neighbors: %s, min_dist: %s"%(nn, mdist))
        ax.scatter(embedding[:, 0], embedding[:, 1], c=color_names)
plt.savefig("N_epoches-500_random_state-1234.png")



####################  Find edges ####################
def umap_and_edges(nn, mdist):
    nn = 10
    mdist = 0.4
    dist_file = "HGSCore_only-CRF_sq_dist.csv"
    dist_df = pd.read_csv(dist_file, sep=",", header=0, index_col=0)
    gns = dist_df.index.values # al row names
    gns = [x.split("|")[1] for x in gns]
    umap_embedding = umap.UMAP(n_neighbors=nn, min_dist=mdist, metric='precomputed', 
                                  random_state=1234, n_epochs=500).fit(dist_df.values)
    
    # Save embedding coordinates
    out_name1 = dist_file.replace(".csv", "_nn%s-mdist%s_coor.csv"%(nn, mdist))
    coor_df = pd.DataFrame(umap_embedding.embedding_)
    coor_df.columns = ["x", "y"]
    coor_df.index = gns
    coor_df = coor_df.rename_axis("gene_name")
    coor_df.to_csv(out_name1)
    
    # Save nearest neighbor information
    out_name2 = dist_file.replace(".csv", "_nn%s-mdist%s_nn.csv"%(nn, mdist))
    with open(out_name2, 'w') as fout:
        wfout = csv.writer(fout, delimiter = ",")
        header = ["idx1", "idx2", "val"]     
        wfout.writerow(header)
        eb_graph = pd.DataFrame(umap_embedding.graph_.todense())
        for i in range(0, eb_graph.shape[1]):
            for j in range(0, eb_graph.shape[1]):
                if eb_graph[i][j] > 0:
                    wfout.writerow([i,j, eb_graph[i][j]])
                    
    out_name3 = dist_file.replace(".csv", "_nn%s-mdist%s_nn_dupr.csv"%(nn, mdist))
    nn_dr = ascii.read(out_name2)
    all_idx = ["%s|%s"%(min([x,y]), max(x,y)) for index, (x, y) in enumerate(zip(list(nn_dr["idx1"]), list(nn_dr["idx2"])))]
    nn_dr["all_idx"] = all_idx
    nn_dr_g = nn_dr.group_by("all_idx")
    keyn = len(nn_dr_g.groups.keys)
    
    with open(out_name3, "w") as fout:
        wfout = csv.writer(fout, delimiter=",")
        header = ["idx1", "idx2", "val"]     
        wfout.writerow(header)
        for x in range(0, keyn):
            tabx = nn_dr_g.groups[x]
            if (len(tabx) > 2):
                print("Error!! Dup > 2!! %s"%tabx["all_idx"][x])
            elif (len(tabx) < 2):
                print("Error!! Dup < 2!! %s"%tabx["all_idx"][x])
            else:
                if tabx["val"][0] == tabx["val"][1]:
                    wfout.writerow(list(tabx[0])[:3])
                else:
                    print("Error!! Duplicate different value %s"%tabx["all_idx"][x])
                
    

######################################## Preprocess distance data ######################################## 
from astropy.io import ascii
import csv
import os
import math

os.chdir("/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/3_CRM_protein_HGSScore")

def filterCRF(inFile):
    inFile = "/Users/yolandatiao/Documents/z_Packages/CRM_protein/HGSCore.csv"
    outFile = inFile.replace(".csv", "_only-CRF.csv")
    refFile = "/Volumes/Yolanda/RNAseq_Compilation/CRM_analysis_2/Info/CRM.csv"
    refTab = ascii.read(refFile)
    all_gn = list(refTab.columns[0]) + list(refTab.columns[1])
    all_gn = [x for x in all_gn if x != "NA" and x != "--"]
    all_gn_up = [x.upper() for x in all_gn]
    
    with open(inFile, "r") as fin:
        with open(outFile, "w") as fout:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            wfout.writerow(next(rfin))
            for row in rfin:
                if (row[0].split("|")[1].upper() in all_gn_up) and (row[1].split("|")[1].upper() in all_gn_up):
                    wfout.writerow(row)
    

def bpTOsquare(inFile):
    inFile = "/Users/yolandatiao/Documents/z_Packages/CRM_protein/HGSCore_only-CRF.csv"
    outFile = inFile.replace(".csv", "_sq.csv")
    inTab = ascii.read(inFile)
    p1List = list(inTab.columns[0])
    p2List = list(inTab.columns[1])
    p1Set = list(set(p1List))
    p2Set = list(set(p2List))
    allSet = list(set(p1Set + p2Set))
    p12List = ["%s||%s"%(x,y) for index, (x, y) in enumerate(zip(p1List, p2List))]
    
    with open(outFile, "w") as fout:
        wfout = csv.writer(fout)
        header = ["x"] + allSet
        wfout.writerow(header)
        for i in allSet:
            newrow = [i]
            for j in allSet:
                score = []
                if i == j:
                    score.append(0)
                elif "%s||%s"%(i, j) in p12List:
                    idx = p12List.index("%s||%s"%(i, j))
                    score.append(inTab[idx][2])
                elif "%s||%s"%(j, i) in p12List:
                    idx = p12List.index("%s||%s"%(j, i))
                    score.append(inTab[idx][2])
                if len(score) == 0:
                    newrow.append(0)
                elif len(score) > 1:
                    print("Duplicate found!! %s : %s"%(i,j) )
                    newrow.append(sum(score)/len(score))
                elif len(score) == 1:
                    newrow.append(score[0])
            wfout.writerow(newrow)

def convertDistance_sqDist(inFile):
    inFile = "HGSCore_only-CRF_sq.csv"
    outFile = inFile.replace(".csv", "_dist.csv")
    with open(inFile, "r") as fin:
        with open(outFile, "w") as fout:
            rfin = csv.reader(fin, delimiter = ",")
            wfout = csv.writer(fout, delimiter = ",")
            wfout.writerow(next(rfin))
            for row in rfin:
                newrow = [row[0]]
                for x in range(1, len(row)):
                    newrow.append(100/math.sqrt((float(row[x])+1)))
                wfout.writerow(newrow)

def convertDistance_Minus(inFile):
    inFile = "HGSCore_only-CRF_sq.csv"
    outFile = inFile.replace(".csv", "_minus_dist.csv")
    with open(inFile, "r") as fin:
        with open(outFile, "w") as fout:
            rfin = csv.reader(fin, delimiter = ",")
            wfout = csv.writer(fout, delimiter = ",")
            wfout.writerow(next(rfin))
            min_n = 500
            for row in rfin:
                newrow = [row[0]]
                for x in range(1, len(row)):
                    newrow.append(430 - float(row[x]))
                    if (500 - float(row[x])) < min_n:
                        min_n =  (430 - float(row[x]))
                wfout.writerow(newrow)
    print(min_n)




