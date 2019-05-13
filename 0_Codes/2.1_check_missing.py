#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 12 20:30:31 2019

@author: yolandatiao
"""


########## Check missing shRNAs in each pool ##########
# Author: Huitian (Yolanda) Diao
# May 12th, 2019

########## Import ##########
import os
import csv
import glob
from operator import itemgetter
from astropy.io import ascii


########## Main ##########
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVitro/DataStats"
os.chdir(wk_dir)

def check_missing(inDir, refFile, outName):
    #inDir = "/Volumes/Yolanda/CRF_Screen/InVitro/1_2_normtoall_ZP"
    #refFile = "/Volumes/Yolanda/CRF_Screen/InVitro/position_Ref.csv"
    refTab = ascii.read(refFile)
    refPos = list(refTab.columns[0])
    for file in glob.glob("%s/*.csv"%inDir):
        file_name = file.split("/")[-1]
        fileTab = ascii.read(file)
        filePos = fileTab.columns[0]
        construct_exist = []
        for Pos in refPos:
            if Pos in filePos:
                construct_exist.append("Y")
            else:
                construct_exist.append("Missing")
        refTab[file_name] = construct_exist
    ascii.write(refTab, outName, format="csv", overwrite=True)


### Check output files...
in_dir = "/Volumes/Yolanda/CRF_Screen/InVitro/1_2_normtoall_ZP"
ref_file = "/Volumes/Yolanda/CRF_Screen/InVitro/position_Ref.csv"
out_name = "1_2_normtoall_ZP_construct_existance_summary.csv"

check_missing(in_dir, ref_file, out_name)


in_dir = "/Volumes/Yolanda/CRF_Screen/InVitro/1_1_shRNAmatched"
ref_file = "/Volumes/Yolanda/CRF_Screen/InVitro/position_Ref.csv"
out_name = "1_1_shRNAmatched_construct_existance_summary.csv"

check_missing(in_dir, ref_file, out_name)



######################################## Find missing value from Meghan's old data... ########################################
##### So much for lossing original data by deleting info from source...
##### Just lost 2 days of my life try to recover some 5 year old mistake...
##### That's why you don't freaking mess with your source data...
'''
def findMissing(inFile, refFile):
    # Create reference
    refTab = ascii.read(refFile)
    refPlate = list(refTab.columns[0])
    refPlate = [x.replace(" ", "") for x in refPlate]
    refWell = list(refTab.columns[2])
    refPos = ["%s_%s"%(x, y) for index, (x, y) in enumerate(zip(refPlate, refWell))]
    refVal = list(refTab.columns[3])
    refConstruct = list(refTab.columns[1])
    
    outFile = inFile.replace(".csv", "whole.csv")
    with open(inFile, "r") as fin:
        with open(outFile, "w") as fout:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            header = next(rfin)
            wfout.writerow(header)
            for row in rfin:
                if row[3] == "":
                    rowPos = row[0]
                    if rowPos in refPos:
                        rowPosIdx = refPos.index(rowPos)
                        if row[2] == refConstruct[rowPosIdx]:
                            row[3] = refVal[rowPosIdx]
                        else:
                            print("Error! shRNA not match with reference:")
                            print(row)
                    else:
                        print ("Error! Position not found!")
                wfout.writerow(row)


wk_dir = "/Volumes/Yolanda/CRF_Screen/InVitro/1_1_shRNAmatched"
ref_dir = "/Volumes/Yolanda/CRF_Screen/InVitro/Megan_originaldata/Meghan_sep"
os.chdir(wk_dir)

for file in glob.glob("Amt_*"):
    ref_file = file.replace("GeoMean", "_GeoMean").replace("Percentage", "_Percentage").replace("_shRNA", "")
    ref_file = "%s/%s"%(ref_dir, ref_file)
    print(file)
    print("--------------------")
    if os.path.isfile(ref_file):
        findMissing(file, ref_file)
    else:
        print("%s not found!" %ref_file)

for file in glob.glob("Amt_*shRNA.csv"):
    ref_file = file.replace("GeoMean", "_GeoMean").replace("Percentage", "_Percentage").replace("_shRNA", "")
    ref_file = "%s/%s"%(ref_dir, ref_file)
    if os.path.isfile(ref_file) == False:
        print(file)
'''

# Not found:
'''
Amt_CD103GeoMean_100U_shRNA.csv
Amt_CD103GeoMean_10U_shRNA.csv
Amt_CD127GeoMean_10U_shRNA.csv
Amt_CXCR3Percentage_100U_shRNA.csv
Amt_CXCR3Percentage_10U_shRNA.csv
Amt_Tim3GeoMean_100U_shRNA.csv
Amt_Tim3GeoMean_10U_shRNA.csv
'''

################################################################################
################################################################################
################################################################################
# Case closed.  Plate 14 panel 1 Well B02 - B09 were missing from the start... #
################################################################################
################################################################################
################################################################################


#####----- For fixing plate 14






