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


in_dir = "/Volumes/Yolanda/CRF_Screen/InVitro/1_2_normtocontrol_ZP/Amt"
ref_file = "/Volumes/Yolanda/CRF_Screen/InVitro/position_Ref.csv"
out_name = "1_2_normtocontrol_ZP_Amt_construct_existance_summary.csv"

check_missing(in_dir, ref_file, out_name)

in_dir = "/Volumes/Yolanda/CRF_Screen/InVitro/1_2_normtocontrol_ZP/AmtGFP"
ref_file = "/Volumes/Yolanda/CRF_Screen/InVitro/position_Ref.csv"
out_name = "1_2_normtocontrol_ZP_AmtGFP_construct_existance_summary.csv"

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
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVitro/1_1_0_shRNAmatched"
os.chdir(wk_dir)

p1_list = ["CD44", "CD25", "CD62L", "CXCR3", "CD103", "CD127"]

o_pos = ["14_B02", "14_B03", "14_B04", "14_B05", "14_B06", "14_B07", "14_B08", "14_B09", "14_B10", "14_B11", 
         "14_C02", "14_C03", "14_C04", "14_C05", "14_C06", "14_C07", "14_C08", "14_C09", "14_C10", "14_C11", 
         "14_D02", "14_D03", "14_D04", "14_D05", "14_D06", "14_D07", "14_D08", "14_D09", "14_D10", "14_D11", 
         "14_E02", "14_E03", "14_E04", "14_E05", "14_E06", "14_E07", "14_E08", "14_E09", "14_E10", "14_E11", 
         "14_F02", "14_F03", "14_F04", "14_F05", "14_F06", "14_F07", "14_F08", "14_F09", "14_F10", "14_F11", 
         "14_G02", "14_G03", "14_G04", "14_G05", "14_G06", "14_G07", "14_G08", "14_G09", "14_G10", "14_G11"]

n_pos = [                                                                      "14_B09", "14_B10", "14_B11", 
         "14_C02", "14_C03", "14_C04", "14_C05", "14_C06", "14_C07", "14_C08", "14_C09", "14_C10", "14_C11", 
         "14_D02", "14_D03", "14_D04", "14_D05", "14_D06", "14_D07", "14_D08", "14_D09", "14_D10", "14_D11", 
         "14_E02", "14_E03", "14_E04", "14_E05", "14_E06", "14_E07", "14_E08", "14_E09", "14_E10", "14_E11", 
         "14_F02", "14_F03", "14_F04", "14_F05", "14_F06", "14_F07", "14_F08", "14_F09", "14_F10", "14_F11", 
         "14_G02", "14_G03", "14_G04", "14_G05", "14_G06", "14_G07", "14_G08", "14_G09", "14_G10", "14_G11",
         "14_B02", "14_B03", "14_B04", "14_B05", "14_B06", "14_B07", "14_B08"]
n_shRNA = ["Prmt7.1", "Prmt2.1", "Rnf20.1", "Ring1.2", "Prmt5.1", "Prmt1.2", "Rnf20.2", "Rbbp7.1", "Rbbp5.2", "Prmt3.2", "Rnf217.1", "Rbbp4.2", "Prmt8.1", "Prmt2.2", "Rnf217.2", "Rbbp4.3", "Prmt6.1", "Prmt5.2", "Rnf217.3", "Rbbp7.2", "Prmt7.2", "Prmt1.3", "Rnf40.2", "Prmt8.2", "Psip1.2", "Prmt1.4", "Ring1.3", "Rbbp5.3", "Prmt6.2", "Prmt2.3", "Rnf20.3", "Rbbp7.3", "Prmt7.3", "Prmt3.3", "Rnf40.3", "Rbbp5.4", "Prmt6.3", "Prmt5.3", "Ring1.4", "Psip1.3", "Prmt8.3", "Prmt5.4", "Rnf20.4", "Rbbp7.4", "Prmt8.4", "Prmt2.4", "Rnf40.4", "Rnf217.4", "Prmt6.4", "Prmt3.4", "Rbbp4.4", "Prmt7.4", "Psip1.4", "Prmt1.1", "Rbbp4.1", "Rbbp5.1", "Psip1.1", "Prmt3.1", "Rnf40.1", "Ring1.1"]

for file in glob.glob("*.csv"):
    new_name = file.replace(".csv", "_fixed.csv")
    marker = file.split("_")[1].replace("Percentage", "").replace("GeoMean", "")
    fileTab = ascii.read(file)
    if ("10U" in file) and (marker in p1_list):
        fileTab[1][0]
        for x in range(0, len(fileTab)):
            if str(fileTab[x][0].split("_")[0]) == '14':
                pos_x = fileTab[x][0]
                pos_x_idx = o_pos.index(pos_x)
                fileTab[x][0] = n_pos[pos_x_idx]
                fileTab[x][2] = n_shRNA[pos_x_idx]
    ascii.write(fileTab, new_name, format="csv", overwrite=True)
            





