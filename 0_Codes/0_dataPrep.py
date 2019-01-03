#!/usr/bin/env python3
# -*- coding: utf-8 -*-

########## Data Preparation for CRF library in vitro screening results ##########
# Author: Huitian (Yolanda) Diao
# Dec 1st, 2018


########## Import ##########
import os
import csv
import glob
from operator import itemgetter
from astropy.io import ascii

########## Self-defined functions ##########
def find_AG(in_dir, header_idx, header_names, factors, fps):
    os.chdir(in_dir)
    file_list = []
    
    for file in glob.glob("*.csv"):
        if ("Amet" not in file) and ("GFP" not in file):
            file_list.append(file)
            
    for file in file_list:
        #os.chdir("/Users/yolandatiao/Desktop/CRF_Screen/Meghan_originaldata/Screen_markers/Screen1_panel1_10U")
        #file = "1 10U p1.csv"
        #header_names = ["Well", "shRNA"]
        #header_idx = [2, 1]
        out_name = file.replace(" ", "_").replace(".csv", "_AmetGFP.csv")
        with open(file, "r") as fin:
            print(" ")
            print(file)
            with open(out_name, "w") as fout:
                rfin = csv.reader(fin, delimiter = ",")
                wfout = csv.writer(fout, delimiter = ",")
                write = False
                write_index = header_idx.copy()
                write_header = header_names.copy()
                #print(header_names)
                #print(write_header)
                for row in rfin:
                    if write == True:
                        wfout.writerow(itemgetter(*write_index)(row))
                    else:
                        for i in range(0, len(row)):
                            if ("Amet+" in row[i]) and ("GFP+" in row[i]):
                                if ("Mean" in row[i]) or ("%" in row[i]):
                                    colname = row[i]
                                    has_key = False
                                    for factor in factors:
                                        if factor in colname:
                                            has_key = True
                                    if has_key == False:
                                        print("Missing Label: %s"%colname)
                                        for x in range(0, len(fps)):
                                            if fps[x] in colname:
                                                colname += "_%s"%factors[x]
                                    write = True
                                    write_header.append(colname)
                                    write_index.append(i)
                        if write == True:
                            print(write_header)
                            wfout.writerow(write_header)

def find_A(in_dir, header_idx, header_names, factors, fps):
    os.chdir(in_dir)
    file_list = []
    
    for file in glob.glob("*.csv"):
        if ("Amet" not in file) and ("GFP" not in file):
            file_list.append(file)
    for file in file_list:
        #os.chdir("/Users/yolandatiao/Desktop/CRF_Screen/Meghan_originaldata/Screen_markers/Screen1_panel1_10U")
        #file = "1 10U p1.csv"
        #header_names = ["Well", "shRNA"]
        #header_idx = [2, 1]
        out_name = file.replace(" ", "_").replace(".csv", "_Amet.csv")
        with open(file, "r") as fin:
            print(" ")
            print(file)
            with open(out_name, "w") as fout:
                rfin = csv.reader(fin, delimiter = ",")
                wfout = csv.writer(fout, delimiter = ",")
                write = False
                write_index = header_idx.copy()
                write_header = header_names.copy()
                #print(header_names)
                #print(write_header)
                for row in rfin:
                    if write == True:
                        wfout.writerow(itemgetter(*write_index)(row))
                    else:
                        for i in range(0, len(row)):
                            if ("Amet+" in row[i]) and ("GFP+" not in row[i]):
                                if ("Mean" in row[i]) or ("%" in row[i]):
                                    colname = row[i]
                                    has_key = False
                                    for factor in factors:
                                        if factor in colname:
                                            has_key = True
                                    if has_key == False:
                                        print("Missing Label: %s"%colname)
                                        for x in range(0, len(fps)):
                                            if fps[x] in colname:
                                                colname += "_%s"%factors[x]
                                    write = True
                                    write_header.append(colname)
                                    write_index.append(i)
                        if write == True:
                            print(write_header)
                            wfout.writerow(write_header)

def compileFactor(in_dir, factor_name, unit_n):
    #in_dir = "/Users/yolandatiao/Desktop/CRF_Screen/Meghan_originaldata/Screen_markers"
    #factor_name = "CXCR3"
    #unit_n = "10U"
    os.chdir(in_dir)
    file_10U = []
    print(factor_name)
    for filename in glob.iglob("*%s/*AmetGFP.csv"%unit_n, recursive=True):
        file_10U.append(filename)
        print(filename)
    print("")
    out_10U_Geo = "AmtGFP_" + factor_name + "GeoMean_%s.csv"%unit_n
    out_10U_P = "AmtGFP_" + factor_name + "Percentage_%s.csv"%unit_n
    with open(out_10U_Geo, "w") as foutg:
        with open(out_10U_P, "w") as foutp:
            wfoutg = csv.writer(foutg, delimiter = ",")
            wfoutp = csv.writer(foutp, delimiter = ",")
            wfoutg.writerow(["Plate", "Well", "GeoMean"])
            wfoutp.writerow(["Plate", "Well", "Percentage"])

            for file in file_10U:
                g_idx = None
                p_idx = None
                plate_n = file.split("/")[-1].split("_")[0]
                if plate_n == "control":
                    plate_n += file.split("_")[0].replace("Screen", "")
                well_idx = None
                with open(file, "r") as fin:
                    rfin = csv.reader(fin, delimiter = ",")
                    header = next(rfin)
                    for i in range(0,len(header)):
                        if factor_name in header[i]:
                            if "GeoMean" in header[i]:
                                g_idx = i
                            if "%" in header[i]:
                                p_idx = i
                        if header[i] == "Well":
                            well_idx = i
                    for row in rfin:
                        newrowg = []
                        newrowg.append(plate_n)
                        newrowp = []
                        newrowp.append(plate_n)
                        if g_idx != None:
                            newrowg.append(row[well_idx])
                            newrowg.append(row[g_idx])
                            wfoutg.writerow(newrowg)
                        if p_idx != None:
                            newrowp.append(row[well_idx])
                            newrowp.append(row[p_idx])
                            wfoutp.writerow(newrowp)

def compileFactor_Amt(in_dir, factor_name, unit_n):
    #in_dir = "/Users/yolandatiao/Desktop/CRF_Screen/Meghan_originaldata/Screen_markers"
    #factor_name = "CXCR3"
    #unit_n = "10U"
    os.chdir(in_dir)
    file_10U = []
    print(factor_name)
    for filename in glob.iglob("*%s/*Amet.csv"%unit_n, recursive=True):
        file_10U.append(filename)
        print(filename)
    print("")
    out_10U_Geo = "Amt_" + factor_name + "GeoMean_%s.csv"%unit_n
    out_10U_P = "Amt_" + factor_name + "Percentage_%s.csv"%unit_n
    with open(out_10U_Geo, "w") as foutg:
        with open(out_10U_P, "w") as foutp:
            wfoutg = csv.writer(foutg, delimiter = ",")
            wfoutp = csv.writer(foutp, delimiter = ",")
            wfoutg.writerow(["Plate", "Well", "GeoMean"])
            wfoutp.writerow(["Plate", "Well", "Percentage"])

            for file in file_10U:
                g_idx = None
                p_idx = None
                plate_n = file.split("/")[-1].split("_")[0]
                if plate_n == "control":
                    plate_n += file.split("_")[0].replace("Screen", "")
                well_idx = None
                with open(file, "r") as fin:
                    rfin = csv.reader(fin, delimiter = ",")
                    header = next(rfin)
                    for i in range(0,len(header)):
                        if factor_name in header[i]:
                            if "GeoMean" in header[i]:
                                g_idx = i
                            if "%" in header[i]:
                                p_idx = i
                        if header[i] == "Well":
                            well_idx = i
                    for row in rfin:
                        newrowg = []
                        newrowg.append(plate_n)
                        newrowp = []
                        newrowp.append(plate_n)
                        if g_idx != None:
                            newrowg.append(row[well_idx])
                            newrowg.append(row[g_idx])
                            wfoutg.writerow(newrowg)
                        if p_idx != None:
                            newrowp.append(row[well_idx])
                            newrowp.append(row[p_idx])
                            wfoutp.writerow(newrowp)                           
                                
def std_Amt(infile):
    outfile = infile.replace("Amt", "standard_Amt")
    with open(infile, "r") as fin:
        with open(outfile, "w") as fout:
            rfin = csv.reader(fin, delimiter = ",")
            wfout = csv.writer(fout, delimiter = ",")
            header = next(rfin)
            wfout.writerow(["plate_well", "shRNA", header[3]])
            for row in rfin:
                new_row = []
                row_plate = row[0]
                row_well = row[2]
                if row_plate != "" and row_well != "":
                    print(row_plate)
                    print(infile)
                    print(row)
                    row_plate_n = ''.join(filter(lambda x: x.isdigit(), row_plate))
                    if row_plate_n != "":
                        row_plate_n = int(row_plate_n)
                    if "control" in row_plate.lower():
                        new_row.append("control%s_%s"%(str(row_plate_n), row_well))
                    else:
                        new_row.append("%s_%s"%(str(row_plate_n), row_well))
                    new_row.append(row[1])
                    new_row.append(row[3])
                    wfout.writerow(new_row)
                    
def shRNAmatch(infile):
    outfile = infile.replace(".csv", "_shRNA.csv")
    ref_file = "/Users/yolandatiao/Desktop/CRF_Screen/position_Ref.csv"
    ref_tab = ascii.read(ref_file)
    pw_list = list(ref_tab["plate_well"])
    shRNA_list = list(ref_tab["shRNA"])
    with open(infile, "r") as fin:
        with open(outfile, "w") as fout:
            rfin = csv.reader(fin, delimiter = ",")
            wfout = csv.writer(fout, delimiter = ",")
            header = next(rfin)
            newheader = ["Position", "Pool", "shRNA"] + [header[2]]
            wfout.writerow(newheader)
            for row in rfin:
                if "control" not in row[0]:
                    plate_n = int(row[0])
                    if plate_n < 8:
                        pool_n = 1
                    elif plate_n <15:
                        pool_n = 2
                    else:
                        pool_n = 3
                else:
                    plate_n = row[0]
                    pool_n = int(''.join(filter(lambda x: x.isdigit(), plate_n)))
                well_n = row[1]
                pw = str(plate_n) + "_" + well_n
                if pw in pw_list:
                    idx = pw_list.index(pw)
                    shRNA = shRNA_list[idx]
                    newrow = [pw, pool_n, shRNA, row[2]]
                    wfout.writerow(newrow)
                                
########## Main ##########
root_dir = "/Users/yolandatiao/Desktop/CRF_Screen/Meghan_originaldata/Screen_markers"

###----- Find all the Amtrine+ GFP+ data
p1_factors = ["CD44", "CD25", "CD62L", "CXCR3", "CD103", "CD127"]
p2_factors = ["PD1", "Lag3", "Tim3"]
p1_fps = ["Alexa 700", "PerCP-Cy5-5", "BV605", "PE-Cy7", "APC", "PE"]
p2_fps = ["BV605", "PerCP-Cy5-5", "PE"]
'''
#--- Screen3_panel1_100U_CXCR3
dir_sp = root_dir + "/Screen3_panel1_100U_CXCR3/Screen3_panel1_100U"
header_sp = ["Well"]
header_idx_sp = [0]
find_AG(dir_sp, header_idx_sp, header_sp, p1_factors, p1_fps)
#--- Screen1_panel1_10U
dir_sp = root_dir + "/Screen1_panel1_10U"
header_sp = ["Well", "shRNA"]
header_idx_sp = [2, 1]
find_AG(dir_sp, header_idx_sp, header_sp, p1_factors, p1_fps)
#--- Screen1_panel1_100U
dir_sp = root_dir + "/Screen1_panel1_100U"
header_sp = ["Well", "shRNA"]
header_idx_sp = [2, 1]
find_AG(dir_sp, header_idx_sp, header_sp, p1_factors, p1_fps)
#--- Screen1_panel2_10U
dir_sp = root_dir + "/Screen1_panel2_10U"
header_sp = ["Well"]
header_idx_sp = [0]
find_AG(dir_sp, header_idx_sp, header_sp, p2_factors, p2_fps)
#--- Screen1_panel2_100U
dir_sp = root_dir + "/Screen1_panel2_100U"
header_sp = ["Well"]
header_idx_sp = [0]
find_AG(dir_sp, header_idx_sp, header_sp, p2_factors, p2_fps)
#--- Screen2_panel1_10U
dir_sp = root_dir + "/Screen2_panel1_10U"
header_sp = ["Well", "shRNA"]
header_idx_sp = [2, 1]
find_AG(dir_sp, header_idx_sp, header_sp, p1_factors, p1_fps)
#--- Screen2_panel1_100U
dir_sp = root_dir + "/Screen2_panel1_100U"
header_sp = ["Well", "shRNA"]
header_idx_sp = [2, 1]
find_AG(dir_sp, header_idx_sp, header_sp, p1_factors, p1_fps)
#--- Screen2_panel2_10U
dir_sp = root_dir + "/Screen2_panel2_10U"
header_sp = ["Well"]
header_idx_sp = [0]
find_AG(dir_sp, header_idx_sp, header_sp, p2_factors, p2_fps)
#--- Screen2_panel2_100U
dir_sp = root_dir + "/Screen2_panel2_100U"
header_sp = ["Well"]
header_idx_sp = [0]
find_AG(dir_sp, header_idx_sp, header_sp, p2_factors, p2_fps)
#--- Screen3_panel1_10U
dir_sp = root_dir + "/Screen3_panel1_10U"
header_sp = ["Well", "shRNA"]
header_idx_sp = [2, 1]
find_AG(dir_sp, header_idx_sp, header_sp, p1_factors, p1_fps)
#--- Screen3_panel1_100U
dir_sp = root_dir + "/Screen3_panel1_100U"
header_sp = ["Well", "shRNA"]
header_idx_sp = [2, 1]
find_AG(dir_sp, header_idx_sp, header_sp, p1_factors, p1_fps)
#--- Screen3_panel2_10U
dir_sp = root_dir + "/Screen3_panel2_10U"
header_sp = ["Well"]
header_idx_sp = [0]
find_AG(dir_sp, header_idx_sp, header_sp, p2_factors, p2_fps)  
#--- Screen3_panel2_100U
dir_sp = root_dir + "/Screen3_panel2_100U"
header_sp = ["Well"]
header_idx_sp = [0]
find_AG(dir_sp, header_idx_sp, header_sp, p2_factors, p2_fps)
'''


###----- Find all the Amtrine+ data
p1_factors = ["CD44", "CD25", "CD62L", "CXCR3", "CD103", "CD127"]
p2_factors = ["PD1", "Lag3", "Tim3"]
p1_fps = ["Alexa 700", "PerCP-Cy5-5", "BV605", "PE-Cy7", "APC", "PE"]
p2_fps = ["BV605", "PerCP-Cy5-5", "PE"]
'''
#--- Screen3_panel1_100U_CXCR3
dir_sp = root_dir + "/Screen3_panel1_100U_CXCR3/Screen3_panel1_100U"
header_sp = ["Well"]
header_idx_sp = [0]
find_A(dir_sp, header_idx_sp, header_sp, p1_factors, p1_fps)
#--- Screen1_panel1_10U
dir_sp = root_dir + "/Screen1_panel1_10U"
header_sp = ["Well", "shRNA"]
header_idx_sp = [2, 1]
find_A(dir_sp, header_idx_sp, header_sp, p1_factors, p1_fps)
#--- Screen1_panel1_100U
dir_sp = root_dir + "/Screen1_panel1_100U"
header_sp = ["Well", "shRNA"]
header_idx_sp = [2, 1]
find_A(dir_sp, header_idx_sp, header_sp, p1_factors, p1_fps)
#--- Screen1_panel2_10U
dir_sp = root_dir + "/Screen1_panel2_10U"
header_sp = ["Well"]
header_idx_sp = [0]
find_A(dir_sp, header_idx_sp, header_sp, p2_factors, p2_fps)
#--- Screen1_panel2_100U
dir_sp = root_dir + "/Screen1_panel2_100U"
header_sp = ["Well"]
header_idx_sp = [0]
find_A(dir_sp, header_idx_sp, header_sp, p2_factors, p2_fps)
#--- Screen2_panel1_10U
dir_sp = root_dir + "/Screen2_panel1_10U"
header_sp = ["Well", "shRNA"]
header_idx_sp = [2, 1]
find_A(dir_sp, header_idx_sp, header_sp, p1_factors, p1_fps)
#--- Screen2_panel1_100U
dir_sp = root_dir + "/Screen2_panel1_100U"
header_sp = ["Well", "shRNA"]
header_idx_sp = [2, 1]
find_A(dir_sp, header_idx_sp, header_sp, p1_factors, p1_fps)
#--- Screen2_panel2_10U
dir_sp = root_dir + "/Screen2_panel2_10U"
header_sp = ["Well"]
header_idx_sp = [0]
find_A(dir_sp, header_idx_sp, header_sp, p2_factors, p2_fps)
#--- Screen2_panel2_100U
dir_sp = root_dir + "/Screen2_panel2_100U"
header_sp = ["Well"]
header_idx_sp = [0]
find_A(dir_sp, header_idx_sp, header_sp, p2_factors, p2_fps)
#--- Screen3_panel1_10U
dir_sp = root_dir + "/Screen3_panel1_10U"
header_sp = ["Well", "shRNA"]
header_idx_sp = [2, 1]
find_A(dir_sp, header_idx_sp, header_sp, p1_factors, p1_fps)
#--- Screen3_panel1_100U
dir_sp = root_dir + "/Screen3_panel1_100U"
header_sp = ["Well", "shRNA"]
header_idx_sp = [2, 1]
find_A(dir_sp, header_idx_sp, header_sp, p1_factors, p1_fps)
#--- Screen3_panel2_10U
dir_sp = root_dir + "/Screen3_panel2_10U"
header_sp = ["Well"]
header_idx_sp = [0]
find_A(dir_sp, header_idx_sp, header_sp, p2_factors, p2_fps)  
#--- Screen3_panel2_100U
dir_sp = root_dir + "/Screen3_panel2_100U"
header_sp = ["Well"]
header_idx_sp = [0]
find_A(dir_sp, header_idx_sp, header_sp, p2_factors, p2_fps)
'''


###----- Compile data
root_dir = "/Users/yolandatiao/Desktop/CRF_Screen/Meghan_originaldata/Screen_markers"
factor_list = ["CD44", "CD25", "CD62L", "CXCR3", "CD103", "PD1", "Lag3", "Tim3"]
#factor_list = ["CD127"]

'''
#--- AmetGFP
for factor in factor_list:
    compileFactor(root_dir, factor, "10U")
    compileFactor(root_dir, factor, "100U")
#--- Amet
for factor in factor_list:
    compileFactor_Amt(root_dir, factor, "10U")
    compileFactor_Amt(root_dir, factor, "100U")
'''

'''
#--- CXCR3 
root_dir = "/Users/yolandatiao/Desktop/CRF_Screen/Meghan_originaldata/Screen_markers/Screen3_panel1_100U_CXCR3"
compileFactor(root_dir, "CXCR3", "100U")
root_dir = "/Users/yolandatiao/Desktop/CRF_Screen/Meghan_originaldata/Screen_markers/Screen3_panel1_100U_CXCR3"
compileFactor_Amt(root_dir, "CXCR3", "100U")
'''

###----- Standardaize data
#--- Amt
'''
wk_dir = "/Users/yolandatiao/Desktop/CRF_Screen/1_Raw/Amt"
os.chdir(wk_dir)

for file in glob.glob("Amt_*"):
    std_Amt(file)
'''

#--- Get all plate-well matching with shRNAs
'''
wk_dir = "/Users/yolandatiao/Desktop/CRF_Screen/1_Raw"
os.chdir(wk_dir)

plate_well_shRNA = []
for file in glob.glob("standard_*"):
    with open(file, "r") as fin:
        rfin = csv.reader(fin, delimiter = ",")
        next(rfin)
        for row in rfin:
            plate_well_shRNA.append("%s--%s" %(row[0], row[1].split("+")[0].replace("sh", "")))
                
plate_well_shRNA_set = list(set(plate_well_shRNA))
len(plate_well_shRNA_set)
plate_well_shRNA_set.sort()

with open("ref.csv", "w") as fout:
    wfout = csv.writer(fout)
    wfout.writerow(["plate_well", "shRNA"])
    for x in plate_well_shRNA_set:
        wfout.writerow(x.split("--"))
'''

###----- Create shRNA matched data sheets
'''
wk_dir = "/Users/yolandatiao/Desktop/CRF_Screen/1_0_Raw"
os.chdir(wk_dir)

for file in glob.glob("*CD127*"):
    shRNAmatch(file)
'''



