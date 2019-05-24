#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 17:26:07 2019

@author: yolandatiao
"""

########################################  GOterm annotation ######################################## 
# Author: Huitian (Yolanda) Diao
# May 20th, 2019

########################################  Import ######################################## 
import os
import csv
from astropy.io import ascii
from astropy.table import Table, join, vstack
import math

########################################  Self-defined functions ########################################
def foundWord(list_x, strx):
    for x in list_x:
        if x in strx:
            return(True)
    return(False)

########################################  Main ########################################
###----- Select CRFs from GO-term reference
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/GO_terms"
os.chdir(wk_dir)

CRM_list_file = "/Volumes/Yolanda/RNAseq_Compilation/CRM_analysis_2/Info/CRM.csv"
CRM_tab = ascii.read(CRM_list_file)
all_CRM = list(CRM_tab.columns[0]) + list(CRM_tab.columns[1])
all_CRM = [x for x in all_CRM if (x != 'NA' and x != "")]
all_CRM = [x.upper() for x in all_CRM]

in_go_file = "/Volumes/Yolanda/CRF_Screen/Ref/mouse_go_ref_20190509.txt"
out_go_file = "mouse_go_CRF.csv"

with open(in_go_file, "r") as fin:
    with open(out_go_file, "w") as fout:
        rfin = csv.reader(fin, delimiter="\t")
        wfout = csv.writer(fout, delimiter=",")
        header = next(rfin)
        wfout.writerow(header)
        for row in rfin:
            if row[2].upper() in all_CRM:
                wfout.writerow(row)

###----- Match GO ID with name
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/GO_terms"
os.chdir(wk_dir)

go_ref_file = "/Volumes/Yolanda/CRF_Screen/Ref/mouse_go_obo_20190509.csv"
go_ref_tab = ascii.read(go_ref_file)
go_ref_tab.columns[0].name = "GO_ID"

crm_go_file = "mouse_go_CRF_simp.csv"
crm_go_tab = ascii.read(crm_go_file)
crm_go_tab.colnames

crm_go_tab = join(crm_go_tab, go_ref_tab, join_type="left", keys = "GO_ID")

ascii.write(crm_go_tab, "mouse_go_CRF_names.csv", format="csv", overwrite=True)


###----- Filter 
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/GO_terms"
os.chdir(wk_dir)

in_file = "mouse_go_CRF_names.csv"
out_file_1 = "CRF_all_processes.csv"

in_tab = ascii.read(in_file)
in_names = list(in_tab["name"])
in_def = list(in_tab["def"])
in_namespace = list(in_tab["namespace"])
in_alldescp = ["%s:::%s:::%s"%(x,y,z) for index, (x, y, z) in enumerate(zip(in_namespace, in_names, in_def))]
in_alldescp_set = list(set(in_alldescp))
with open(out_file_1, "w") as fout:
    wfout = csv.writer(fout, delimiter=",")
    wfout.writerow(["namespace", "name", "def"])
    for descp in in_alldescp_set:
        wfout.writerow(descp.split(":::"))

#--- Find all words appeared in the process names
# And count their appearance
in_file2 = "CRF_all_processes.csv"
out_file_2 = "CRF_all_processes_name_word_count.csv"

in_tab2 = ascii.read(in_file2)
in_names2 = list(in_tab2["name"])
in_names2_allchar = []
for x in in_names2:
    in_names2_allchar = in_names2_allchar + x.split(" ")
in_names2_allchar_set = list(set(in_names2_allchar))
in_names2_allchar_set_count = [in_names2_allchar.count(x) for x in in_names2_allchar_set]

with open(out_file_2, "w") as fout:
    wfout = csv.writer(fout, delimiter=",")
    wfout.writerow(["word", "count"])
    for index, (x, y) in enumerate(zip(in_names2_allchar_set, in_names2_allchar_set_count)):
        wfout.writerow([x, y])

#--- Filter selected processes
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/GO_terms"
os.chdir(wk_dir)

ref_file = "CRF_all_processes_name_word_count_anno.csv"
ref_tab = ascii.read(ref_file)
use_list = list(ref_tab["use"])
word_list = list(ref_tab["word"])
word_use_list = [y for index, (x, y) in enumerate(zip(use_list, word_list)) if x != ""]

in_file = "mouse_go_CRF_names.csv"
out_file_3 = "mouse_go_CRF_names_slt.csv"
with open(in_file, "r") as fin:
    with open(out_file_3, "w") as fout:
        rfin = csv.reader(fin, delimiter=",")
        wfout = csv.writer(fout, delimiter=",")
        header = next(rfin)
        wfout.writerow(header)
        for row in rfin:
            if foundWord(word_use_list, row[3]):
                wfout.writerow(row)

#--- Filter out processes that are too general:
# GO:0000122	 biological_process	negative regulation of transcription by RNA polymerase II
# GO:0000083 biological_process	regulation of transcription involved in G1/S transition of mitotic cell cycle
# GO:0033613	 molecular_function	activating transcription factor binding
# GO:0043044	 biological_process	ATP-dependent chromatin remodeling
# GO:0000785	 cellular_component	chromatin
# GO:0031497 biological_process	chromatin assembly
# GO:0006333 biological_process	chromatin assembly or disassembly
# GO:0003682	 molecular_function	chromatin binding
# GO:0031490 molecular_function	chromatin DNA binding
# GO:0006325 biological_process	chromatin organization
# GO:0048096	 biological_process	chromatin-mediated maintenance of transcription
# GO:0005694 cellular_component	chromosome
# GO:0051276 biological_process	chromosome organization
# GO:0070192 biological_process	chromosome organization involved in meiotic cell cycle
# GO:0000775 cellular_component	chromosome, centromeric region
# GO:0000781	 cellular_component	chromosome, telomeric region
# GO:0000794	 cellular_component	condensed nuclear chromosome
# GO:0000780	 cellular_component	condensed nuclear chromosome, centromeric region
# GO:0008340	 biological_process	determination of adult lifespan
# GO:0003700	 molecular_function	DNA-binding transcription factor activity
# GO:0000792 cellular_component	heterochromatin
# GO:0042393	 molecular_function	histone binding
# GO:0000122	 biological_process	negative regulation of transcription by RNA polymerase II
# GO:0045892	 biological_process	negative regulation of transcription, DNA-templated
# GO:0016604	 cellular_component	nuclear body
# GO:0000790	 cellular_component	nuclear chromatin
# GO:0008134	 molecular_function	transcription factor binding
# GO:0016740	 molecular_function	transferase activity

wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/GO_terms"
os.chdir(wk_dir)

flt_file = "Filter_list.txt"
flt_tab = ascii.read(flt_file)
flt_ids = list(flt_tab['GO_ID'])

in_file = "mouse_go_CRF_names_slt.csv"
out_file_4 = "mouse_go_CRF_names_slt_flt.csv"
with open(in_file, "r") as fin:
    with open(out_file_4, "w") as fout:
        rfin = csv.reader(fin, delimiter=",")
        wfout = csv.writer(fout, delimiter=",")
        header = next(rfin)
        wfout.writerow(header)
        for row in rfin:
            if row[1] not in flt_ids:
                wfout.writerow(row)
                
##########-------------------- Find the best cluster annotation for all CRFs
wk_dir = "/Volumes/Yolanda/CRF_Screen/InVivo/2_Protein/GO_terms"
os.chdir(wk_dir)


###----- Annotate complexes and count complex appearance
ref_file = "mouse_go_CRF_names_slt_flt.csv"
ref_tab = ascii.read(ref_file)
ref_tab_gn = list(ref_tab["DB_Object_Symbol"])
ref_tab_name = list(ref_tab["name"])

all_crf_file = "/Volumes/Yolanda/RNAseq_Compilation/CRM_analysis_2/Info/CRM.csv"
out_file_5 = "CRM_complexes.csv"

anno_keys = ["methyl", "acetyl", "phospho", "helicase", "SMAD protein complex assembly",
             "MAPK", "RNA polymerase II", "oxygenase", "deiminase", "ribosyltransferase",
             "complex", "mRNA 5'-splice site recognition", "ubiquit", "succinylation", "malonylation"]

all_complex_list = []

with open(all_crf_file, "r") as fin:
    with open(out_file_5, "w") as fout:
        rfin = csv.reader(fin, delimiter=",")
        wfout = csv.writer(fout, delimiter=",")
        wfout.writerow(next(rfin))
        for row in rfin:
            gn = row[0]
            newrow = row[:2]
            complex_list = []
            for idx, (refGN, refGoname) in enumerate(zip(ref_tab_gn, ref_tab_name)):
                if refGN == gn:
                    appendx = False
                    for key in anno_keys:
                        if key in refGoname:
                            appendx = True
                    if appendx:
                        complex_list.append(refGoname)
            complex_list = list(set(complex_list))
            newrow = newrow + complex_list
            all_complex_list += complex_list
            wfout.writerow(newrow)

all_complex_set = list(set(all_complex_list))
all_complex_set_count = [all_complex_list.count(x) for x in all_complex_set]

with open("complex_count.txt", "w") as fout:
    wfout = csv.writer(fout, delimiter = "\t")
    wfout.writerow(["complex", "count"])
    for index, (x, y) in enumerate(zip(all_complex_set, all_complex_set_count)):
        wfout.writerow([x,y])

###----- Annotate complexes and count complex appearance
all_crf_file = "CRM_complexes.csv"
out_file_6 = "CRM_complexes_rank.csv"

complex_ref = "complex_count_rank_use_true.txt"
complex_ref_tab = ascii.read(complex_ref)
complex_ref_tab_name = list(complex_ref_tab["complex"])

with open(all_crf_file, "r") as fin:
    with open(out_file_6, "w") as fout:
        rfin = csv.reader(fin, delimiter=",")
        wfout = csv.writer(fout, delimiter=",")
        wfout.writerow(next(rfin))
        for row in rfin:
            newrow = row[:2]
            complex_list = row[2:]
            for name in complex_ref_tab_name:
                if name in complex_list:
                    newrow.append(name)
            wfout.writerow(newrow)


###----- Create complex v.s. gene file
in_file = "CRM_complexes_rank.csv"
out_file_7 = "CRM_complexes_count.csv"

complex_count_tab = ascii.read("complex_count_rank_use_true.txt")
complexes = list(complex_count_tab["complex"])

with open(in_file, "r") as fin:
    with open(out_file_7, "w") as fout:
        rfin = csv.reader(fin, delimiter=",")
        wfout = csv.writer(fout, delimiter=",")
        header = ["gene_name"] + complexes
        wfout.writerow(header)
        next(rfin)
        for row in rfin:
            new_row = [row[0]]
            row_complexes = row[2:]
            for i in complexes:
                if i in row_complexes:
                    new_row.append("Yes")
                else:
                    new_row.append("--")
            wfout.writerow(new_row)



