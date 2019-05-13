#!/bin/bash

cd /Volumes/Yolanda/CRF_Screen/InVivo/1_0_Raw/2_flt/Exp56

nu_list=(1 17 3 18 20 6 7 8 9 10 11 12 13 14 15)
name_list=(P1-7_Input P1-7_Q1 P1-7_Q2 P1-7_Q3 P1-7_Q4 P8-14_Input P8-14_Q1 P8-14_Q2 P8-14_Q3 P8-14_Q4 P15-21_Input P15-21_Q1 P15-21_Q2 P15-21_Q3 P15-21_Q4)

for i in $(seq 0 14)
do
  mv BC_${nu_list[i]}_count_flt.csv ${name_list[i]}.csv
done