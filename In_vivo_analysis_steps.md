## In vitro screen data analysis
### 0. Data compilation <br/>
* __0.1__ Count collection <br/>
* __0.2__ File inspection and filter (target constructs v.s. control constructs v.s. others according to reference file) <br/>
* __0.3__ Rename files from barcode number to sample name <br/>

>__Script__: [0_0_count_collection.py](0.1_Codes_Invivo/0_0_count_collection.py); [0_1_file_inspection_filter.py](0.1_Codes_Invivo/0_1_file_inspection_filter.py); [0_2_rename.sh](0.1_Codes_Invivo/0_2_rename.sh) <br/>
__Source files__: /InVivo/1_0_Raw/CRF_Screen_raw/ (Not uploaded) <br/>
__Output files__: <br/>
[Counts](InVivo/1_0_Raw/1_count) <br/>
[Non target filtered counts](InVivo/1_0_Raw/2_flt) <br/>
[Exp35 reads summary](InVivo/1_0_Raw/Exp35_reads_summary.csv) <br/>

### 1. Data normalization and conversion <br/> 
* __1.0__ Filter outliers by Z score <br/>
* __1.1__ Combine Exp35 and Exp56 data <br/>
* __1.2__ Filter bench contaminant constructs <br/>
* __1.3__ Data conversion: <br/>
    * Calculate percentages of counts in each group <br/>
    * Combine percentage of `counts * million` in each file to create total distribution <br/> 
    * Calculate percentile of counts in total distribution (negative binomial distribution) <br/>
    * Calculate gate percentile shift (e.g. for shX: Q4-Q1, Q3-average_of_others, etc.) <br/>
    * Compile results of different pools <br/>
    * Calculate z score <br/>
    
>__Script__: [Normalization and conversion](0.1_Codes_Invivo/1_0_shRNAlibrary_analysis_0513_Exp35Exp56_CombinedRawCount_nbPctg.py); [Z-Score](0.1_Codes_Invivo/1_1_Zscore_of_Pctg_shift.py) <br/>
__Source files__: [Non target filtered counts](InVivo/1_0_Raw/2_flt) <br/>
__Output files__: <br/>
[Outlier filtered counts](InVivo/1_1_Norm/20190513_Exp35Exp56_nbPctl-All/0_fltOutlier) <br/>
[Exp35 Exp56 combined - outlier filtered counts](InVivo/1_1_Norm/20190513_Exp35Exp56_nbPctl-All/1_Exp35Exp56_combined) <br/>
[Exp35 Exp56 combined - outlier and contaminant filtered counts](InVivo/1_1_Norm/20190513_Exp35Exp56_nbPctl-All/2_flt_comtaminants) <br/>
[Gate comparisons - seperated each pool](InVivo/1_1_Norm/20190513_Exp35Exp56_nbPctl-All/3_gate_comparisons_bypool) <br/>
[Gate comparisons - compiled + Z-Score](InVivo/1_1_Norm/20190513_Exp35Exp56_nbPctl-All/4_gate_comparisons_combined)


    
