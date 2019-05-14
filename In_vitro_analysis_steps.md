## Data Analysis Steps
### 0. Data compilation <br/> 
* __0.1__ Identify Amtrine+ / Amtrine+GFP+ columns (transduced population) in original spreadsheet <br/>
* __0.2__ Extract GeoMean / Percentage data for each marker for Amtrine+ / Amtrine+GFP+ <br/> 
* __0.3__ Write output files for each individual marker including: __plate and well position__ and __Flow readout (MFI or percentage)__ <br/> 
* __0.4__ Match shRNA names to position <br/>

>__Script__: [0_dataPrep.py](0_Codes/0_dataPrep.py) <br/> 
__Source files__: [Screen_markers](InVitro/Megan_originaldata/Screen_markers)  <br/> 
__Output files__: [1_0_Raw](InVitro/1_0_Raw); [1_1_shRNAmatched](InVitro/1_1_shRNAmatched) <br/>

### 1. Data normalization <br/> 
* __1.1__ Calculate percentile and Z-score of each sample <br/> 
   * _1.1.1_ Approach 1: Calculate normal distribution percentile and Z-score based on total population <br/> 
   * _1.1.2_ Approach 2: Calculate normal distribution percentile and Z-score based on control plate <br/> 
* __1.2__ Calculate stats of un-normalized data / normalized data (mean and standard deviation) <br/> 
>__Script__: [1_dataConversion.R](0_Codes/1_dataConversion.R) <br/> 
__Source files__: [1_1_shRNAmatched](InVitro/1_1_shRNAmatched) <br/> 
__Output files__: [1_2_normtoall_ZP](InVitro/1_2_normtoall_ZP); [1_2_normtocontrol_ZP/Amt](InVitro/1_2_normtocontrol_ZP/0_Amt); [1_2_normtocontrol_ZP/Amt](InVitro/1_2_normtocontrol_ZP/0_AmtGFP) <br/>

### 2. Per gene shRNA effect for each marker (by t-test)
* __2.1__ Test consistency of effect of shRNAs targeting the same gene for each marker <br/> 
   * _2.1.1_ Calculate p-value of _shRNAs targeting the same gene_ v.s. _all samples_ with two sample t-test (<- from percentile calculated in last step)  <br/> 
   * _2.1.2_ Calculate average Z score and percentile of _shRNA targeting the same gene_  <br/> 
* __2.2__ Compile data into spreadsheets to summarize _average Z score_, _percentile_ and _p-value_ of each gene for different markers
* __2.3__ Plot gene average Z-score heatmaps for each marker
>__Script__: [2_convertedDataAnalysis.R](0_Codes/2_convertedDataAnalysis.R) <br/> 
__Source files__: [1_2_normtocontrol_ZP/Sep](InVitro/1_2_normtocontrol_ZP/Sep) <br/> 
__Output files__: [2_0_t-test_by_gene/0_Amt_sep](InVitro/2_0_t-test_by_gene/0_Amt_sep); [2_0_t-test_by_gene/0_AmtGFP_sep](InVitro/2_0_t-test_by_gene/0_AmtGFP_sep); [2_0_t-test_by_gene/1_Amt_compile](InVitro/2_0_t-test_by_gene/1_Amt_compile); [2_0_t-test_by_gene/1_AmtGFP_compile](InVitro/2_0_t-test_by_gene/1_AmtGFP_compile)<br/>

### 3. Per gene shRNA effect for all markers (by t-sne)
* __3.1__ Cluster shRNAs by Z-score for each marker with t-sne clustering
* __3.2__ Examine clustering of shRNAs targeting the same gene with the t-sne clustering result
>__Script__: [3_Scikit_learn_byshRNA.py](0_Codes/3_Scikit_learn_byshRNA.py); [3_tsne_cbs-plot_byshRNA.R](0_Codes/3_tsne_cbs-plot_byshRNA.R) <br/> 
__Source files__: [1_2_normtocontrol_ZP/Compiled_* ](InVitro/1_2_normtocontrol_ZP) <br/> 
__Output files__: [1_2_normtocontrol_ZP/Compiled_* /per*_ clusterbyshRNA](InVitro/1_2_normtocontrol_ZP) (per*: perplexity) <br/>

### 4. Gene knockdown effect clustering (by t-sne)
__Aim__: <br/>
* __4.1__ Cluster genes by average Z-score for each marker with t-sne clustering (per*: perplexity)
* __4.2__ Plot genes which significantly up/down regulate markers when knocked down (by marker)
>__Script__: [3_Scikit_learn_byGene.py](0_Codes/3_Scikit_learn_byGene.py); [3_tsne_cbs-plot_byGene.R](0_Codes/3_tsne_cbs-plot_byGene.R) <br/> 
__Source files__: [2_0_t-test_by_gene/1_Amt_compile](InVitro/2_0_t-test_by_gene/1_Amt_compile) <br/> 
__Output files__: [t-sne coordinates](InVitro/2_0_t-test_by_gene/1_Amt_compile/Amt_normbycontrolZP_t-test.by.geneavg_z-score_tsne_per5.csv); [bubble plot source](nVitro/2_0_t-test_by_gene/2_Amt_compile_cluster/0_bbplot_source); [annotated bubble plots](InVitro/2_0_t-test_by_gene/2_Amt_compile_cluster/1_bbplot_anno); [non-annotated bubble plots](InVitro/2_0_t-test_by_gene/2_Amt_compile_cluster/1_bbplot_plain) <br/>
