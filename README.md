# CRF_Screen
![Cell](https://i.pinimg.com/originals/7b/70/de/7b70dee0342490ca05c8f2e72b1d9cbc.jpg)

#### Contributors: <br/>
> Runqiang Chen: Conducted in vivo screen <br/>
> Megan Fridrick: Conducted in vitro screen <br/>
> Huitian Diao (Yolanda): Binformatics <br/>
> Matthew Pipkin: Advisor

## Manuscript progress
[Chromatin remodeling atlas - GSuite](https://drive.google.com/drive/folders/1lQkHaRpWIaQ0S_ir95EF-ODGVYXxwww4?usp=sharing)
[Chromatin remodeling atlas - Dropbox](https://www.dropbox.com/sh/lrswxf2msgenqcj/AADE3R-FuQcxOk59wkrtzQ5Ja?dl=0)

## Data analysis description <br/>
* [In vitro analysis](In_vitro_analysis_steps.md)
* [In vivo analysis](In_vivo_analysis_steps.md)

## External data storage <br/>
* __In vitro screen analysis__: <br/>
[Ishtar folder](/Volumes/pipkinlab/Lab Members/Former members/Megan Frederick/MAF Screen Analysis)
* __In vivo screen analysis__: <br/>
[Drop box link](https://www.dropbox.com/sh/8vgdzzyc4w94z5q/AACK0rPymadt8TE3ot2Z0rzBa?dl=0)
* __Raw data__: <br/>
  * In vitro: <br/>
  [Ishtar folder](/Volumes/pipkinlab/Lab Members/Former members/Megan Frederick/MAF Screen Data) <br/>
  * In vivo: <br/>
  [Ishtar folder](/Volumes/pipkinNGS/exp035) <br/>
  [Ishtar folder](/Volumes/pipkinNGS/exp056) 

## Experiment design </br>
* __In vitro screen__: <br/>

| Fluorephore | Panel 1 | Panel 2|
| --- | --- | --- |
| AF700 | CD44 | n.a. |
| PerCP-Cy5.5 | CD25 | Lag3 |
| BV605 | CD62L | PD1 |
| PE-Cy7 | CXCR3 | n.a. |
| APC | CD103 | n.a. |
| PE | CD127 | Tim3 |

* __In vivo screen__: <br/>
[Exp35 barcode](InVivo/sample_barcode_exp35.csv) </br>
[Exp56 barcode](InVivo/sample_barcode_exp56.csv) </br>
[Target constructs](InVivo/shRNA_control.csv) </br>
[Control constructs](InVivo/shRNA_ref.csv) </br>

## Fixed Glitches
### In vitro screen - Missing wells
1. InVitro/Megan_originaldata/Screen_markers/Screen3_panel1_10U/16 10U p1.csv <br/>
Fixed shifted cells: **G10 G11** left shifted 2 cells
2. InVitro/Megan_originaldata/Screen_markers/Screen3_panel1_100U/19 100U p1.cs <br/>
Fixed shifted cells: **G10 G11** left shifted 2 cells
3. InVitro/Megan_originaldata/Screen_markers/Screen2_panel1_10U/14 10U p1.csv <br/>
The missing wells are actually **B2 - B8** instead of **G5 - G11**

### In vivo screen - Miss placed constructs
1. Conflict references: <br/>
__Incorrect refrence:__  <br/>
/InVivo/originalAnlysis/CR in vivo Replicate 2/shRNA Pools.xlsx <br/>
__Correct refrence:__  <br/>
[Dropbox folder](https://www.dropbox.com/sh/nukl01bprlml2e9/AABasCVaJv_7AodnUUD2h6r3a?dl=0)

## References:
1. [BioGrid protein interaction](https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-3.5.172/)

| Species | Taxonomy ID |
| --- | --- |
| Homo sapiens | 9606 |
| Mus musculus | 10090 |
| Rattus norvegicus | 10116 |



