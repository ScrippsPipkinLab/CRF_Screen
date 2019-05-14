# CRF_Screen
![Cell](https://i.pinimg.com/originals/7b/70/de/7b70dee0342490ca05c8f2e72b1d9cbc.jpg)

#### Contributors: <br/>
> Runqiang Chen: Conducted in vivo screen <br/>
> Megan Fridrick: Conducted in vitro screen <br/>
> Huitian Diao (Yolanda): Binformatics <br/>
> Matthew Pipkin: Advisor

## External data storage <br/>
* __In vitro screen analysis__: <br/>
[Ishtar folder](/Volumes/pipkinlab/Lab Members/Former members/Megan Frederick/MAF Screen Analysis)
* __In vivo screen analysis__: <br/>
[Drop box link](https://www.dropbox.com/sh/8vgdzzyc4w94z5q/AACK0rPymadt8TE3ot2Z0rzBa?dl=0)
* __Raw data__: <br/>
  * In vitro: <br/>
  /Volumes/pipkinlab/Lab Members/Former members/Megan Frederick/MAF Screen Data <br/>
  * In vivo: <br/>
  /Volumes/pipkinNGS/exp035 <br/>
  /Volumes/pipkinNGS/exp056 

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
[Exp35 barcode](InVivo/sample_barcode_exp35.csv)
[Exp56 barcode](InVivo/sample_barcode_exp56.csv)
[Target constructs](InVivo/shRNA_control.csv)
[Control constructs](InVivo/shRNA_ref.csv)

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



