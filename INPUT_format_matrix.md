### The cell type peak & signal & epigenetic-state file list:
The 1st column is the cell type label;
The 2nd column is the cell-type specific peak bed absolute-file-path;
The 3rd column is the cell-type specific signal bigWig absolute-file-path;
The 4th column is the cell-type specific functional epigenetic state bigBed absolute-file-path.

```
>>> cat input_data/peak_signal_state_list.txt 
HSC	atac_pk/HSC.merge.sort.tabNBPfdrpk.bed	atac_sig/HSC.merge.bw	function_label/HSC_100258.state.bb
LMPP	atac_pk/LMPP.merge.sort.tabNBPfdrpk.bed	atac_sig/LMPP.merge.bw	function_label/LMPP_100268.state.bb
MPP	atac_pk/MPP.merge.sort.tabNBPfdrpk.bed	atac_sig/MPP.merge.bw	function_label/MPP_100272.state.bb
CMP	atac_pk/CMP.merge.sort.tabNBPfdrpk.bed	atac_sig/CMP.merge.bw	function_label/CMP_100246.state.bb
MEP	atac_pk/MEP.merge.sort.tabNBPfdrpk.bed	atac_sig/MEP.merge.bw	function_label/MEP_Donor2596.state.bb
ERY	atac_pk/ERY.merge.sort.tabNBPfdrpk.bed	atac_sig/ERY.merge.bw	function_label/ERY_S002R5.state.bb
GMP	atac_pk/GMP.merge.sort.tabNBPfdrpk.bed	atac_sig/GMP.merge.bw	function_label/GMP_100256.state.bb
MONp	atac_pk/MONp.merge.sort.tabNBPfdrpk.bed	atac_sig/MONp.merge.bw	function_label/MONp_Prim_mon_C.state.bb
CLP	atac_pk/CLP.merge.sort.tabNBPfdrpk.bed	atac_sig/CLP.merge.bw	function_label/CLP_100266.state.bb
B	atac_pk/B.merge.sort.tabNBPfdrpk.bed	atac_sig/B.merge.bw	function_label/B_B15_50.state.bb
NK	atac_pk/NK.merge.sort.tabNBPfdrpk.bed	atac_sig/NK.merge.bw	function_label/NK_S005YG.state.bb
TCD4	atac_pk/TCD4.merge.sort.tabNBPfdrpk.bed	atac_sig/TCD4.merge.bw	function_label/T_CD4_S008H1.state.bb
TCD8	atac_pk/TCD8.merge.sort.tabNBPfdrpk.bed	atac_sig/TCD8.merge.bw	function_label/T_CD8_C0066PH1.state.bb
```

### The Cell Type Differentiation Tree
Each row in the cell type differentiation tree file represents one edge in the cell type differentiation tree. 
The format of the file is as follows: 
1. The first column represents the progenitor cell type
2. The second column represents the differentiated cell type
Here's a sample of the file content: 
```
>>> cat input_data/cd_tree.txt
HSC,LMPP
HSC,MPP
HSC,CMP
CMP,MEP
MEP,ERY
CMP,GMP
GMP,MONp
HSC,CLP
CLP,NK
CLP,B
CLP,TCD4
CLP,TCD8
```

## The functional state color list
Each row in the functional state color list represents the color assigned to a specific epigenetic state label. 
The following is the structure of the list:
1. Column 1: Epigenetic state label
2. Column 2: RGB color of the state
Here's a sample of the file content:
```
>>> cat input_data/function_color_list.txt 
0	255,255,255
1	232,232,232
2	195,230,195
3	174,174,245
4	253,253,161
5	240,199,226
6	250,250,0
7	95,95,95 
8	91,187,90
9	231,142,253
10	253,222,174
11	0,0,211
12	255,44,43
13	188,0,248
14	253,0,0
15	246,0,0
16	248,211,0
17	245,243,0
18	0,0,202
19	190,0,91 
20	249,154,0
21	211,68,106 
22	216,76,211 
23	228,0,0
24	168,0,235
```
