### The Index matrix
- The 1st column is the cCRE chromosome;
- The 2nd column is the cCRE start_position;
- The 3rd column is the cCRE end_position;
- The 4th column is the cCRE IDs.
- The following columns are the binary index of each cCRE at each cell-type based on the cCRE presence and absence calls inside a certain cell-type.
```
>>> head input_data/Input_matrix/snapshot_test_run_merge.index.matrix.txt
#chr	start	end	cCREsID	HSC	LMPP	MPP	CMP	MEP	ERY	GMP	MONp	CLP	B	NK	TCD4	TCD8
chr1	817200	817400	1	0	0	0	0	0	0	1	0	0	0	0	0	0
chr1	827400	827600	2	1	0	1	1	0	1	1	1	0	1	1	1	1
chr1	842000	842400	3	0	0	0	0	0	0	0	0	0	0	0	1	1
chr1	842800	843200	4	0	0	0	0	0	1	0	0	0	0	0	0	0
chr1	858000	858200	5	0	0	0	0	0	0	0	0	0	0	0	1	0
chr1	865600	866000	6	0	0	0	0	0	0	0	0	0	1	0	0	0
chr1	869800	870000	7	0	0	0	0	0	1	0	0	0	1	1	1	1
chr1	897400	897600	8	1	0	1	0	0	0	0	0	0	0	0	0	0
chr1	904600	904800	9	0	0	0	0	0	0	0	0	0	0	0	1	1
```

### The Signal matrix
- The 1st column is the cCRE chromosome;
- The 2nd column is the cCRE start_position;
- The 3rd column is the cCRE end_position;
- The 4th column is the cCRE IDs.
- The following columns are the numeric signal of each cCRE at each cell-type.
```
head input_data/Input_matrix/snapshot_test_run_merge.signal.matrix.txt
#chr	start	end	cCREsID	HSC	LMPP	MPP	CMP	MEP	ERY	GMP	MONp	CLP	B	NK	TCD4	TCD8
chr1	817200	817400	1	0.609961	0.109752	0.344694	0.897595	0.344694	0.109752	2.03892	1.41998	0.109752	0.344694	0.344694	0.344694	0.227223
chr1	827400	827600	2	2.03892	0.727432	1.87953	3.30794	0.685021	4.08286	2.03892	6.24662	0.344694	6.71504	4.06742	3.83032	6.54589
chr1	842000	842400	3	0.109752	0.109752	0.109752	0.109752	0.109752	0.227223	0.168488	0.514857	1.12738	0.37104	0.227223	2.85771	1.90568
chr1	842800	843200	4	0.168488	0.109752	0.397386	0.109752	0.227223	2.80997	0.227223	0.227223	0.227223	0.397386	0.168488	0.456122	0.285959
chr1	858000	858200	5	0.227223	0.397386	0.109752	1.14458	1.11017	0.397386	0.344694	0.109752	0.227223	1.24981	0.344694	1.93183	0.727432
chr1	865600	866000	6	0.109752	0.609961	0.109752	0.109752	0.344694	0.109752	0.359857	0.109752	0.109752	4.10623	0.227223	0.168488	0.109752
chr1	869800	870000	7	0.727432	0.397386	0.514857	1.41998	0.344694	2.77387	0.974417	1.86875	0.685021	2.94155	2.15493	6.02499	4.20888
chr1	897400	897600	8	2.17881	1.35716	2.17881	0.727432	0.344694	0.397386	0.514857	0.109752	0.344694	0.344694	0.227223	0.227223	0.109752
chr1	904600	904800	9	0.609961	0.344694	0.685021	0.514857	0.109752	0.974417	0.974417	0.856946	0.227223	1.41998	1.87953	3.93349	4.72233
```

### The Functional Epigenetic State matrix
- The 1st column is the cCRE chromosome;
- The 2nd column is the cCRE start_position;
- The 3rd column is the cCRE end_position;
- The 4th column is the cCRE IDs.
- The following columns are the Functional Epigenetic State label of each cCRE at each cell-type .
```
head input_data/Input_matrix/snapshot_test_run_merge.function.matrix.txt
#chr	start	end	cCREsID	HSC	LMPP	MPP	CMP	MEP	ERY	GMP	MONp	CLP	B	NK	TCD4	TCD8
chr1	817200	817400	1	0	0	0	0	0	0	0	21	0	0	0	6	0
chr1	827400	827600	2	15	15	15	15	15	15	15	15	15	15	15	15	15
chr1	842000	842400	3	0	0	0	0	0	5	0	0	5	0	0	5	19
chr1	842800	843200	4	2	2	2	2	2	6	2	0	2	2	2	6	6
chr1	858000	858200	5	0	0	0	0	13	0	0	0	0	9	0	9	9
chr1	865600	866000	6	3	4	3	3	3	3	3	3	3	16	0	3	0
chr1	869800	870000	7	9	9	9	9	9	13	9	13	9	24	13	13	13
chr1	897400	897600	8	18	18	18	3	11	11	3	3	3	3	11	11	11
chr1	904600	904800	9	11	11	11	11	11	12	11	24	11	24	24	13	24
```


### The Cell Type Differentiation Tree
Each row in the cell type differentiation tree file represents one edge in the cell type differentiation tree. 
The format of the file is as follows: 
- The first column represents the progenitor cell type
- The second column represents the differentiated cell type
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
- Column 1: Epigenetic state label
- Column 2: RGB color of the state
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
