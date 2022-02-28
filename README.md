# Snapshot
### Epigenetic modification of chromatin plays a pivotal role in regulating gene expression during cell differentiation. The scale and complexity of epigenetic data pose significant challenges for biologists to identify the regulatory events controlling cell differentiation. Here, we present a new method, called Snapshot, that uses epigenetic data to generate a hierarchical visualization for DNA regions with epigenetic features segregating along any given cell differentiation hierarchy of interest. Different hierarchies of cell types may be used to highlight the epigenetic history specific to any particular cell lineage. We demonstrate the utility of Snapshot using data from the VISION project, an international project for ValIdated Systematic IntegratiON of epigenomic data in mouse and human hematopoiesis.

## Reference
#### Xiang, G., Giardine, B., An, L., Sun, C., Keller, C., Heuston, E., Bodine, D., Hardison, R. and Zhang, Y., 2018. Snapshot: clustering and visualizing epigenetic history during cell differentiation. bioRxiv, p.291880.
######(https://www.biorxiv.org/content/early/2018/04/09/291880)


## The overall workflow of Snapshot

#### Figure 1
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/f1.png" width="800"/>

##### The overall workflow of snapshot. (a) In the 1st step, we use the binarized presence/absence status of ccREs at each location across all cell types to create a ccRE index to represent the unique combinatorial pattern of ccREs at the location. All ccREs with the same index are grouped into one cluster called an index-set (IS). (b) We select the threshold to filter the less prevalent ISs based on the distribution of the number of ccREs per IS (red dash line). (c) We visualize the ISs and their ccREs in a heatmap. The blue dash box represents one example of IS: 0_1_0_0_1. The filtered ISs and their ccREs are shown in the end of the heatmap (black dash box). (d) We apply the QDA model to correct ccRE’s index by borrowing information from the signals of the ccRE across all cell types. After the QDA step, some of the ccREs grouped in the rare ISs can be rescued. (e) Finally, we visualize the ccRE clusters in a heatmap. Each row in the heatmap is the ccRE pattern for each index-set, and each column is a cell type. The ccRE patterns are sorted by their indices in the heatmap. By our definition of the ccRE index, the ccRE patterns are separated if they have different ccRE status in a cell type; conversely they are clustered together if they have similar ccRE status in a cell type. These different plots will together highlight the epigenetic activity across cell types and the associated functional annotations and their enrichments for each index-set, which enhances the interpretation of the functional roles of each ccRE cluster during cell differentiation.

## Example of Snapshot output: 
### Hematopoietic cell differentiation in VISION (ValIdated Systematic IntegratiON of hematopoietic epigenomes) project
#### Figure 2
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/f2.png" width="800"/>

##### The heatmap of index-sets. (a) The heatmap of index-set colored by the average ATAC-seq signal in each cell type. (b) The heatmap of index-set colored by the most frequent functional annotation in each cell type. (c) The density plot of the number of genomic region covered by the index-set. (d) The color code and epigenetic composition of functional annotation used in (b).

#### Figure 3
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/f3.png" width="800"/>

##### The data visualization for index-set-149 and corresponding GO analysis and MEME-ChIP TF binding motif analysis. (a) The hematopoietic cell differentiation tree colored by the average ATAC-seq signal in each cell type of the index-set-149. The violin plot represents the distribution of ATAC-seq signal in each cell type of the index-set-149 is in below.  (b) The same cell differentiation tree colored by the most frequent functional annotation in each cell type of the index-set-149. The two most frequence functional annotation in erythroblasts lineage. The bar plot based on the proportion of each functional annotation in each cell type of the index-set-149 is below the cell differentiation tree. (c) The index-set-149 relevant GO term. (d) The index-set-149 significantly enriched TF binding motif from MEME-ChIP analysis.



## Dependence:
#### Python
###### numpy
###### sklearn
#### R
###### ggplot2; pheatmap; igraph; networkD3
#### bedtools

## Install Snapshot
#### (1) clone the github repository 
```
git clone https://github.com/guanjue/snapshot.git
```
#### (2) set conda environment named as "snapshot". 
##### Details about how to install conda can be found in (https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)
```
>>> cat 00_INSTALL_conda_ENV.sh
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels mvdbeek

conda create -n snapshot python=3 r=3.6 bedtools ucsc_tools numpy scikit-learn r-ggplot2 r-pheatmap r-igraph r-networkD3
conda activate snapshot
```

## Input data
##### The cell type peak & signal file list: 1st column is the cell type label; 2nd column is the cell-type specific peak bed file path; 3rd column is the cell-type specific signal bigWig file path
```
>>> head input_data/peak_signal_list.txt 
LSK	atac_pk/LSK.pk.bed	atac_sig/LSK.atac.sig.bw
CMP	atac_pk/CMP.pk.bed	atac_sig/CMP.atac.sig.bw
MEP	atac_pk/MEP.pk.bed	atac_sig/MEP.atac.sig.bw
GMP	atac_pk/GMP.pk.bed	atac_sig/GMP.atac.sig.bw
```

##### The cell type functional state file list: 1st column is the cell type label; 2nd column is the epigenetic state bedgraph file path
```
>>> head input_data/function_list.txt 
LSK	function_label/LSK.ideas.bedgraph
CMP	function_label/CMP.ideas.bedgraph
MEP	function_label/MEP.ideas.bedgraph
GMP	function_label/GMP.ideas.bedgraph
```

##### The cell type differentiation tree: Each row represent one edge in the ell type differentiation tree. The 1st cell type is the progenitor cell type and the 2nd cell type is the differentiated cell type
```
>>> head input_data/cd_tree.txt 
LSK,CMP
CMP,MEP
CMP,GMP
```

##### The functional state color list: 1st column is the epigenetic state label; 2nd column is the epigenetic state RGB color
```
>>> head input_data/function_color_list.txt 
0	255,255,255
1	180,180,180
2	25,160,25
3	126,126,240
4	253,253,157
5	240,185,254
6	253,213,154
7	0,0,212
8	0,146,0
9	250,248,0
```

## RUN Snapshot
##### (1) User need to change the script_folder, input_folder, output_folder, in 'run_Snapshot.sh' file. 
##### The "min_number_per_index_set" is the only parameter user need to decide for Snapshot. It is minimum number of peak per index-set. The index-set with lower number of peaks will be merged into the last X_X_X_... index-set 
```
>>> cat run_Snapshot.sh
##################################
script_folder='/Users/universe/Documents/2022_Independent/snapshot/bin/'
input_folder='/Users/universe/Documents/2022_Independent/00_Independent_analysis/snapshot_test_data/input_data/'
output_folder='/Users/universe/Documents/2022_Independent/00_Independent_analysis/snapshot_test_data/output_result/'
master_peak_bed='/Users/universe/Documents/2022_Independent/00_Independent_analysis/snapshot_test_data/input_data/atac_pk/cCRE.Pool.Merged.bed'

peak_signal_list_file='peak_signal_list.txt'
IDEAS_state_200bp_bed_files_list_file='function_list.txt'
IDEAS_state_color_list_file='function_color_list.txt'
cell_type_tree_file='cd_tree.txt'

output_name='snapshot_test_run'
min_number_per_index_set=10


### run snapshot (CORE!!!)
echo 'run snapshot :o'
cd $input_folder
time python $script_folder'snapshot_v1.py' -p $peak_signal_list_file \
-n $output_name -t $min_number_per_indexset \
-f $IDEAS_state_200bp_bed_files_list_file \
-c $IDEAS_state_color_list_file \
-e $cell_type_tree_file \
-i $input_folder -o $output_folder -s $script_folder \
-m $master_peak_bed
echo 'complete :)'
```
##### (2) use 'runall_commandline.sh' script to run Snapshot
```
time bash run_Snapshot.sh
```


## Output results for test data
### All output files will be to the 'output_folder'

## The heatmap for index set
##### Average atac-seq signal heatmap (left). Most abundant functional state heatmap (right).
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/snapshot_test_run.meansig.png" width="350"/> <img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/snapshot_test_run.indexset_fun.png" width="350"/> 

##### Functional state epigenetic patterns.
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/functional_state_epigenetic_pattern.png" width="350"/>

## The cell differentiation tree for index set 6
##### Average signal tree (left). Most abundant functional state tree (right).
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/6.peak_signal_list.txt1_0_0_1.tree.png" width="400"/> <img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/6.function_list.txt1_0_0_1.tree.png" width="400"/> 

##### Cell type differentiation mean signal violin plot & functional state bar plot
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/6.1_0_0_1.violin.png" width="400"/> <img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/6.1_0_0_1.bar.png" width="400"/> 


##### Merged peak file (bed format)
```
>>> head output_result/snapshot_test_run.sort.bed 
chr1	3445639	3446478	1
chr1	3531951	3532124	2
chr1	3670451	3671268	3
chr1	3672091	3672710	4
chr1	3915538	3915756	5
chr1	4247201	4247354	6
chr1	4332543	4332767	7
chr1	4351941	4352297	8
chr1	4405847	4406057	9
chr1	4412515	4412820	10
```

##### Index set mean signal matrix (bed format)
```
>>> head output_result/snapshot_test_run.meansig.txt 
0_0_0_1	1.1942193919753075	1.2288909759259259	0.5143000398765433	2.2558045432098757
0_0_1_0	1.1964813092783506	1.4681236288659794	4.559413402061856	0.8977813195876289
0_1_0_0	2.3508093636363636	5.0732800000000005	2.2966802727272726	2.4899857272727273
0_1_0_1	1.2273509024999998	3.1057296249999986	0.5968926647500001	4.654063825
0_1_1_0	1.4231207	4.145504333333332	6.617822999999999	1.2878787333333332
1_0_0_0	3.55514962962963	1.6250158148148148	1.213055225925926	1.6773348888888886
1_0_0_1	3.7572421052631584	1.427374894736842	0.6267781684210526	7.191793684210525
1_0_1_0	2.780161818181819	2.0187935000000006	2.9209631363636364	1.9096209545454548
1_1_0_0	5.526250833333332	3.49659475	0.6552302083333333	1.249175055555556
1_1_0_1	6.688364330508474	7.341211720338984	0.9215638576271191	8.187531830508476
```

##### Index signal matrix (bed format)
```
>>> head head output_result/snapshot_test_run.sig.txt 
158	0_0_0_1	0.461536	0.97629	0.420611	1.17384
122	0_0_0_1	1.05508	1.59315	0.332113	3.5154
124	0_0_0_1	0.570329	1.321	0.197305	5.40436
590	0_0_0_1	0.495322	0.82881	0.202169	1.4323
126	0_0_0_1	1.05409	0.395961	0.25931	2.14795
589	0_0_0_1	2.23673	1.88691	0.730118	3.32344
587	0_0_0_1	1.57064	0.631929	0.13965	0.517572
129	0_0_0_1	0.5076	0.613323	0.229215	1.54835
130	0_0_0_1	1.38835	2.22171	0.0815782	3.49794
121	0_0_0_1	2.41654	3.42985	0.510925	5.30627
```

##### Index set most abundant functional state matrix (bed format)
```
>>> head output_result/snapshot_test_run.indexset_fun.txt 
0_0_0_1	0	0	0	20
0_0_1_0	0	0	20	0
0_1_0_0	0	20	11	20
0_1_0_1	0	20	0	12
0_1_1_0	0	20	12	1
1_0_0_0	20	0	0	0
1_0_0_1	20	4	0	12
1_0_1_0	20	0	0	0
1_1_0_0	12	20	0	0
1_1_0_1	12	12	0	12
```

##### Index functional state matrix (bed format)
```
>>> head output_result/snapshot_test_run.fun.txt         
158	0_0_0_1	7	4	0	0
122	0_0_0_1	0	11	7	11
124	0_0_0_1	4	4	0	20
590	0_0_0_1	0	0	0	0
126	0_0_0_1	0	0	0	4
589	0_0_0_1	25	13	0	20
587	0_0_0_1	7	7	7	7
129	0_0_0_1	10	0	0	0
130	0_0_0_1	4	20	0	20
121	0_0_0_1	20	20	0	15
```






