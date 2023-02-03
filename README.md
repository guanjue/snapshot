# Snapshot
### Epigenetic modification of chromatin plays a pivotal role in regulating gene expression during cell differentiation. The scale and complexity of epigenetic data pose significant challenges for biologists to identify the regulatory events controlling cell differentiation. Here, we present a new method, called Snapshot, that uses epigenetic data to generate a hierarchical visualization for DNA regions with epigenetic features segregating along any given cell differentiation hierarchy of interest. Different hierarchies of cell types may be used to highlight the epigenetic history specific to any particular cell lineage. We demonstrate the utility of Snapshot using data from the VISION project, an international project for ValIdated Systematic IntegratiON of epigenomic data in mouse and human hematopoiesis.

## Reference
#### Xiang, G., Giardine, B., An, L., Sun, C., Keller, C., Heuston, E., Bodine, D., Hardison, R. and Zhang, Y. Snapshot: clustering and visualizing epigenetic history during cell differentiation. bioRxiv, p.291880.


## The overall workflow of Snapshot

#### Figure 1
![logo](https://raw.githubusercontent.com/guanjue/snapshot/master/test_data/example/f1.png)

##### The overall workflow of snapshot. (a) In the 1st step, we use the binarized presence/absence status of ccREs at each location across all cell types to create a ccRE index to represent the unique combinatorial pattern of ccREs at the location. All ccREs with the same index are grouped into one cluster called an index-set (IS). (b) We select the threshold to filter the less prevalent ISs based on the distribution of the number of ccREs per IS (red dash line). (c) We visualize the ISs and their ccREs in a heatmap. The blue dash box represents one example of IS: 0_1_0_0_1. The filtered ISs and their ccREs are shown in the end of the heatmap (black dash box). (d) We apply the QDA model to correct ccRE’s index by borrowing information from the signals of the ccRE across all cell types. After the QDA step, some of the ccREs grouped in the rare ISs can be rescued. (e) Finally, we visualize the ccRE clusters in a heatmap. Each row in the heatmap is the ccRE pattern for each index-set, and each column is a cell type. The ccRE patterns are sorted by their indices in the heatmap. By our definition of the ccRE index, the ccRE patterns are separated if they have different ccRE status in a cell type; conversely they are clustered together if they have similar ccRE status in a cell type. These different plots will together highlight the epigenetic activity across cell types and the associated functional annotations and their enrichments for each index-set, which enhances the interpretation of the functional roles of each ccRE cluster during cell differentiation.

## Example of Snapshot output: 
### Hematopoietic cell differentiation in VISION (ValIdated Systematic IntegratiON of hematopoietic epigenomes) project
#### Figure 2
![logo](https://raw.githubusercontent.com/guanjue/snapshot/master/test_data/example/f2.png)

##### The heatmap of index-sets. (a) The heatmap of index-set colored by the average ATAC-seq signal in each cell type. (b) The heatmap of index-set colored by the most frequent functional annotation in each cell type. (c) The density plot of the number of genomic region covered by the index-set. (d) The color code and epigenetic composition of functional annotation used in (b).

#### Figure 3
![logo](https://raw.githubusercontent.com/guanjue/snapshot/master/test_data/example/f3.png)

##### The data visualization for index-set-149 and corresponding GO analysis and MEME-ChIP TF binding motif analysis. (a) The hematopoietic cell differentiation tree colored by the average ATAC-seq signal in each cell type of the index-set-149. The violin plot represents the distribution of ATAC-seq signal in each cell type of the index-set-149 is in below.  (b) The same cell differentiation tree colored by the most frequent functional annotation in each cell type of the index-set-149. The two most frequence functional annotation in erythroblasts lineage. The bar plot based on the proportion of each functional annotation in each cell type of the index-set-149 is below the cell differentiation tree. (c) The index-set-149 relevant GO term. (d) The index-set-149 significantly enriched TF binding motif from MEME-ChIP analysis.



# Snapshot Installation Guide

## Dependencies
The following packages and tools are required to run Snapshot:
- Python:
  - numpy
  - sklearn
- R:
  - ggplot2
  - pheatmap
  - igraph
  - networkD3
- bedtools

## Installation Steps
1. Clone the Github repository:
```
git clone https://github.com/guanjue/snapshot.git
```
2. Set up a conda environment named "snapshot":
```
conda create -n snapshot python=3 r=3.6 bedtools ucsc_tools numpy scikit-learn r-ggplot2 r-pheatmap r-igraph r-networkD3 r-data.table r-mclust r-dplyr r-lsa r-cba r-RColorBrewer r-tidyverse
conda activate snapshot
```

Note: Detailed instructions on how to install conda can be found in the [conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html).


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

conda create -n snapshot_test2 python=3 r=3.6 bedtools ucsc_tools numpy scikit-learn r-ggplot2 r-pheatmap r-igraph r-networkD3 r-data.table r-mclust r-dplyr r-lsa r-cba r-RColorBrewer r-tidyverse 
conda activate snapshot_test2
```

## Input data
###### The cell type peak & signal & epigenetic-state file list: 
###### The 1st column is the cell type label; 
###### The 2nd column is the cell-type specific peak bed absolute-file-path; 
###### The 3rd column is the cell-type specific signal bigWig absolute-file-path; 
###### The 4th column is the cell-type specific functional epigenetic state bigBed absolute-file-path.
```
>>> head input_data/peak_signal_state_list.txt 
LSK	atac_pk/LSK.pk.bed	atac_sig/ideasVisionV20p8NormLskAtac.bw	function_label/ideasVisionV20p8SegLsk.bb
CMP	atac_pk/CMP.pk.bed	atac_sig/ideasVisionV20p8NormCmpAtac.bw	function_label/ideasVisionV20p8SegCmp.bb
MEP	atac_pk/MEP.pk.bed	atac_sig/ideasVisionV20p8NormMepAtac.bw	function_label/ideasVisionV20p8SegMep.bb
CFUE	atac_pk/CFUE.pk.bed	atac_sig/ideasVisionV20p8NormCfueAtac.bw	function_label/ideasVisionV20p8SegCfue.bb
ERYad	atac_pk/ERYad.pk.bed	atac_sig/ideasVisionV20p8NormEryadAtac.bw	function_label/ideasVisionV20p8SegEryad.bb
ERYfl	atac_pk/ERYfl.pk.bed	atac_sig/ideasVisionV20p8NormEryflAtac.bw	function_label/ideasVisionV20p8SegEryfl.bb
G1E	atac_pk/G1E.pk.bed	atac_sig/ideasVisionV20p8NormG1eAtac.bw	function_label/ideasVisionV20p8SegG1e.bb
ER4	atac_pk/ER4.pk.bed	atac_sig/ideasVisionV20p8NormEr4Atac.bw	function_label/ideasVisionV20p8SegEr4.bb
CFUMK	atac_pk/CFUMK.pk.bed	atac_sig/ideasVisionV20p8NormCfumAtac.bw	function_label/ideasVisionV20p8SegCfum.bb
iMK	atac_pk/iMK.pk.bed	atac_sig/ideasVisionV20p8NormImkAtac.bw	function_label/ideasVisionV20p8SegImk.bb
.......
```

##### The cell type differentiation tree: Each row represent one edge in the ell type differentiation tree. The 1st cell type is the progenitor cell type and the 2nd cell type is the differentiated cell type
```
>>> head input_data/cd_tree.txt 
LSK,CMP
CMP,MEP
MEP,CFUE
CFUE,ERYad
ERYad,ERYfl
MEP,CFUMK
CFUMK,iMK
MEP,G1E
G1E,ER4
CMP,GMP
......
```

##### The functional state color list: 1st column is the epigenetic state label; 2nd column is the epigenetic state RGB color
```
>>> head input_data/function_color_list.txt 
0	255,255,255
1	0,145,0
2	144,144,144
3	49,49,231
4	251,251,75
5	250,230,4
6	253,222,175
7	231,60,251
8	0,148,0
9	207,77,155
10	245,0,0
......
```

## RUN Snapshot
##### (1) User need to change the script_folder, input_folder, output_folder, in 'Step1_run_Snapshot.sh' file. 
##### The file path or folder path should be replace by the real absolute-path to the target folder

##### This script is a pipeline that runs the "snapshot_v2.py" program, which is used for analyzing genomic data. The script specifies various required and optional parameters and input files, which are passed as arguments to the "snapshot_v2.py" program. The input files are specified as absolute file paths and include:

```
peak_signal_list_file - A file that lists the peak signals to be analyzed.
IDEAS_state_color_list_file - A file that provides the color scheme for the functional states of genomic regions.
cell_type_tree_file - A file that represents the hierarchical structure of the cell types being analyzed.
genome_size_file - A file that lists the sizes of the chromosomes in the genome being analyzed.
input_folder - The path to the folder that contains the input data.
output_folder - The path to the folder where the output of the analysis will be stored.
script_folder - The path to the folder that contains the "snapshot_v2.py" script.
```
In addition to these required parameters, the script also specifies several optional parameters, including:
```
master_peak_bed - A file that provides information about the peaks to be analyzed.
min_number_per_indexset - The minimum number of cCRES points required for an IndexSet to include.
QDA_round_num - The number of rounds of the "Quantitative Domain Analysis" (QDA) algorithm to be run.
normalization_method - The method to be used for normalizing the data.
have_function_state_files - A flag indicating whether the input data includes functional state information.
index_matrix_txt - A file that provides information about the data points to be analyzed.
signal_matrix_txt - A file that provides the signal data to be analyzed.
function_state_matrix_txt - A file that provides information about the functional states of the genomic regions to be analyzed.
```
Finally, the script changes to the "input_folder" directory and runs the "snapshot_v2.py" program with the specified parameters and input files. The output of the analysis will be stored in the "output_folder". After the program completes, the script returns to the previous directory and outputs "complete :)".


```
>>> cat Step1_run_Snapshot.sh
##################################
### required parameters or input files
output_name='snapshot_test_run_merge'
peak_signal_list_file='peak_signal_state_list.merge.txt'
IDEAS_state_color_list_file='function_color_list.txt'
cell_type_tree_file='cd_tree.txt'
genome_size_file='/Users/universe/Downloads/Snapshot_test/input_data_hg38/hg38.chrom.1_22XY.sizes'

### required folder path
input_folder='/Users/universe/Downloads/Snapshot_test/input_data_hg38/'
output_folder='/Users/universe/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_merge/'
script_folder='/Users/universe/Documents/projects/snapshot/bin/'

### optional parameters or input files
master_peak_bed='/Users/universe/Downloads/Snapshot_test/input_data_hg38/snapshot_test_run_merge.bedinfo.bed'
min_number_per_indexset=100
QDA_round_num=1
normalization_method=S3norm
have_function_state_files=F

### input matrix
index_matrix_txt='/Users/universe/Downloads/Snapshot_test/input_data_hg38/used_input_matrix/snapshot_test_run_merge.index.matrix.txt'
signal_matrix_txt='/Users/universe/Downloads/Snapshot_test/input_data_hg38/used_input_matrix/snapshot_test_run_merge.signal.matrix.txt'
function_state_matrix_txt='/Users/universe/Downloads/Snapshot_test/input_data_hg38/used_input_matrix/snapshot_test_run_merge.function.matrix.txt'


### run snapshot (CORE!!!)
echo 'run snapshot :o'
cd $input_folder


time python3 $script_folder/snapshot_v2.py -p $peak_signal_list_file \
-n $output_name -t $min_number_per_indexset \
-f $genome_size_file \
-c $IDEAS_state_color_list_file \
-e $cell_type_tree_file \
-i $input_folder -o $output_folder -s $script_folder \
-m $master_peak_bed -q $QDA_round_num -z $normalization_method \
-b $have_function_state_files \
-a $index_matrix_txt \
-r $signal_matrix_txt \
-g $function_state_matrix_txt


cd ..
echo 'complete :)'

```

##### (2) use 'Step0_get_input_data_for_testing_run.sh' script to download input data for testing run
```
cd /snapshot/test_data/
time bash Step0_get_input_data_for_testing_run.sh
```

##### (3) use 'runall_commandline.sh' script to run Snapshot. The testing run can be finish within 20min (Macbook pro, 2.4 GHz Quad-Core Intel Core i5; 8GB RAM)
```
time bash Step1_run_Snapshot.sh
```


## Output results for test data
### All output files will be to the 'output_folder'

## The heatmap for index set
##### Average ATAC-seq signals patterns in all Index-Sets.
![logo](https://raw.githubusercontent.com/guanjue/snapshot/master/test_data/example/snapshot_test_run.meansig.png)

##### Most Abundant Functional Epigenetic States patterns in all Index-Sets.
![logo](https://raw.githubusercontent.com/guanjue/snapshot/master/test_data/example/snapshot_test_run.indexset_fun.png)

##### Cell type differentiation system (Index-Set 31): Cell differentiation Tree of Average signals
![logo](https://raw.githubusercontent.com/guanjue/snapshot/master/test_data/example/31.peak_signal_state_list.tree.signal.png)

##### Cell type differentiation system (Index-Set 31): Cell differentiation Tree of Functional Epigenetic States
![logo](https://raw.githubusercontent.com/guanjue/snapshot/main/test_data/example/31.peak_signal_state_list.tree.state.png)

##### Cell type differentiation system (Index-Set 31): Violin of peak signals
![logo](https://raw.githubusercontent.com/guanjue/snapshot/master/test_data/example/31.violin.png)

##### Cell type differentiation system (Index-Set 31): Barplot of Functional Epigenetic States
![logo](https://raw.githubusercontent.com/guanjue/snapshot/master/test_data/example/31.bar.png)


##### Merged peak file (bed format)
```
>>> head output_result/snapshot_test_run.sort.bed 
chr1	3446000	3446200	1
chr1	3451600	3452000	2
chr1	3914600	3915400	3
chr1	4414800	4415000	4
chr1	4416800	4417000	5
chr1	4418000	4418200	6
chr1	4424800	4425000	7
chr1	4490600	4490800	8
chr1	4496400	4496600	9
chr1	4571800	4572000	10
......
```

##### Index set mean signal matrix (bed format)
```
>>> head output_result/snapshot_test_run.meansig.txt 
0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	2.7203165450529787	3.1321496870415695	3.447007176430317	3.166636248981257	3.552099701898949	3.878185635696823	1.4819312050529738	1.7342681281988508	3.3944306135289284	2.7201608405053013	3.0900341948655226	3.2680497713121426	2.2836158302363523	1.2009829315403382	2.6714230573757036	0.8960898602282016	0.8786743454767777
0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_1_1	1.621244671099654	1.6409139173539502	1.3145519881786947	0.8115011299656362	0.8404619017697598	0.5566642946735354	0.48029953384879853	0.6203494304123708	1.3355942798969065	1.6299852580756033	1.3657937759450174	1.2376911058419244	1.923485744158075	1.998705881443296	1.732353492955328	11.136037485223389	11.323727326804129
0_0_0_0_0_0_0_0_0_0_0_0_0_1_0_0_0	0.7655575370357036	0.8487737590320805	0.7240723858022319	0.544760997469177	0.587946304191074	0.44096724214616756	0.3034562685433762	0.4000610984323545	0.7191299221277548	1.0813688480200843	0.8120614006359809	1.3240255077631793	0.6157234188708521	6.308377507913562	0.5194036706833899	0.40839617461087263	0.5462352972925977
0_0_0_0_0_0_0_0_0_0_0_0_0_1_0_0_1	1.630883616570327	1.62141689479769	1.1632295371868975	0.7996866307899804	0.8271003921387284	0.5240952834296697	0.5384812385356467	0.5470913254335256	1.291946048747591	3.2170100712909417	1.5444558283236993	3.318453157996146	1.932036053872834	13.108413119460499	1.9525637666666669	1.9164648273603078	7.640808234296725
0_0_0_0_0_0_0_0_0_0_0_0_0_1_0_1_1	2.3792834576888886	2.1911267263999954	1.3348013142577793	0.7720584791822235	0.8296216599466655	0.6163781108444436	0.5498947246222227	0.6311941117333376	1.713560748711111	3.155713040177777	1.8039701929777778	2.6633003584337787	2.202068096835559	11.958543297777787	2.9996011282666664	9.450972157155567	9.93109183484443
0_0_0_0_0_0_0_0_0_0_0_0_1_0_0_0_0	1.0011610849844872	0.9538348655118912	0.7408574598138579	0.6121845035780783	0.6614750357445695	0.449885089710446	0.3671608014994848	0.4407027254395075	0.8831982698552219	1.0661161927094087	1.1324645551189223	1.6663684419017564	7.369715429989662	0.49081203412615776	13.063455216804595	0.4003475963288561	0.29709514389866326
0_0_0_0_0_0_0_0_0_0_0_0_1_0_0_1_1	1.1236453496782168	0.9294576571782167	0.7909972775495048	0.47438374594059424	0.523503076039604	0.39734370433167876	0.4203407258663374	0.42161150358910876	0.8437542070544551	1.3763707129950502	0.9310128323019788	1.2495086113564355	5.532267155445545	0.9878571967821775	6.183162135148522	10.555262302970316	10.365179264108914
0_0_0_0_0_0_0_0_0_0_0_0_1_1_0_1_1	2.6014175861685227	2.2045709594594545	1.2658092203815567	0.7042049232114471	0.8088810879491259	0.4923063357710608	0.41330659220985755	0.505679246422893	1.6570122882352911	3.2779114435612033	1.8122156941176475	4.138298289189187	9.942824578696339	9.628220720190775	10.855374585373601	9.840905967567567	10.101142152305247
0_0_0_0_0_0_0_0_0_0_0_1_0_0_0_0_0	1.3830720146739126	1.4681040086956503	1.1355101451630438	1.0165573032608692	1.1735018167119562	0.846829470108693	0.5238043187500003	0.6801225817934781	1.4602143374999998	3.075292687499999	2.381737922826086	14.848748043478256	1.8427717603260858	0.7328498288043491	0.9578335339673909	0.4317764978260879	0.3456542766304343
0_0_0_0_0_0_0_0_0_0_0_1_0_1_0_0_0	1.2912370071770327	1.4225890755980857	1.1894378943540667	0.8568112155980858	0.9286863127751195	0.5038397952153117	0.7489491789473682	0.7766797397129189	1.222042716267943	3.3206504306220097	1.953926073205741	8.564044961722491	1.2924173251674644	13.53591655502393	0.8503540727272724	0.4893061593301443	0.6661528889952149
......
```

##### Index signal matrix (bed format)
```
>>> head head output_result/snapshot_test_run.sig.txt 
43168	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	7.07691	6.18274	9.12767	3.82776	5.58093	5.89785	2.03239	3.37089	8.27097	6.85553	8.26579	2.67891	3.26334	4.66176	5.90397	2.9205	3.08933
95287	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	0.637461	0.236972	0.690395	0.993535	0.785778	2.39158	2.31543	0.391602	0.681384	0.343735	0.597157	0.182025	1.70661	0.101864	0.0801293	0.0283914	0.0280334
7680	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	1.3759	1.90162	2.36344	1.33513	0.416245	0.438437	0.369149	0.188444	3.2301	2.1975	2.29221	0.999266	0.0330778	1.43727	0.542269	0.255617	0.0280334
39986	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	4.00286	8.25309	8.1995	1.16145	2.1636	6.34651	3.62213	4.13962	6.88729	3.99581	6.16218	3.3524	2.71671	1.17237	1.55326	2.03524	1.9255
48720	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	3.2115	2.02981	1.96464	1.25283	1.20164	4.01863	2.49397	2.96509	0.259521	0.396701	0.410677	1.48194	1.08404	0.300992	0.63197	0.0283914	0.187807
14028	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	2.79874	3.03816	3.44399	5.77083	4.17465	5.16375	1.57511	2.43745	6.25566	4.17458	5.19544	6.41478	4.9301	1.34674	2.9481	0.829629	0.444549
14029	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	0.619248	2.72965	1.57756	1.10768	1.2897	4.97009	0.0971869	0.778748	1.99531	1.70862	1.80994	1.20135	4.12433	0.480249	0.179645	0.0283914	0.0280334
91357	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	7.65608	9.42799	9.92788	5.64864	5.32062	5.35844	4.48675	6.38107	8.89295	3.55431	5.25201	3.95081	4.8496	0.391709	0.902306	1.92658	0.444549
63083	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	4.41733	3.28725	2.09059	1.40121	1.55031	1.55343	0.140052	0.268098	4.16032	2.93077	4.49774	2.60125	5.01934	1.45313	10.9039	2.61408	3.26873
2294	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	4.81572	6.29913	2.71223	2.88013	2.54477	5.55415	3.33966	2.61847	5.79641	4.22963	4.85441	7.19835	1.77697	0.101864	1.1775	0.245224	0.187807
......
```

##### Index set most abundant functional state matrix (bed format)
```
>>> head output_result/snapshot_test_run.indexset_fun.txt 
0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_1_1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	9	9
0_0_0_0_0_0_0_0_0_0_0_0_0_1_0_0_0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0_0_0_0_0_0_0_0_0_0_0_0_0_1_0_0_1	0	0	0	0	0	0	0	0	0	0	0	0	0	9	0	0	9
0_0_0_0_0_0_0_0_0_0_0_0_0_1_0_1_1	0	0	0	0	0	0	0	0	0	0	0	0	0	9	0	9	9
0_0_0_0_0_0_0_0_0_0_0_0_1_0_0_0_0	0	0	0	0	0	0	0	0	0	0	0	0	9	0	9	0	0
0_0_0_0_0_0_0_0_0_0_0_0_1_0_0_1_1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	9	9
0_0_0_0_0_0_0_0_0_0_0_0_1_1_0_1_1	0	0	0	0	0	0	0	0	0	0	0	9	9	9	9	9	9
0_0_0_0_0_0_0_0_0_0_0_1_0_0_0_0_0	0	0	0	0	0	0	0	0	0	0	0	9	0	0	0	0	0
0_0_0_0_0_0_0_0_0_0_0_1_0_1_0_0_0	0	0	0	0	0	0	0	0	0	0	0	9	0	9	0	0	0
......
```

##### Index functional state matrix (bed format)
```
>>> head output_result/snapshot_test_run.fun.txt         
43168	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	10	10	10	11	11	21	10	10	15	10	10	10	15	15	15	18	15
95287	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	0	0	0	0	0	0	7	0	0	0	0	0	0	0	0	0	0
7680	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	11	0	15	0	0	0	0	0	15	0	11	18	18	0	0	0	0
39986	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	13	13	13	7	7	13	26	26	13	13	13	13	7	7	7	7	7
48720	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	0	0	0	0	0	0	4	0	0	0	0	7	0	0	0	0	0
14028	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	19	18	11	9	9	21	19	19	21	21	22	22	21	18	19	0	18
14029	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	7	7	7	7	7	13	7	7	7	7	7	7	13	7	7	7	7
91357	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	13	13	13	13	13	13	13	13	13	13	13	13	13	7	7	7	13
63083	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	13	13	13	13	13	13	7	7	13	13	13	13	13	13	26	13	13
2294	0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0	13	13	7	7	7	13	7	7	13	13	13	13	7	7	7	7	7
......
```






