The details about the input files can be found at following links:
- [Input Matrix format](https://github.com/guanjue/snapshot/blob/main/INPUT_format_matrix.md).
- [Input format for Raw data](https://github.com/guanjue/snapshot/blob/main/INPUT_format_raw.md).


```
##################################
### required parameters or input files
output_name	snapshot_test_run_merge
peak_signal_list_file: peak_signal_state_list.merge.txt
IDEAS_state_color_list_file: function_color_list.txt
cell_type_tree_file: cd_tree.txt
genome_size_file: chromosome size. (e.g. /Users/universe/Downloads/Snapshot_test/input_data_hg38/hg38.chrom.1_22XY.sizes)

### required folder path
input_folder: working directory. (e.g. /Users/universe/Downloads/Snapshot_test/input_data_hg38/)
output_folder: output directory. (e.g. /Users/universe/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_merge/)
script_folder: Snapshot script directory. (/Users/universe/Documents/projects/snapshot/bin/)

### optional parameters or input files
master_peak_bed: user provided master peak list. (F/bedfile.bed)
min_number_per_indexset: the minimum number of cCREs required for a abundant Index-Set. (numeric value e.g. 100)
normalization_method: different internal normalization methods. (F/S3norm/scale_quantile/QTnorm)
have_function_state_files: check if the user have functional epigenetic state data for the analysis. (T/F)

### input matrix
index_matrix_txt: user provided binary index matrix for each cCRE in each cell-type (F/snapshot_test_run_merge.index.matrix.txt)
signal_matrix_txt: user provided signal matrix for each cCRE in each cell-type (F/snapshot_test_run_merge.signal.matrix.txt)
function_state_matrix_txt: user provided functional epigenetic state matrix for each cCRE in each cell-type (F/snapshot_test_run_merge.function.matrix.txt)


# ‘F’: means user do NOT have the data / the procedure is NOT be included in the analysis
```