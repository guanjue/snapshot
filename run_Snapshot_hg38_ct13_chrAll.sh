##################################
### required parameters or input files
output_name='snapshot_test_run_merge'
peak_signal_list_file='peak_signal_state_list.merge.txt'
IDEAS_state_color_list_file='function_color_list.txt'
cell_type_tree_file='cd_tree.txt'
genome_size_file='/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38.chrom.1_22XY.sizes'

### required folder path
input_folder='/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/'
output_folder='/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_merge/'
script_folder='/Users/guanjuexiang/Documents/projects/snapshot/bin/'

### optional parameters or input files
master_peak_bed='/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/snapshot_test_run_merge.bedinfo.bed'
min_number_per_indexset=100
QDA_round_num=100
normalization_method=S3norm

### run snapshot (CORE!!!)
echo 'run snapshot :o'
cd $input_folder


time python3 $script_folder/snapshot_v2.py -p $peak_signal_list_file \
-n $output_name -t $min_number_per_indexset \
-f $genome_size_file \
-c $IDEAS_state_color_list_file \
-e $cell_type_tree_file \
-i $input_folder -o $output_folder -s $script_folder \
-m $master_peak_bed -q $QDA_round_num -z $normalization_method

cd ..
echo 'complete :)'


