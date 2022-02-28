##################################
### required parameters or input files
output_name='snapshot_test_run'
peak_signal_list_file='peak_signal_state_list.txt'
IDEAS_state_color_list_file='function_color_list.txt'
cell_type_tree_file='cd_tree.txt'
genome_size_file='mm10.chrom.1_19XY.sizes'

### required folder path
input_folder='/Users/universe/Documents/2022_Independent/00_Independent_analysis/snapshot_test_data/input_data/'
output_folder='/Users/universe/Documents/2022_Independent/00_Independent_analysis/snapshot_test_data/output_result/'
script_folder='/Users/universe/Documents/2022_Independent/snapshot/bin/'

### optional parameters or input files
master_peak_bed='/Users/universe/Documents/2022_Independent/snapshot/test_data/input_data/S3V2_IDEAS_mm10_ccre2.cCRE.M.notall0.bed'
min_number_per_indexset=100
QDA_round_num=1

### run snapshot (CORE!!!)
echo 'run snapshot :o'
cd $input_folder

time python3 $script_folder'snapshot_v1.py' -p $peak_signal_list_file \
-n $output_name -t $min_number_per_indexset \
-f $genome_size_file \
-c $IDEAS_state_color_list_file \
-e $cell_type_tree_file \
-i $input_folder -o $output_folder -s $script_folder \
-m $master_peak_bed -q $QDA_round_num

echo 'complete :)'
