##################################
script_folder='/Users/universe/Documents/2022_Independent/snapshot/bin/'
input_folder='/Users/universe/Documents/2022_Independent/00_Independent_analysis/snapshot_test_data/input_data/'
output_folder='/Users/universe/Documents/2022_Independent/00_Independent_analysis/snapshot_test_data/output_result/'
master_peak_bed='/Users/universe/Documents/2022_Independent/00_Independent_analysis/snapshot_test_data/input_data/S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.bed'

peak_signal_list_file='peak_signal_list.txt'
IDEAS_state_200bp_bed_files_list_file='function_list.txt'
IDEAS_state_color_list_file='function_color_list.txt'
cell_type_tree_file='cd_tree.txt'

output_name='snapshot_test_run'
min_number_per_indexset=10


### run snapshot (CORE!!!)
echo 'run snapshot :o'
cd $input_folder

time python $script_folder'snapshot_v1.py' -p $peak_signal_list_file \
-n $output_name -t $min_number_per_indexset \
-f $IDEAS_state_200bp_bed_files_list_file \
-c $IDEAS_state_color_list_file \
-e $cell_type_tree_file \
-i $input_folder -o $output_folder -s $script_folder \
#-m $master_peak_bed

echo 'complete :)'
