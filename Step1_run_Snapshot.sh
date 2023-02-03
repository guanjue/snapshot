### read parameter setting info file
input_info_file=$1
param(){
	cat $input_info_file | grep -w $1 | cut -f 2
}

##################################
### required parameters or input files
output_name=$(param "output_name")
peak_signal_list_file=$(param "peak_signal_list_file")
IDEAS_state_color_list_file=$(param "IDEAS_state_color_list_file")
cell_type_tree_file=$(param "cell_type_tree_file")
genome_size_file=$(param "genome_size_file")

### required folder path
input_folder=$(param "input_folder")
output_folder=$(param "output_folder")
script_folder=$(param "script_folder")

### optional parameters or input files
master_peak_bed=$(param "master_peak_bed")
min_number_per_indexset=$(param "min_number_per_indexset")
normalization_method=$(param "normalization_method")
have_function_state_files=$(param "have_function_state_files")

### Use matrices as Input
index_matrix_txt=$(param "index_matrix_txt")
signal_matrix_txt=$(param "signal_matrix_txt")
function_state_matrix_txt=$(param "function_state_matrix_txt")

echo 'Required parameters or input files'
echo $output_name
echo $peak_signal_list_file
echo $IDEAS_state_color_list_file
echo $cell_type_tree_file
echo $genome_size_file
echo 'Required folder path'
echo $input_folder
echo $output_folder
echo $script_folder
echo 'optional parameters or input files'
echo $master_peak_bed
echo $min_number_per_indexset
echo $normalization_method
echo $have_function_state_files
echo 'Use matrices as Input settings'
echo $index_matrix_txt
echo $signal_matrix_txt
echo $function_state_matrix_txt


### run snapshot (CORE!!!)
echo 'run snapshot :o'
cd $input_folder

time python3 $script_folder/snapshot_v2.py -p $peak_signal_list_file \
-n $output_name -t $min_number_per_indexset \
-f $genome_size_file \
-c $IDEAS_state_color_list_file \
-e $cell_type_tree_file \
-i $input_folder -o $output_folder -s $script_folder \
-m $master_peak_bed -z $normalization_method \
-b $have_function_state_files \
-a $index_matrix_txt \
-r $signal_matrix_txt \
-g $function_state_matrix_txt


echo 'complete :)'
