##################################
script_folder='/Users/universe/Documents/2018_BG/snapshot/bin/'
input_folder='/Users/universe/Documents/2018_BG/snapshot/test_data/input_data/'
output_folder='/Users/universe/Documents/2018_BG/snapshot/test_data/output_result/'


### run snapshot (CORE!!!)
echo 'run snapshot :o'
cd $input_folder
time python $script_folder'snapshot_v.0.4.py' -p peak_list.txt -n snapshot -t 1 -s signal_list.txt -l F -z F -x 0.01 -f function_list.txt -m mostfreq -c function_color_list.txt -e cd_tree.txt -i $input_folder -o $output_folder -b $script_folder -q 1
echo 'complete :)'
