##########################################
###### Generate input_data folders ######
cd input_data/
mkdir function_label

##########################################
###### Download IDEAS state matrix ######
wget https://usevision.org/data/hg38/IDEASstates/ideasJointMay2021/S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state
input_IDEAS_state_file='S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state'
#
### Extract IDEAS state bed files per cell type from the state matrix 
output_file='function_label/CMP.J_IDEAS.bed'
target_cell_col=12
cat $input_IDEAS_state_file | awk -F ' ' -v OFS='\t' -v used_col=$target_cell_col '{if ($2=="CHR") print "#"$2,$3,$4, $used_col; else print $2,$3,$4, $used_col}' > $output_file
#
output_file='function_label/GMP.J_IDEAS.bed'
target_cell_col=18
cat $input_IDEAS_state_file | awk -F ' ' -v OFS='\t' -v used_col=$target_cell_col '{if ($2=="CHR") print "#"$2,$3,$4, $used_col; else print $2,$3,$4, $used_col}' > $output_file
#
output_file='function_label/HSC.J_IDEAS.bed'
target_cell_col=20
cat $input_IDEAS_state_file | awk -F ' ' -v OFS='\t' -v used_col=$target_cell_col '{if ($2=="CHR") print "#"$2,$3,$4, $used_col; else print $2,$3,$4, $used_col}' > $output_file
#
output_file='function_label/MEP.J_IDEAS.bed'
target_cell_col=30
cat $input_IDEAS_state_file | awk -F ' ' -v OFS='\t' -v used_col=$target_cell_col '{if ($2=="CHR") print "#"$2,$3,$4, $used_col; else print $2,$3,$4, $used_col}' > $output_file
#
output_file='function_label/ERY.J_IDEAS.bed'
target_cell_col=16
cat $input_IDEAS_state_file | awk -F ' ' -v OFS='\t' -v used_col=$target_cell_col '{if ($2=="CHR") print "#"$2,$3,$4, $used_col; else print $2,$3,$4, $used_col}' > $output_file
##########################################



