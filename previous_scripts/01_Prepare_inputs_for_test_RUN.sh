##########################################
###### Generate input_data folders ######
cd input_data/
mkdir atac_pk
mkdir atac_sig
mkdir function_label
##########################################


##########################################
###### Download cCRE bed file (Master peak list) ######
wget https://usevision.org/data/hg38/IDEASstates/ideasJointMay2021/S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.bed
##########################################


##########################################
###### Download cell-type peak bed file ######
wget https://usevision.org/data/mm10/atac/842.homerPeaks.noBlack.noChrM.bed
wget https://usevision.org/data/mm10/atac/843.homerPeaks.noBlack.noChrM.bed
wget https://usevision.org/data/mm10/atac/987.homerPeaks.noBlack.noChrM.bed
wget https://usevision.org/data/mm10/atac/844.homerPeaks.noBlack.noChrM.bed
### move bigwig files into atac_sig folder
mv *.homerPeaks.noBlack.noChrM.bed atac_pk/
##########################################


##########################################
###### Download ATAC-seq signal bigwig file ######
wget https://usevision.org/data/hg38/atac/CMP_100246.bw
wget https://usevision.org/data/hg38/atac/GMP_100256.bw
wget https://usevision.org/data/hg38/atac/HSC_100258.bw
wget https://usevision.org/data/hg38/atac/MEP_Donor2596.bw
wget https://usevision.org/data/hg38/atac/ERY_S002R5.bw
### move bigwig files into atac_sig folder
mv *.bw atac_sig/
##########################################


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



