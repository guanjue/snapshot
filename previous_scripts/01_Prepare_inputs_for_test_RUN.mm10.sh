##########################################
###### Generate input_data folders ######
cd input_data/
mkdir atac_pk
mkdir atac_sig
mkdir function_label
##########################################


##########################################
###### Download cCRE bed file (Master peak list) ######
wget wget https://usevision.org/data/hg38/IDEASstates/ideasJointMay2021/S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state
https://usevision.org/data/mm10/ideasJointMay2021/S3V2_IDEAS_mm10_ccre2.cCRE.M.notall0.bed
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
wget https://usevision.org/data/mm10/IDEAShumanHem/bws_RC/CMP_842.ATAC.S3V2.bedgraph.bw
wget https://usevision.org/data/mm10/IDEAShumanHem/bws_RC/GMP_843.ATAC.S3V2.bedgraph.bw
wget https://usevision.org/data/mm10/IDEAShumanHem/bws_RC/LSK_987.ATAC.S3V2.bedgraph.bw
wget https://usevision.org/data/mm10/IDEAShumanHem/bws_RC/MEP_844.ATAC.S3V2.bedgraph.bw
### move bigwig files into atac_sig folder
mv *.ATAC.S3V2.bedgraph.bw atac_sig/
##########################################


##########################################
###### Download IDEAS state matrix ######
wget https://usevision.org/data/mm10/ideasJointMay2021/S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state
input_IDEAS_state_file='S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state'
#
### Extract IDEAS state bed files per cell type from the state matrix 
output_file='function_label/CMP.J_IDEAS.bed'
target_cell_col=11
cat $input_IDEAS_state_file | awk -F ' ' -v OFS='\t' -v used_col=$target_cell_col '{if ($2=="CHR") print "#"$2,$3,$4, $used_col; else print $2,$3,$4, $used_col}' > $output_file
#
output_file='function_label/GMP.J_IDEAS.bed'
target_cell_col=20
cat $input_IDEAS_state_file | awk -F ' ' -v OFS='\t' -v used_col=$target_cell_col '{if ($2=="CHR") print "#"$2,$3,$4, $used_col; else print $2,$3,$4, $used_col}' > $output_file
#
output_file='function_label/LSK.J_IDEAS.bed'
target_cell_col=22
cat $input_IDEAS_state_file | awk -F ' ' -v OFS='\t' -v used_col=$target_cell_col '{if ($2=="CHR") print "#"$2,$3,$4, $used_col; else print $2,$3,$4, $used_col}' > $output_file
#
output_file='function_label/MEP.J_IDEAS.bed'
target_cell_col=25
cat $input_IDEAS_state_file | awk -F ' ' -v OFS='\t' -v used_col=$target_cell_col '{if ($2=="CHR") print "#"$2,$3,$4, $used_col; else print $2,$3,$4, $used_col}' > $output_file
##########################################



