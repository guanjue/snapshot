### get VISION ATAC/DNase peak bed file
cd input_dat/atac_pk/
while read ct link filename
do
	wget $link
	mv $filename $ct'.pk.bed'
done < peak.link.txt
cd -

### get VISION ATAC/DNase signal bigWig file
cd input_dat/atac_sig/
for link in $(cat signal_rc_bw.link.txt)
do
	wget $link
done 
cd -

### get IDEAS epigenetic state bigBed file
cd input_dat/function_label/
for link in $(cat statebb.link.txt)
do
	wget $link
done 
cd -

### get master cCRE list
cd input_dat/
wget https://usevision.org/data/mm10/ideasJointMay2021/S3V2_IDEAS_mm10_ccre2.cCRE.M.notall0.bed
cd -


