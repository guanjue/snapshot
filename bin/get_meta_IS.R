library(pheatmap)
library(tidyverse)

####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
sig_mat = args[1]
signal_input_list = args[2]
sig_cCRE = args[3]
outputname = args[4]
cCRE_list = args[5]
function_matrix = args[6]
have_function_state_files = args[7]
#ct_cluster_num = as.numeric(args[5])
#meta_IS_num = as.numeric(args[6])


#binary_mat = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_merge/snapshot_test_run_merge.index_binary_mat.txt'
#sig_mat = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_merge/snapshot_test_run_merge.meansig.txt'
#signal_input_list = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/peak_signal_state_list.merge.txt'
#sig_cCRE = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_merge/snapshot_test_run_merge.sig.txt'
#ct_cluster_num = 8

#binary_mat = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_rep1/snapshot_test_run_rep1.index_binary_mat.txt'
#sig_mat = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_rep1/snapshot_test_run_rep1.meansig.txt'
#signal_input_list = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/peak_signal_state_list.rep1.txt'
#sig_cCRE = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_rep1/snapshot_test_run_rep1.sig.txt'
#ct_cluster_num = 8

#binary_mat = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_rep2/snapshot_test_run_rep2.index_binary_mat.txt'
#sig_mat = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_rep2/snapshot_test_run_rep2.meansig.txt'
#signal_input_list = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/peak_signal_state_list.rep2.txt'
#sig_cCRE = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_rep2/snapshot_test_run_rep2.sig.txt'
#ct_cluster_num = 8
#meta_IS_num = 20


#sig_mat = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_merge/snapshot_test_run_merge.meansig.txt'
#signal_input_list = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/peak_signal_state_list.merge.txt'
#sig_cCRE = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_merge/snapshot_test_run_merge.sig.txt'
#outputname = 'snapshot_test_run_merge'
#cCRE_list = 'hg38_outputs/hg38_chrAll_analysis_merge/snapshot_test_run_merge.sort.bed'

####################################
### read binary & signal matrix file
d_ISmeansig = read.table(sig_mat, header=F)
d_sig_cCRE = read.table(sig_cCRE, header=F)
cCRE_info = read.table(cCRE_list, header=F)

### read colnames file
colname_file = read.table(signal_input_list, header=F, sep='\t')
colname_file[,1] = apply(colname_file, 1, function(x) as.character(x[1]) )
### add colnames
colnames(d_ISmeansig)[-1] = colname_file[,1]
### cluster cell types
input_mat = d_ISmeansig[,-1]
####################################


####################################
### hclust ISs
IS_dist = dist(input_mat)
d_sig_hclust = hclust(IS_dist)
####################################


####################################
### determine meta-IS set cluster number
AIC_cluster_dist_all_IS = c()
for (k in 2:(dim(input_mat)[1])){
	#for (k in 2:50){
	IS_cluster = cutree(d_sig_hclust, k)
	within_cluster_dist_i = c()
	for (j in 1:k){
		### get input_mat_j
		d_sig_j = input_mat[IS_cluster==j,]
		### get center
		if (!is.null(dim(d_sig_j))){
			d_sig_j_center = colMeans(d_sig_j)
			dist_within = apply(d_sig_j, 1, function(x) sum((x-d_sig_j_center)^2) )
		} else{
			d_sig_j_center = d_sig_j
			dist_within = 0
		}
		within_cluster_dist_i = c(within_cluster_dist_i, dist_within )
	}
	totss = sum(within_cluster_dist_i)
	aic_k = totss + 2*k*dim(input_mat)[2]
	AIC_cluster_dist_all_IS = c(AIC_cluster_dist_all_IS, aic_k)
}

pdf(paste0(outputname, '.IS_cluster.AIC.pdf'), width=5, height=5)
plot(2:(dim(input_mat)[1]), AIC_cluster_dist_all_IS, xlab='IS_SET_Number')
lines(2:(dim(input_mat)[1]), AIC_cluster_dist_all_IS)
meta_IS_num = (2:(dim(input_mat)[1]))[which.min(AIC_cluster_dist_all_IS)]
print(paste0('meta_IS_num: ', meta_IS_num))
abline(v=meta_IS_num)
dev.off()
####################################


####################################
### final merge Meta-IS clusters
d_sig_merge_metaIS = cutree(hclust(dist(input_mat)), meta_IS_num)
d_sig_merge_meta_cluster = c()
for (k in unique(d_sig_merge_metaIS)){
	d_sig_merge_meta_cluster = rbind(d_sig_merge_meta_cluster, colMeans(input_mat[d_sig_merge_metaIS==k,]))
}
####################################


pdf(paste0(outputname, '.meta_cluster_ave.pdf'), height=15)
pheatmap(input_mat[order(d_sig_merge_metaIS),], cluster_col=F, cluster_row=F,cex=1.8, show_rownames=F)
dev.off()


library(RColorBrewer)
n <- 120
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

pdf(paste0(outputname, '.meta_cluster_label.pdf'), width=3)
pheatmap(cbind(cbind(d_sig_merge_metaIS)[order(d_sig_merge_metaIS),]), cluster_col=F, cluster_row=F, color=col_vector)
dev.off()

###
IStoISmeta = cbind(as.data.frame(d_ISmeansig[,1]), d_sig_merge_metaIS)
colnames(IStoISmeta) = c('ISid','metaISid')
colnames(d_sig_cCRE)[2] = 'ISid'

### merged by ISid
shared=merge(d_sig_cCRE, IStoISmeta, by='ISid')
shared_sig = shared[,-c(1,2,dim(shared)[2])]
shared_sig_mean = c()
for (i in 1:max(shared[,dim(shared)[2]])){
	shared_sig_mean = rbind(shared_sig_mean, colMeans(shared_sig[shared[,dim(shared)[2]]==i,]))
}
colnames(shared_sig_mean) = colname_file[,1]
pdf(paste0(outputname, '.meta_cluster_cCRE_ave_merge.pdf'), width = 6)
breaksList = seq(0, (max(shared_sig_mean)), by = 0.01)
my_colorbar=colorRampPalette(c('white', 'red'))(n = length(breaksList))
rownames(shared_sig_mean) = 1:dim(shared_sig_mean)[1]
pheatmap(shared_sig_mean, cluster_col=F, cluster_row=F,cex=1.5, color=my_colorbar, breaks = breaksList)
dev.off()


### write MetaISid matrix
Final_output_mat = shared[order(as.matrix(shared[,2])),]
write.table(Final_output_mat, paste0(outputname, '.metaISid.mat.txt'), quote=F, sep='\t', col.names=F, row.names=F)

### write signal matrix for violin plot
write.table(cbind(Final_output_mat[,c(2,dim(Final_output_mat)[2], 3:(dim(Final_output_mat)[2]-1) )]), paste0(outputname, '.MetaIS.forviolin.sig.txt'), quote=F, sep='\t', col.names=F, row.names=F)

### write MetaIS meansignal matrix
MetaIS_id_pre_cCRE_vec = Final_output_mat[,dim(Final_output_mat)[2]]
unique_MetaIS = unique(MetaIS_id_pre_cCRE_vec)
unique_MetaIS = unique_MetaIS[order(unique_MetaIS)]
MetaIS_meansig_mat = c()
for (metaIS_i in unique_MetaIS){
	MetaIS_meansig_mat_i = Final_output_mat[MetaIS_id_pre_cCRE_vec==metaIS_i,-c(1,2, dim(Final_output_mat)[2])]
	MetaIS_meansig_mat = rbind(MetaIS_meansig_mat, colMeans(MetaIS_meansig_mat_i))
}
MetaIS_meansig_mat = cbind(unique_MetaIS, MetaIS_meansig_mat)
write.table(MetaIS_meansig_mat, paste0(outputname, '.metaISid.meansigmat.txt'), quote=F, sep='\t', col.names=F, row.names=F)

### write MetaIS Functional State matrix
if (have_function_state_files!='F'){
	function_state_mat = read.table(function_matrix, header=F, sep='\t')
	function_state_mat_withMetaISID = cbind(as.data.frame(function_state_mat[,4]), MetaIS_id_pre_cCRE_vec, function_state_mat[,-c(1:4)] )
	write.table(function_state_mat_withMetaISID, paste0(outputname, '.metaIS_all.fun.txt'), quote=F, sep='\t', col.names=F, row.names=F)
}


### reorgnize matrix for output
Final_output_mat = shared[order(as.matrix(shared[,2])),]
Final_output_mat = cbind(cCRE_info[,-4], Final_output_mat)
Final_output_mat = cbind(Final_output_mat[,c(1:3,5,4)], Final_output_mat[,dim(Final_output_mat)[2]], Final_output_mat[,-c(1:5, dim(Final_output_mat)[2])])
colnames(Final_output_mat) = c('#chr', 'start', 'end', 'cCREsID', 'IndexSetID', 'MetaISID', colname_file[,1])
### 
write.table(Final_output_mat, paste0(outputname, '.IS_metaISid.mat.final.txt'), quote=F, sep='\t', col.names=T, row.names=F)


###### Write Bed files for each Meta-ISs
output_dir_MetaISID = paste0(outputname, '_MetaISs_bed_files')
dir.create(output_dir_MetaISID, showWarnings=F)
for (MetaISID_i in unique(Final_output_mat$MetaISID)){
	bed_info_MetaIS_i = Final_output_mat[Final_output_mat$MetaISID==MetaISID_i,c(1:4,6)]
	write.table(bed_info_MetaIS_i, paste0(output_dir_MetaISID, '/MetaIS.', MetaISID_i, '.bed'), quote=F, sep='\t', col.names=F, row.names=F)
}


