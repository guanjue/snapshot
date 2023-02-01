library(pheatmap)
library(tidyverse)

####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
sig_mat = args[1]
signal_input_list = args[2]
sig_cCRE = args[3]
outputname = args[4]
ct_cluster_num = as.numeric(args[5])
meta_IS_num = as.numeric(args[6])


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

### read binary & signal matrix file
d_sig = read.table(sig_mat, header=F)
d_sig_cCRE = read.table(sig_cCRE, header=F)

### read colnames file
colname_file = read.table(signal_input_list, header=F, sep='\t')
colname_file[,1] = apply(colname_file, 1, function(x) as.character(x[1]) )
### add colnames
colnames(d_sig)[-1] = colname_file[,1]

### cluster cell types
input_mat = d_sig[,-1]
#ct_dist = dist(t(input_mat))
#d_sig_hclust = hclust(ct_dist)
ct_dist = as.dist(1-cor(input_mat))
d_sig_hclust = hclust(ct_dist)
ct_cluster = cutree(d_sig_hclust, ct_cluster_num)
print(cbind(colnames(d_sig)[-1], ct_cluster)[order(ct_cluster),])

### determine cell-type cluster number
within_cluster_dist_all = c()
for (k in 2:(dim(input_mat)[2])){
ct_cluster = cutree(d_sig_hclust, k)
within_cluster_dist_i = c()
for (j in 1:k){
### get input_mat_j
d_sig_j = input_mat[,ct_cluster==j]
### get center
if (!is.null(dim(d_sig_j))){
d_sig_j_center = rowMeans(d_sig_j)
} else{
	d_sig_j_center = d_sig_j
}
within_cluster_dist_i = c(within_cluster_dist_i, 1-cor(cbind(d_sig_j_center), d_sig_j) )
}
within_cluster_dist_all = c(within_cluster_dist_all, mean(within_cluster_dist_i))
}

pdf(paste0(outputname, '.ct_cluster_num_vs_within_cluster_dist.pdf'))
plot(2:(dim(input_mat)[2]), within_cluster_dist_all)
lines(2:(dim(input_mat)[2]), within_cluster_dist_all)
abline(v=ct_cluster_num)
dev.off()


### get new signal matrix
d_sig_merge = matrix(0, nrow=dim(d_sig)[1], ncol=ct_cluster_num)
for (k in 1:ct_cluster_num){
	d_sig_k = input_mat[,ct_cluster==k]
	if (!is.null(dim(d_sig_k))){
		d_sig_merge[,k] = rowMeans(d_sig_k)
	} else{
		d_sig_merge[,k] = d_sig_k
	}
}

### hclust
IS_dist = dist(d_sig_merge)
d_sig_hclust = hclust(IS_dist)



### determine meta-IS set cluster number
AIC_cluster_dist_all_IS = c()
for (k in 2:(dim(d_sig_merge)[1])){
#for (k in 2:50){
IS_cluster = cutree(d_sig_hclust, k)
within_cluster_dist_i = c()
for (j in 1:k){
### get input_mat_j
d_sig_j = d_sig_merge[IS_cluster==j,]
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
#bic_k = totss + log(dim(d_sig_merge)[1])*k*dim(d_sig_merge)[2]
aic_k = totss + 2*k*dim(d_sig_merge)[2]
AIC_cluster_dist_all_IS = c(AIC_cluster_dist_all_IS, aic_k)
}

pdf(paste0(outputname, '.IS_cluster.AIC.pdf'), width=5, height=5)
plot(2:(dim(input_mat)[1]), AIC_cluster_dist_all_IS, xlab='IS_SET_Number')
lines(2:(dim(input_mat)[1]), AIC_cluster_dist_all_IS)
meta_IS_num = (2:(dim(d_sig_merge)[1]))[which.min(AIC_cluster_dist_all_IS)]
print(paste0('meta_IS_num: ', meta_IS_num))
abline(v=meta_IS_num)
dev.off()



### final merge Meta-IS clusters
d_sig_merge_metaIS = cutree(hclust(dist(d_sig_merge)), meta_IS_num)

d_sig_merge_meta_cluster = c()
for (k in unique(d_sig_merge_metaIS)){
	d_sig_merge_meta_cluster = rbind(d_sig_merge_meta_cluster, colMeans(input_mat[d_sig_merge_metaIS==k,]))
}

pdf(paste0(outputname, '.meta_cluster_ave.pdf'), height=15)
pheatmap(input_mat[order(d_sig_merge_metaIS),], cluster_col=F, cluster_row=F,cex=1.8, show_rownames=F)
dev.off()

pdf(paste0(outputname, '.meta_cluster_ave_merge.pdf'), width = 5)
d_sig_merge_meta_cluster_plot = d_sig_merge_meta_cluster
breaksList = seq(0, 16, by = 0.01)
my_colorbar=colorRampPalette(c('white', 'red'))(n = length(breaksList))
pheatmap(d_sig_merge_meta_cluster_plot, cluster_col=F, cluster_row=F,cex=1.5, color=my_colorbar, breaks = breaksList)
dev.off()

library(RColorBrewer)
n <- 120
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

pdf(paste0(outputname, '.meta_cluster_label.pdf'), width=3)
pheatmap(cbind(cbind(d_sig_merge_metaIS)[order(d_sig_merge_metaIS),]), cluster_col=F, cluster_row=F, color=col_vector)
dev.off()


IStoISmeta = cbind(as.data.frame(d_sig[,1]), d_sig_merge_metaIS)
colnames(IStoISmeta) = c('ISid','metaISid')
colnames(d_sig_cCRE)[2] = 'ISid'

shared=merge(d_sig_cCRE, IStoISmeta, by='ISid')

print(head(shared))

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


write.table(shared[order(as.matrix(shared[,2])),], paste0(outputname, '.metaISid.mat.txt'), quote=F, sep='\t', col.names=F, row.names=F)






