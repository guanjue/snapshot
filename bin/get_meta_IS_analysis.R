library(pheatmap)
library(tidyverse)
library(mclust)

binary_mat = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_merge/snapshot_test_run_merge.index_binary_mat.txt'
sig_mat = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_merge/snapshot_test_run_merge.meansig.txt'
signal_input_list = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/peak_signal_state_list.merge.txt'
sig_cCRE = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_merge/snapshot_test_run_merge.sig.txt'
ct_cluster_num = 8

binary_mat = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_rep1/snapshot_test_run_rep1.index_binary_mat.txt'
sig_mat = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_rep1/snapshot_test_run_rep1.meansig.txt'
signal_input_list = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/peak_signal_state_list.rep1.txt'
sig_cCRE = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_rep1/snapshot_test_run_rep1.sig.txt'
ct_cluster_num = 8

binary_mat = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_rep2/snapshot_test_run_rep2.index_binary_mat.txt'
sig_mat = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_rep2/snapshot_test_run_rep2.meansig.txt'
signal_input_list = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/peak_signal_state_list.rep2.txt'
sig_cCRE = '/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_rep2/snapshot_test_run_rep2.sig.txt'
ct_cluster_num = 8

### read binary & signal matrix file
d_binary = read.table(binary_mat, header=F)
d_sig = read.table(sig_mat, header=F)
d_sig_cCRE = read.table(sig_cCRE, header=F)

### read colnames file
colname_file = read.table(signal_input_list, header=F, sep='\t')
colname_file[,1] = apply(colname_file, 1, function(x) as.character(x[1]) )
### add colnames
colnames(d_binary)[-1] = colname_file[,1]
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
for (k in 2:(dim(input_mat)[2]-1)){
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

pdf('ct_cluster_num_vs_within_cluster_dist.pdf')
plot(2:(dim(input_mat)[2]-1), within_cluster_dist_all)
lines(2:(dim(input_mat)[2]-1), within_cluster_dist_all)
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


d_sig_merge_metaIS = cutree(hclust(dist(d_sig_merge)), 20)


d_sig_merge_meta_cluster = c()
for (k in unique(d_sig_merge_metaIS)){
	d_sig_merge_meta_cluster = rbind(d_sig_merge_meta_cluster, colMeans(input_mat[d_sig_merge_metaIS==k,]))
}

pdf('meta_cluster_ave.pdf', height=15)
pheatmap(input_mat[order(d_sig_merge_metaIS),], cluster_col=F, cluster_row=F,cex=1.8, show_rownames=F)
dev.off()

pdf('meta_cluster_ave_merge.pdf', width = 5)
d_sig_merge_meta_cluster_plot = d_sig_merge_meta_cluster
breaksList = seq(0, 16, by = 0.01)
my_colorbar=colorRampPalette(c('white', 'red'))(n = length(breaksList))
pheatmap(d_sig_merge_meta_cluster_plot, cluster_col=F, cluster_row=F,cex=1.5, color=my_colorbar, breaks = breaksList)
dev.off()

library(RColorBrewer)
n <- 120
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

pdf('meta_cluster_label.pdf', width=3)
pheatmap(cbind(cbind(d_sig_merge_metaIS)[order(d_sig_merge_metaIS),]), cluster_col=F, cluster_row=F, color=col_vector)
dev.off()


IStoISmeta = cbind(as.data.frame(d_sig[,1]), d_sig_merge_metaIS)
colnames(IStoISmeta) = c('ISid','metaISid')
colnames(d_sig_cCRE)[2] = 'ISid'

shared=merge(d_sig_cCRE, IStoISmeta, by='ISid')

shared_rep2 = shared

shared_rep1 = shared


used_id = sample(dim(shared_rep1)[1], 5000)
adjustedRandIndex(shared_rep1[order(shared_rep1[,2]),16][used_id], shared_rep2[order(shared_rep2[,2]),16][used_id])
adjustedRandIndex(shared_rep1[order(shared_rep1[,2]),1][used_id], shared_rep2[order(shared_rep2[,2]),1][used_id])

#d0 = read.table('snapshot_test_run_rep2.signal.matrix.txt', header=F)[,-c(1:4)]

d0 = read.table('snapshot_test_run_merge.sig.txt', header=F)[,-c(1:2)]

ref_colmean = median(colMeans(d0))
d0 = apply(d0,2,function(x) x/mean(x)*ref_colmean )

d0s = d0[sample(dim(d0)[1], 5000),]

d0s_dist = 1-cor(t(d0s))

d0s_hclust = hclust(dist(d0s))

d0s_hclust = hclust(as.dist(d0s_dist))

d0s_hclust_cutree = cutree(d0s_hclust, 20)

d0_cluster = c()
for (k in unique(d0s_hclust_cutree)){
	d0_cluster = rbind(d0_cluster, colMeans(d0s[d0s_hclust_cutree==k,]))
}

colnames(d0_cluster) = colnames(d_sig_merge_meta_cluster)
pheatmap(d0_cluster, cluster_col=F)


d0s = d0[sample(dim(d0)[1], 5000),]
d0s_kmeans = kmeans(d0s, 20)
d0_KMcluster = c()
for (k in unique(d0s_kmeans$cluster)){
	d0_KMcluster = rbind(d0_KMcluster, colMeans(d0s[d0s_kmeans$cluster==k,]))
}

colnames(d0_KMcluster) = colnames(d_sig_merge_meta_cluster)
pheatmap(d0_KMcluster, cluster_col=F)



kmeans_cor <- function(data, k, max_iter) {
  set.seed(2019)
  #print('# Determine the number of data points and features')
  n <- nrow(data)
  p <- ncol(data)
  
  #print('# Initialize the centroids randomly')
  centroids <- data[sample(1:n, k),]
  
  #print('# Initialize the cluster assignments')
  clusters <- rep(0, n)
  
  #print('# Iterate over the maximum number of iterations')
  for (i in 1:max_iter) {
    print(i)
    clusters_pre = clusters
    # Assign each data point to the closest centroid
    distances_mat = 1-cor(t(data), t(centroids))
    clusters = apply(distances_mat, 1, which.min)    
    # Update the centroids
    for (j in 1:k) {
      centroids[j,] <- colMeans(data[clusters == j,])
    }
    #print(table(clusters))
    if (sum(clusters==clusters_pre)==n){break}
  }
  
  #print('# Return the cluster assignments')
  output = list()
  output$centroids = centroids
  output$clusters = clusters
  return(output) 
}


d0s = d0[sample(dim(d0)[1], 5000),]
d0s_kmeans = kmeans(d0s, 20, 10000)
d0_KMcluster = c()
for (k in unique(d0s_kmeans$cluster)){
	d0_KMcluster = rbind(d0_KMcluster, colMeans(d0s[d0s_kmeans$cluster==k,]))
}

colnames(d0_KMcluster) = colnames(d_sig_merge_meta_cluster)
pheatmap(d0_KMcluster, cluster_col=F)


get_km_2step = function(d0s, km_n1, km_n2){
	d0s = d1[used_id,]
d0s_kmeans = kmeans(d0s, km_n1, 10000)
d0_KMcluster = c()
d0s_kmeans_cluster = d0s_kmeans$cluster
for (k in unique(d0s_kmeans_cluster)){
	d0_KMcluster = rbind(d0_KMcluster, colMeans(d0s[d0s_kmeans$cluster==k,]))
}
###
d0s_kmeans_meta = kmeans(d0_KMcluster, km_n2, 10000)
IStoISmeta = cbind(unique(d0s_kmeans_cluster), d0s_kmeans_meta$cluster)
colnames(IStoISmeta) = c('ISid','metaISid')
ISvec = cbind(1:dim(d0s)[1], d0s_kmeans$cluster)
colnames(ISvec) = c('pkid', 'ISid')
cCRE_ISmetaID = merge(ISvec, IStoISmeta, by='ISid', sort=F)
return(cCRE_ISmetaID[order(cCRE_ISmetaID[,2]),-2])
}

set.seed(2019)
d1aa = read.table('/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_rep1/snapshot_test_run_rep1.sig.txt', header=F)
d2aa = read.table('/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/hg38_outputs/hg38_chrAll_analysis_rep2/snapshot_test_run_rep2.sig.txt', header=F)

d1 = d1aa[order(d1aa[,1]),-c(1,2)]
d2 = d2aa[order(d2aa[,1]),-c(1,2)]


km_1step = c()

for (i in 1:30){
	print(i)
	used_id = sample(dim(d1)[1], 5000)
cCRE_ISmetaID_1 = get_km_2step(d1[used_id,], 20, 10)
cCRE_ISmetaID_2 = get_km_2step(d2[used_id,], 20, 10)
km_1step = rbind(km_1step, c(adjustedRandIndex(cCRE_ISmetaID_1[,1], cCRE_ISmetaID_2[,1]), adjustedRandIndex(cCRE_ISmetaID_1[,2], cCRE_ISmetaID_2[,2]) ))
}




km_2step = c()

for (i in 1:30){
	print(i)
	used_id = sample(dim(d1)[1], 5000)
cCRE_ISmetaID_1 = get_km_2step(d1[used_id,], 83, 20)
cCRE_ISmetaID_2 = get_km_2step(d2[used_id,], 143, 20)
km_2step = rbind(km_2step, c(adjustedRandIndex(cCRE_ISmetaID_1[,1], cCRE_ISmetaID_2[,1]), adjustedRandIndex(cCRE_ISmetaID_1[,2], cCRE_ISmetaID_2[,2]) ))
}




set.seed(2019)
d0s = d0[sample(dim(d0)[1], 5000),]
d0s_kmeans = kmeans(d0s, 105, 10000)
d0_KMcluster = c()
for (k in unique(d0s_kmeans$cluster)){
	d0_KMcluster = rbind(d0_KMcluster, colMeans(d0s[d0s_kmeans$cluster==k,]))
}

colnames(d0_KMcluster) = colnames(d_sig_merge_meta_cluster)
pdf('try_KM.105.pdf')
pheatmap(d0_KMcluster, cluster_col=F)
dev.off()


d0s_kmeans_meta = kmeans(d0_KMcluster, 20, 10000)

d0_KMcluster_meta = c()
for (k in unique(d0s_kmeans_meta$cluster)){
	if (sum(d0s_kmeans_meta$cluster==k)>1){
	d0_KMcluster_meta = rbind(d0_KMcluster_meta, colMeans(d0_KMcluster[d0s_kmeans_meta$cluster==k,]))
	} else{
		d0_KMcluster_meta = rbind(d0_KMcluster_meta, (d0_KMcluster[d0s_kmeans_meta$cluster==k,]))
	}
}

pdf('try_KM_meta.20.pdf')
pheatmap(d0_KMcluster_meta, cluster_col=F)
dev.off()






