library(data.table)
library(mclust)
library(pheatmap)
library(cba)


### test optimal leaf
d1_Index = as.data.frame(fread('hg38_outputs/hg38_chrAll_analysis_merge/snapshot_test_run_merge.index.matrix.txt'))[,-c(1:4)]
d1_Index = d1_Index[rowSums(d1_Index)>0,]

ct_list = read.table('/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/peak_signal_state_list.merge.txt', header=F)[,1]
ct_list = apply(cbind(as.data.frame(ct_list)), 1, toString)
colnames(d1_Index) = ct_list

set.seed(2019)
used_rows = sample(dim(d1_Index)[1], 10000)
d1_Index_sample = d1_Index[used_rows,]

### Hclust
dist_d1_sp0 = dist(d1_Index_sample)
hclust_result = hclust(dist_d1_sp0)
hc <- hclust_result
co <- order.optimal(dist_d1_sp0, hclust_result$merge)
ho <- hc
ho$merge <- co$merge
ho$order <- co$order
d1_Index_sample_hclust_raw = d1_Index_sample[hclust_result$order ,]
d1_Index_sample_hclust_opt = d1_Index_sample[ho$order ,]

Index_set_num = 69

### Kmeans
set.seed(2019)
km = kmeans(d1_Index_sample, Index_set_num)
d1_Index_sample_km = d1_Index_sample[order(km$cluster),]

### Mclust
set.seed(2019)
mc = Mclust(d1_Index_sample, Index_set_num)
d1_Index_sample_mclust = d1_Index_sample[order(mc$classification),]

### Kmeans
set.seed(2019)
km19 = kmeans(d1_Index_sample, 19)
d1_Index_sample_km19 = d1_Index_sample[order(km19$cluster),]

### Mclust
set.seed(2019)
mc19 = Mclust(d1_Index_sample, 19)
d1_Index_sample_mclust19 = d1_Index_sample[order(mc19$classification),]

### Snapshot IS for plotting
index_set_binary = read.table('snapshot_test_run_merge.metaISid.mat.txt', header=F)
###
index_set_binary_Index = t(apply(index_set_binary, 1, function(x) unlist(strsplit(x[1], '_')) ))
index_set_binary_Index = apply(index_set_binary_Index, 2, as.numeric)
###
X_mat = index_set_binary_Index[index_set_binary[,1]=='X_X_X_X_X_X_X_X_X_X_X_X_X',]
X_mat_OD_Index = as.matrix(d1_Index[index_set_binary[,1]=='X_X_X_X_X_X_X_X_X_X_X_X_X',])
index_set_binary_Index[index_set_binary[,1]=='X_X_X_X_X_X_X_X_X_X_X_X_X',] = X_mat_OD_Index
###
index_set_binary_Index = apply(index_set_binary_Index, 2, as.numeric)
index_set_binary_Index_sample = index_set_binary_Index[used_rows,]
###
index_set_binary_sample = index_set_binary[used_rows,1]
d1_Index_sample_SP_plot = index_set_binary_Index_sample[order(index_set_binary_sample),]

### Snapshot IS
index_set_binary = read.table('snapshot_test_run_merge.metaISid.mat.txt', header=F)
index_set_binary_Index = t(apply(index_set_binary, 1, function(x) unlist(strsplit(x[1], '_')) ))
index_set_binary_Index = apply(index_set_binary_Index, 2, as.numeric)
index_set_binary_Index_sample = index_set_binary_Index[used_rows,]
###
index_set_binary_sample = index_set_binary[used_rows,1]
d1_Index_sample_SP = index_set_binary_Index_sample[order(index_set_binary_sample),]



get_entropy = function(d1_Index_sample, height){
###
d1_Index_sample_height_width = c()
#for (j in 1:(dim(d1_Index_sample)[2]-width+1)){
for (j in 2:dim(d1_Index_sample)[2]){
print(j)
### initialize pattern label vector for column j
output_length = (dim(d1_Index_sample)[1]-height+1)
d1_Index_sample_SP_ij_label_j_vec = rep('0', output_length)
###
for (i in 1:(dim(d1_Index_sample)[1]-height+1)){
### extract pattern
d1_Index_sample_SP_ij = d1_Index_sample[ i:(i+height-1) , 1:j]
### convert matrix to label
d1_Index_sample_SP_ij_label = paste(as.numeric(as.matrix(d1_Index_sample_SP_ij)), collapse='_')
d1_Index_sample_SP_ij_label_j_vec[i] = d1_Index_sample_SP_ij_label
}
### get pattern counts & p
d1_Index_sample_SP_ij_label_count = table(d1_Index_sample_SP_ij_label_j_vec)
d1_Index_sample_SP_ij_label_count_p = d1_Index_sample_SP_ij_label_count / sum(d1_Index_sample_SP_ij_label_count)
### get entropy
entropy_j = -sum(d1_Index_sample_SP_ij_label_count_p*log(d1_Index_sample_SP_ij_label_count_p))
d1_Index_sample_height_width = c(d1_Index_sample_height_width, entropy_j)
}
return(d1_Index_sample_height_width)
}

set.seed(2019)
#scanner_height = 100
#for (scanner_height in c(3, 5, 10, 50, 100, 200)){
for (scanner_height in c(5, 10, 50)){

SP_entropy_10 = get_entropy(d1_Index_sample_SP, scanner_height)
km_entropy_10 = get_entropy(d1_Index_sample_km, scanner_height)
mclust_entropy_10 = get_entropy(d1_Index_sample_mclust, scanner_height)
#km_entropy_19_10 = get_entropy(d1_Index_sample_km19, scanner_height)
#mclust_entropy_19_10 = get_entropy(d1_Index_sample_mclust19, scanner_height)
hclust_raw_entropy_10 = get_entropy(d1_Index_sample_hclust_raw, scanner_height)
hclust_opt_entropy_10 = get_entropy(d1_Index_sample_hclust_opt, scanner_height)

shuffle_entropy_10 = get_entropy(d1_Index_sample_SP[sample(dim(d1_Index_sample_SP)[1], dim(d1_Index_sample_SP)[1]),], scanner_height)

Shannon_Entropy_mat_10 = cbind(SP_entropy_10, km_entropy_10, mclust_entropy_10, km_entropy_19_10, mclust_entropy_19_10, hclust_raw_entropy_10, hclust_opt_entropy_10, shuffle_entropy_10)
colnames(Shannon_Entropy_mat_10) = c('Snapshot', 'Kmeans', 'Mclust', 'Kmeans19', 'Mclust19', 'Hclust', 'Hclust_optimal_leaf', 'Shuffle')

pdf(paste0('Shannon_Entropy_mat.',scanner_height,'.pdf'), width=4, height=4)
plot(2:dim(d1_Index_sample)[2], Shannon_Entropy_mat_10[,1], ylim=c(min(Shannon_Entropy_mat_10), max(Shannon_Entropy_mat_10)), xlab='Number of Cell-Type for SE calculation', ylab='Shannon Entropy (SE)')
lines(2:dim(d1_Index_sample)[2], Shannon_Entropy_mat_10[,1], col='dodgerblue', lwd=2)
###
points(2:dim(d1_Index_sample)[2], Shannon_Entropy_mat_10[,2], col='green')
lines(2:dim(d1_Index_sample)[2], Shannon_Entropy_mat_10[,2], col='green')
###
points(2:dim(d1_Index_sample)[2], Shannon_Entropy_mat_10[,3], col='orange')
lines(2:dim(d1_Index_sample)[2], Shannon_Entropy_mat_10[,3], col='orange')
###
#points(2:dim(d1_Index_sample)[2], Shannon_Entropy_mat_10[,4], col='springgreen4')
#lines(2:dim(d1_Index_sample)[2], Shannon_Entropy_mat_10[,4], col='springgreen4')
###
#points(2:dim(d1_Index_sample)[2], Shannon_Entropy_mat_10[,5], col='darkorange3')
#lines(2:dim(d1_Index_sample)[2], Shannon_Entropy_mat_10[,5], col='darkorange3')
###
points(2:dim(d1_Index_sample)[2], Shannon_Entropy_mat_10[,7], col='purple')
lines(2:dim(d1_Index_sample)[2], Shannon_Entropy_mat_10[,7], col='purple')
###
points(2:dim(d1_Index_sample)[2], Shannon_Entropy_mat_10[,8], col='black')
lines(2:dim(d1_Index_sample)[2], Shannon_Entropy_mat_10[,8], col='black')
dev.off()
}



library(pheatmap)

my_colorbar=colorRampPalette(c('white', 'black'))(n = length(c(1:10)))

png('hc.Nk.cCRE.optimal.heatmap.png', width=300)
pheatmap(d1_Index_sample[ho$order,], cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()

png('hc.Nk.cCRE.hclust.heatmap.png', width=300)
pheatmap(d1_Index_sample[hclust_result$order,], cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()

png('km.Nk.cCRE.heatmap.png', width=300)
pheatmap(d1_Index_sample[order(km$cluster),], cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()

png('km19.Nk.cCRE.heatmap.png', width=300)
pheatmap(d1_Index_sample[order(km19$cluster),], cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()

png('mclust.Nk.cCRE.heatmap.png', width=300)
pheatmap(d1_Index_sample[order(mc$classification),], cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()

png('mclust19.Nk.cCRE.heatmap.png', width=300)
pheatmap(d1_Index_sample[order(mc19$classification),], cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()

png('Snapshot.Nk.cCRE.optimal.heatmap.png', width=300)
#labels_s = apply(d1_Index[used_rows,],1, function(x) paste(x, collapse='_'))
pheatmap(d1_Index_sample_SP_plot, cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()



###### reorder cCREs by Index within each clusters
### Kmeans 19
d1_Index_sample_KM = c()
for (K in 1:max(km$cluster)){
	d1_Index_sample_KM_K = d1_Index_sample[km$cluster==K,]
	index_set_binary_sample_K = index_set_binary_sample[km$cluster==K]
	### Index reorder
	d1_Index_sample_KM_K_reorder = d1_Index_sample_KM_K[order(index_set_binary_sample_K),]
	d1_Index_sample_KM = rbind(d1_Index_sample_KM, d1_Index_sample_KM_K_reorder)
}

png('km.Nk.cCRE.heatmap.index_reorder.png', width=300)
pheatmap(d1_Index_sample_KM, cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()

### Kmeans 19
d1_Index_sample_KM = c()
for (K in 1:max(km19$cluster)){
	d1_Index_sample_KM_K = d1_Index_sample[km19$cluster==K,]
	index_set_binary_sample_K = index_set_binary_sample[km19$cluster==K]
	### Index reorder
	d1_Index_sample_KM_K_reorder = d1_Index_sample_KM_K[order(index_set_binary_sample_K),]
	d1_Index_sample_KM = rbind(d1_Index_sample_KM, d1_Index_sample_KM_K_reorder)
}

png('km19.Nk.cCRE.heatmap.index_reorder.png', width=300)
pheatmap(d1_Index_sample_KM, cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()

### Mclust
d1_Index_sample_Mclust = c()
for (K in 1:max(mc19$classification)){
	d1_Index_sample_KM_K = d1_Index_sample[mc19$classification==K,]
	index_set_binary_sample_K = index_set_binary_sample[mc19$classification==K]
	### Index reorder
	d1_Index_sample_KM_K_reorder = d1_Index_sample_KM_K[order(index_set_binary_sample_K),]
	d1_Index_sample_Mclust = rbind(d1_Index_sample_Mclust, d1_Index_sample_KM_K_reorder)
}

png('mclust.Nk.cCRE.heatmap.index_reorder.png', width=300)
pheatmap(d1_Index_sample_Mclust, cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()


### Hclust optimal
ho_cutree = cutree(ho, k=19)
ho_cutree_order_uniq = unique(ho_cutree[ho$order])
ho_cutree_order = (ho_cutree[ho$order])
ho_cutree_vec = c()

d1_Index_sample_Hclust_opt = c()
d1_Index_sample_Hclust_opt_OD = c()
for (K in 1:max(ho_cutree_order)){
	d1_Index_sample_KM_K = d1_Index_sample[ho$order,][ho_cutree_order==K,]
	index_set_binary_sample_K = index_set_binary_sample[ho$order][ho_cutree_order==K]
	### Index reorder
	d1_Index_sample_KM_K_reorder = d1_Index_sample_KM_K[order(index_set_binary_sample_K),]
	d1_Index_sample_Hclust_opt = rbind(d1_Index_sample_Hclust_opt, d1_Index_sample_KM_K_reorder)
	### opt reorder
	d1_Index_sample_Hclust_opt_OD = rbind(d1_Index_sample_Hclust_opt_OD, d1_Index_sample_KM_K)
	### 
	ho_cutree_vec = c(ho_cutree_vec, ho_cutree_order[ho_cutree_order==K])
}

png('hclust_opt.Nk.cCRE.heatmap.index_reorder.png', width=300)
pheatmap(d1_Index_sample_Hclust_opt, cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()

png('hclust_opt.Nk.cCRE.heatmap.opt_reorder.png', width=300)
pheatmap(d1_Index_sample_Hclust_opt_OD, cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()


### plot cluster labels
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

###
png('km19.Nk.cCRE.heatmap.index_reorder.labels.png', width=20)
pheatmap(cbind(km19$cluster[order(km19$cluster)]), cluster_cols=F, cluster_rows=F, show_rownames=F, color=col_vector, legend=F, show_colnames=F)
dev.off()
png('km.Nk.cCRE.heatmap.index_reorder.labels.png', width=20)
pheatmap(cbind(km$cluster[order(km$cluster)]), cluster_cols=F, cluster_rows=F, show_rownames=F, color=col_vector, legend=F, show_colnames=F)
dev.off()
###
png('mclust.Nk.cCRE.heatmap.index_reorder.labels.png', width=20)
pheatmap(cbind(mc19$classification[order(mc19$classification)]), cluster_cols=F, cluster_rows=F, show_rownames=F, color=col_vector, legend=F, show_colnames=F)
dev.off()
###
png('hclust_opt.Nk.cCRE.heatmap.index_reorder.labels.png', width=20)
pheatmap(cbind(ho_cutree_vec), cluster_cols=F, cluster_rows=F, show_rownames=F, color=col_vector, legend=F, show_colnames=F)
dev.off()

