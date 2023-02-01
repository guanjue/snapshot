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
dist_d1_sp0 = dist(d1_Index[used_rows,])
hclust_result = hclust(dist_d1_sp0)
hc <- hclust_result
co <- order.optimal(dist_d1_sp0, hclust_result$merge)
ho <- hc
ho$merge <- co$merge
ho$order <- co$order

set.seed(2019)
km = kmeans(d1_Index[used_rows,], 69)

library(pheatmap)

my_colorbar=colorRampPalette(c('white', 'black'))(n = length(c(1:10)))

png('hc.Nk.cCRE.optimal.heatmap.png', width=300)
pheatmap(d1_Index[used_rows,][ho$order,], cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()

png('hc.Nk.cCRE.hclust.heatmap.png', width=300)
pheatmap(d1_Index[used_rows,][hclust_result$order,], cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()

png('km.Nk.cCRE.heatmap.png', width=300)
pheatmap(d1_Index[used_rows,][order(km$cluster),], cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()


index_set_binary = read.table('snapshot_test_run_merge.metaISid.mat.txt', header=F)
labels_s = index_set_binary[used_rows,1]
png('Snapshot.Nk.cCRE.optimal.heatmap.png', width=300)
#labels_s = apply(d1_Index[used_rows,],1, function(x) paste(x, collapse='_'))
pheatmap(d1_Index[used_rows,][order(labels_s),], cluster_cols=F, cluster_rows=F, show_rownames=F, color=my_colorbar, legend=F, show_colnames=F)
dev.off()











