library(data.table)
library(mclust)
library(pheatmap)
library(cba)

### snapshot
d1_sp0 = as.data.frame(fread('snapshot_test_run_rep1.metaISid.mat.txt'))
d2_sp0 = as.data.frame(fread('snapshot_test_run_rep2.metaISid.mat.txt'))

set.seed(2019)
ari100_sp = c()
ari20_sp = c()
for (i in 1:30){
	print(i)
used_id = sample(dim(d1_sp0)[1], dim(d1_sp0)[1])
used_id = sample(dim(d1_sp0)[1], 20000)
ari20_sp = c(ari20_sp, adjustedRandIndex(d1_sp0[used_id, 16], d2_sp0[used_id, 16]) )
ari100_sp = c(ari100_sp, adjustedRandIndex(d1_sp0[used_id, 1], d2_sp0[used_id, 1]) )
}

table(d1_sp0[,16])
table(d2_sp0[,16])
summary(ari20_sp)


### km snapshot sig
d1_Index = as.data.frame(fread('hg38_outputs/hg38_chrAll_analysis_rep1/snapshot_test_run_rep1.index.matrix.txt'))[,-c(1:4)]
d2_Index = as.data.frame(fread('hg38_outputs/hg38_chrAll_analysis_rep2/snapshot_test_run_rep2.index.matrix.txt'))[,-c(1:4)]

set.seed(2019)
Knum1 = length(unique(d1_sp0[used_id, 16]))
Knum2 = length(unique(d2_sp0[used_id, 16]))
d1_Signal_scale = d1_sp0[,-c(1,2,16)]
d2_Signal_scale = d2_sp0[,-c(1,2,16)]
ari20_km = c()
for (i in 1:30){
	print(i)
used_rows = sample(dim(d1_Signal_scale)[1], dim(d1_Signal_scale)[1])
#used_rows = sample(dim(d1_sp0)[1], 20000)
km1 = kmeans(d1_Index[used_rows,], centers=Knum1)
km2 = kmeans(d2_Index[used_rows,], centers=Knum2)
ari_b = adjustedRandIndex(km1$cluster, km2$cluster)
###
km1 = kmeans(d1_Signal_scale[used_rows,], centers=Knum1)
km2 = kmeans(d2_Signal_scale[used_rows,], centers=Knum2)
ari_s = adjustedRandIndex(km1$cluster, km2$cluster)
###
km1 = kmeans(cbind(d1_Index, d1_Signal_scale)[used_rows,], centers=Knum1)
km2 = kmeans(cbind(d2_Index, d2_Signal_scale)[used_rows,], centers=Knum2)
ari_bs = adjustedRandIndex(km1$cluster, km2$cluster)
ari20_km = rbind(ari20_km, c(ari_b, ari_s, ari_bs))
print(summary(ari20_km))
}


# Snapshot 30000 9.591

ptm <- proc.time()
used_rows = sample(dim(d1_sp0)[1], 30000)
km1 = kmeans(d1_Signal_scale[used_rows,], centers=Knum1)
proc.time() - ptm
0.425

ptm <- proc.time()
set.seed(2019)
used_rows = sample(dim(d1_sp0)[1], 30000)
dist_d1_sp0 = dist(d1_Signal_scale[used_rows,])
hclust_result = hclust(dist_d1_sp0)
hclust1 = cutree(hclust_result, k=Knum1)
proc.time() - ptm
26.284




set.seed(2019)
used_rows = sample(dim(d1_sp0)[1], 5000)
dist_d1_sp0 = dist(d1_Signal_scale[used_rows,])
hclust_result = hclust(dist_d1_sp0)
hc <- hclust_result
co <- order.optimal(dist_d1_sp0, hclust_result$merge)
ho <- hc
ho$merge <- co$merge
ho$order <- co$order


library(pheatmap)

my_colorbar=colorRampPalette(c('white', 'red'))(n = length(breaksList))

png('hc.20k.cCRE.hclust.heatmap.png')
pheatmap(hc, main="hclust")
dev.off()

png('hc.20k.cCRE.optimal.heatmap.png')
pheatmap(d1_Signal_scale[used_rows,][ho$order,], cluster_cols=F, cluster_rows=F)
dev.off()





ptm <- proc.time()
used_rows = sample(dim(d1_sp0)[1], 30000)
mclust1 = Mclust(d1_Signal_scale[used_rows,], Knum1)$classification
proc.time() - ptm
1398.387


running_time_vec = rbind(c(9.591, 0.425, 26.284, 1398.387))
colnames(running_time_vec) = c('Snapshot', 'K-means', 'Hclust', 'Mclust')
pdf('runtime.barplot.pdf', width=4, height=3)
par(mar = c(5, 4, 1, 1))
barplot(running_time_vec, xlab='', ylab='Running time (second)', cex.axis=1, width=0.1, las=2, log='')
box()
dev.off()



### test optimal leaf
d1_Index = as.data.frame(fread('hg38_outputs/hg38_chrAll_analysis_merge/snapshot_test_run_merge.index.matrix.txt'))[,-c(1:4)]
d1_Index = d1_Index[rowSums(d1_Index)>0,]

ct_list = read.table('/Users/guanjuexiang/Downloads/Snapshot_test/input_data_hg38/peak_signal_state_list.merge.txt', header=F)[,1]
ct_list = apply(cbind(as.data.frame(ct_list)), 1, toString)
colnames(d1_Index) = ct_list

set.seed(2019)
used_rows = sample(dim(d1_Index)[1], 5000)
dist_d1_sp0 = dist(d1_Index[used_rows,])
hclust_result = hclust(dist_d1_sp0)
hc <- hclust_result
co <- order.optimal(dist_d1_sp0, hclust_result$merge)
ho <- hc
ho$merge <- co$merge
ho$order <- co$order


library(pheatmap)

png('hc.20k.cCRE.optimal.heatmap.png')
pheatmap(d1_Index[used_rows,][ho$order,], cluster_cols=F, cluster_rows=F, show_rownames=F)
dev.off()

png('Snapshot.20k.cCRE.optimal.heatmap.png')
labels_s = apply(d1_Index[used_rows,],1, function(x) paste(x, collapse='_'))
pheatmap(d1_Index[used_rows,][order(labels_s),], cluster_cols=F, cluster_rows=F, show_rownames=F)
dev.off()



### DPGP
#conda create --name dpgp python GPy pandas numpy scipy matplotlib



d1 = read.table('snapshot_test_run_merge.index.matrix.allpk.txt', header=F)

non0row = apply(d1, 1, function(x) sum(as.numeric(x[-c(1:4)]))>0 )

d1a = d1[non0row,]
d1a[,4] = 1:dim(d1a)[1]
write.table(d1a, 'snapshot_test_run_merge.index.matrix.txt', quote=F, sep='\t', col.names=F, row.names=F)

d1 = read.table('snapshot_test_run_merge.signal.matrix.allpk.txt', header=F)
d1a = d1[non0row,]
d1a[,4] = 1:dim(d1a)[1]
write.table(d1a, 'snapshot_test_run_merge.signal.matrix.txt', quote=F, sep='\t', col.names=F, row.names=F)

d1 = read.table('snapshot_test_run_merge.function.matrix.allpk.txt', header=F)
d1a = d1[non0row,]
d1a[,4] = 1:dim(d1a)[1]
write.table(d1a, 'snapshot_test_run_merge.function.matrix.txt', quote=F, sep='\t', col.names=F, row.names=F)






