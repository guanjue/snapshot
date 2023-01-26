library(data.table)
library(mclust)

d_1e_2 = as.data.frame(fread('snapshot_test_run_merge.metaISid.mat.1e_2.txt'))
d_1e_1 = as.data.frame(fread('snapshot_test_run_merge.metaISid.mat.1e_1.txt'))
d_02 = as.data.frame(fread('snapshot_test_run_merge.metaISid.mat.02.txt'))

d_1e_2_90 = as.data.frame(fread('snapshot_test_run_merge.metaISid.mat.1e_2.90.txt'))
d_1e_2_85 = as.data.frame(fread('snapshot_test_run_merge.metaISid.mat.1e_2.85.txt'))


ari_mat = c()
set.seed(2019)

for (i in 1:30){
	print(i)
### used rows
used_id = sample(dim(d_1e_2)[1], 10000)
ari1 = adjustedRandIndex(d_1e_1[used_id,16], d_1e_2[used_id,16])
ari2 = adjustedRandIndex(d_02[used_id,16], d_1e_2[used_id,16])
ari3 = adjustedRandIndex(d_1e_2_90[used_id,16], d_1e_2[used_id,16])
ari4 = adjustedRandIndex(d_1e_2_85[used_id,16], d_1e_2[used_id,16])
###
km_1 = kmeans(d_1e_1[used_id,-c(1:2)], length(unique(d_1e_1[used_id,1])))
km_2 = kmeans(d_1e_2[used_id,-c(1:2)], length(unique(d_1e_2[used_id,1])))
ari5 = adjustedRandIndex(km_1$cluster, km_2$cluster)
###
km_1 = kmeans(d_02[used_id,-c(1:2)], length(unique(d_02[used_id,1])))
km_2 = kmeans(d_1e_2[used_id,-c(1:2)], length(unique(d_1e_2[used_id,1])))
ari6 = adjustedRandIndex(km_1$cluster, km_2$cluster)
###
km_1 = kmeans(d_1e_2_90[used_id,-c(1:2)], length(unique(d_1e_2_90[used_id,1])))
km_2 = kmeans(d_1e_2[used_id,-c(1:2)], length(unique(d_1e_2[used_id,1])))
ari7 = adjustedRandIndex(km_1$cluster, km_2$cluster)
###
km_1 = kmeans(d_1e_2_85[used_id,-c(1:2)], length(unique(d_1e_2_85[used_id,1])))
km_2 = kmeans(d_1e_2[used_id,-c(1:2)], length(unique(d_1e_2[used_id,1])))
ari8 = adjustedRandIndex(km_1$cluster, km_2$cluster)
###
ari_mat = rbind(ari_mat, c(ari1, ari2, ari3, ari4, ari5, ari6, ari7, ari8))
print(ari_mat)
}


colnames(ari_mat) = c('FDR:1e-1/1e-2', 'FDR:0.2/1e-2', 'TOP:10%/5%', 'TOP:15%/5%', 'FDR:1e-1/1e-2 NC', 'FDR:0.2/1e-2 NC', 'TOP:10%/5% NC', 'TOP:15%/5% NC')

pdf('ARI.difFDR.boxplot.pdf', width=4, height=4)
par(mar = c(9, 4, 1, 1))
boxplot(ari_mat[,c(1:3, 5:7)], cex.axis=1, las=2, ylab='Adjusted Rand Index', ylim=c(0,1))
dev.off()

