library(data.table)
library(mclust)
library(pheatmap)

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
used_rows = sample(dim(d1_sp0)[1], 30000)
hclust1 = cutree(hclust(dist(d1_Signal_scale[used_rows,])), k=Knum1)
proc.time() - ptm
26.284

ptm <- proc.time()
used_rows = sample(dim(d1_sp0)[1], 30000)
mclust1 = Mclust(d1_Signal_scale[used_rows,], Knum1)$classification
proc.time() - ptm
1301.633


running_time_vec = rbind(c(9.591, 0.425, 26.284, 1301.633))
colnames(running_time_vec) = c('Snapshot', 'K-means', 'Hclust', 'Mclust')
pdf('runtime.barplot.pdf', width=4, height=4)
barplot(running_time_vec, xlab='', ylab='Running time', cex.axis=1, width=0.5, las=2)
box()
dev.off()

### Hclust snapshot sig
set.seed(2019)
ptm <- proc.time()
Knum1 = length(unique(d1_sp0[used_id, 16]))
Knum2 = length(unique(d2_sp0[used_id, 16]))
d1_Signal_scale = d1_sp0[,-c(1,2,16)]
d2_Signal_scale = d2_sp0[,-c(1,2,16)]
ari20_hclust = c()
for (i in 1:30){
	print(i)
used_rows = sample(dim(d1_Index)[1], 20000)
hclust1 = cutree(hclust(dist(d1_Index[used_rows,])), k=Knum1)
hclust2 = cutree(hclust(dist(d2_Index[used_rows,])), k=Knum2)
ari_b = adjustedRandIndex(hclust1, hclust2)
###
#used_rows = sample(dim(d1_Signal_scale)[1], 5000)
hclust1 = cutree(hclust(dist(d1_Signal_scale[used_rows,])), k=Knum1)
hclust2 = cutree(hclust(dist(d2_Signal_scale[used_rows,])), k=Knum2)
ari_s = adjustedRandIndex(hclust1, hclust2)
###
#used_rows = sample(dim(d1_Signal_scale)[1], 5000)
hclust1 = cutree(hclust(dist(cbind(d1_Index, d1_Signal_scale)[used_rows,])), k=Knum1)
hclust2 = cutree(hclust(dist(cbind(d2_Index, d2_Signal_scale)[used_rows,])), k=Knum2)
ari_bs = adjustedRandIndex(hclust1, hclust2)
ari20_hclust = rbind(ari20_hclust, c(ari_b, ari_s, ari_bs))
print(summary(ari20_hclust))
}
proc.time() - ptm


### 
pdf('ari20.boxplot.snapshot_sig.pdf')
par(mar = c(10, 4, 4, 4))
ari_mat = cbind(cbind(as.numeric(ari20_sp)), ari20_hclust, ari20_km)
colnames(ari_mat) = c('Snapshot', 'Hclust_Index', 'Hclust_Sig', 'Hclust_Both', 'KM_Index', 'KM_Sig', 'KM_Both')
ari_mat_plot = ari_mat[,c(2,5, 3,6, 4,7, 1)]
ari_mat_plot = ari_mat[,c(2:7, 1)]
boxplot(ari_mat_plot, las=2, cex.axis=1.5)
dev.off()
#    user   system  elapsed 
#1907.192  258.481 2321.164


### Mclust snapshot sig
set.seed(2019)
ptm <- proc.time()
Knum1 = length(unique(d1_sp0[used_id, 16]))
Knum2 = length(unique(d2_sp0[used_id, 16]))
d1_Signal_scale = d1_sp0[,-c(1,2,16)]
d2_Signal_scale = d2_sp0[,-c(1,2,16)]
ari20_Mclust = c()
for (i in 1:30){
	print(i)
used_rows = sample(dim(d1_Index)[1], 10000)
mclust1 = Mclust(d1_Index[used_rows,]+matrix(runif(dim(d1_Index)[1]*dim(d1_Index)[2], -0.01, 0.01), dim(d1_Index)[1]), Knum1)$classification
mclust2 = Mclust(d2_Index[used_rows,]+matrix(runif(dim(d1_Index)[1]*dim(d1_Index)[2], -0.01, 0.01), dim(d1_Index)[1]), Knum1)$classification
ari_b = adjustedRandIndex(mclust1, mclust2)
###
#used_rows = sample(dim(d1_Signal_scale)[1], 5000)
mclust1 = Mclust(d1_Signal_scale[used_rows,]+matrix(runif(dim(d1_Index)[1]*dim(d1_Index)[2], -0.01, 0.01), dim(d1_Index)[1]), Knum1)$classification
mclust2 = Mclust(d2_Signal_scale[used_rows,]+matrix(runif(dim(d1_Index)[1]*dim(d1_Index)[2], -0.01, 0.01), dim(d1_Index)[1]), Knum1)$classification
ari_s = adjustedRandIndex(mclust1, mclust2)
###
#used_rows = sample(dim(d1_Signal_scale)[1], 5000)
mclust1 = Mclust(cbind(d1_Index+matrix(runif(dim(d1_Index)[1]*dim(d1_Index)[2], -0.01, 0.01), dim(d1_Index)[1]), d1_Signal_scale+matrix(runif(dim(d1_Index)[1]*dim(d1_Index)[2], -0.01, 0.01), dim(d1_Index)[1]))[used_rows,], Knum1)$classification
mclust2 = Mclust(cbind(d2_Index+matrix(runif(dim(d1_Index)[1]*dim(d1_Index)[2], -0.01, 0.01), dim(d1_Index)[1]), d2_Signal_scale+matrix(runif(dim(d1_Index)[1]*dim(d1_Index)[2], -0.01, 0.01), dim(d1_Index)[1]))[used_rows,], Knum1)$classification
ari_bs = adjustedRandIndex(mclust1, mclust2)
ari20_Mclust = rbind(ari20_Mclust, c(ari_b, ari_s, ari_bs))
print(summary(ari20_Mclust))
}
proc.time() - ptm

### 
pdf('ari20.boxplot.snapshot_sig.pdf')
par(mar = c(10, 4, 4, 4))
ari_mat = cbind(cbind(as.numeric(ari20_sp)), ari20_hclust, ari20_km)
colnames(ari_mat) = c('Snapshot', 'Hclust_Index', 'Hclust_Sig', 'Hclust_Both', 'KM_Index', 'KM_Sig', 'KM_Both')
ari_mat_plot = ari_mat[,c(2,5, 3,6, 4,7, 1)]
ari_mat_plot = ari_mat[,c(2:7, 1)]
boxplot(ari_mat_plot, las=2, cex.axis=1.5)
dev.off()









### km raw sig
d1_Index = as.data.frame(fread('hg38_outputs/hg38_chrAll_analysis_rep1/snapshot_test_run_rep1.index.matrix.txt'))[,-c(1:4)]
d2_Index = as.data.frame(fread('hg38_outputs/hg38_chrAll_analysis_rep2/snapshot_test_run_rep2.index.matrix.txt'))[,-c(1:4)]

d1_Signal = as.data.frame(fread('hg38_outputs/hg38_chrAll_analysis_rep1/snapshot_test_run_rep1.signal.matrix.txt'))[,-c(1:4)]
d2_Signal = as.data.frame(fread('hg38_outputs/hg38_chrAll_analysis_rep2/snapshot_test_run_rep2.signal.matrix.txt'))[,-c(1:4)]

set.seed(2019)
Knum1 = length(unique(d1_sp0[used_id, 16]))
Knum2 = length(unique(d2_sp0[used_id, 16]))

mean_ref = mean(as.matrix(d2_Signal))
d1_Signal_scale = apply(d1_Signal, 2, function(x) x/mean(x)*mean_ref )
d2_Signal_scale = apply(d2_Signal, 2, function(x) x/mean(x)*mean_ref )

d1_Signal_scale = scale(d1_Signal)
d2_Signal_scale = scale(d2_Signal)

ari20_km_raw = c()
for (i in 1:30){
	print(i)
#used_rows = sample(dim(d1_Signal)[1], dim(d1_Signal)[1])
used_rows = sample(dim(d1_Signal)[1], 10000)
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
ari20_km_raw = rbind(ari20_km_raw, c(ari_b, ari_s, ari_bs))
print(summary(ari20_km_raw))
}



set.seed(2019)
Knum1 = length(unique(d1_sp0[used_id, 16]))
Knum2 = length(unique(d2_sp0[used_id, 16]))
ari20_hclust_raw = c()
for (i in 1:30){
	print(i)
used_rows = sample(dim(d1_Index)[1], 10000)
hclust1 = cutree(hclust(dist(d1_Index[used_rows,])), k=Knum1)
hclust2 = cutree(hclust(dist(d2_Index[used_rows,])), k=Knum2)
ari_b = adjustedRandIndex(hclust1, hclust2)
###
#used_rows = sample(dim(d1_Signal_scale)[1], 5000)
hclust1 = cutree(hclust(dist(d1_Signal_scale[used_rows,])), k=Knum1)
hclust2 = cutree(hclust(dist(d2_Signal_scale[used_rows,])), k=Knum2)
ari_s = adjustedRandIndex(hclust1, hclust2)
###
#used_rows = sample(dim(d1_Signal_scale)[1], 5000)
hclust1 = cutree(hclust(dist(cbind(d1_Index, d1_Signal_scale)[used_rows,])), k=Knum1)
hclust2 = cutree(hclust(dist(cbind(d2_Index, d2_Signal_scale)[used_rows,])), k=Knum2)
ari_bs = adjustedRandIndex(hclust1, hclust2)
ari20_hclust_raw = rbind(ari20_hclust_raw, c(ari_b, ari_s, ari_bs))
print(summary(ari20_hclust_raw))
}





pdf('ari20.boxplot.1.raw.pdf')
boxplot(cbind(cbind(as.numeric(ari20_sp)),ari20_hclust_raw, ari20_km_raw))
dev.off()



pdf('ari100.boxplot.1.pdf')
boxplot(cbind(cbind(as.numeric(ari100_sp)), ari100_hclust, ari100_km))
dev.off()




