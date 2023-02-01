library(data.table)
library(mclust)
library(pheatmap)
library(cba)


### read index
d1_Index = as.data.frame(fread('hg38_outputs/hg38_chrAll_analysis_merge/snapshot_test_run_merge.index.matrix.txt'))[,-c(1:4)]
Index_set_binary_OD = apply(as.matrix(d1_Index), 1, function(x) paste0(x, collapse='_'))
Index_set_binary0 = read.table('snapshot_test_run_merge.metaISid.mat.txt', header=F)
###
Index_set_binary1 = apply(Index_set_binary0, 1, function(x) as.character(x[1]) )

### replace X by original Index
Index_set_binary = Index_set_binary1
Index_set_binary[Index_set_binary1=='X_X_X_X_X_X_X_X_X_X_X_X_X'] = Index_set_binary_OD[Index_set_binary1=='X_X_X_X_X_X_X_X_X_X_X_X_X']

### hist

pdf('Index_set_binary.count.hist.pdf')
hist_info = hist(log2(table(Index_set_binary)), breaks=20, plot=F)
plot(hist_info$breaks[-1], hist_info$counts, type='h', cex.axis=2, xlab='', ylab='', log='y', ylim=c(1,1000), xlim=c(0,13))
abline(v=log2(173), col='red', lty=2, lwd=3)
dev.off()

pdf('Index_set_binary_OD.count.hist.pdf')
hist_info = hist(log2(table(Index_set_binary_OD)), breaks=20, plot=F)
plot(hist_info$breaks[-1], hist_info$counts, type='h', cex.axis=2, xlab='', ylab='', log='y', ylim=c(1,1000), xlim=c(0,13))
abline(v=log2(173), col='red', lty=2, lwd=3)
dev.off()


