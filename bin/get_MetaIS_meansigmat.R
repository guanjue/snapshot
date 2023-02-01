library(data.table)

d = as.data.frame(fread('snapshot_test_run_merge.metaISid.mat.txt'))

d1 = cbind(d[,c(2,16,3:15)])
write.table(d1, 'snapshot_test_run_merge.MetaIS.forviolin.sig.txt', quote=F, col.names=F, row.names=F, sep='\t')

mean_sig_mat = c()

for (i in 1:max(d[,16])){
	mean_sig_mat_i = d[d[,16]==i ,-c(1,2,16)]
	mean_sig_mat = rbind(mean_sig_mat, colMeans(mean_sig_mat_i))
}

write.table(cbind(1:max(d[,16]), mean_sig_mat), 'snapshot_test_run_merge.metaISid.meansigmat.txt', quote=F, col.names=F, row.names=F, sep='\t')


