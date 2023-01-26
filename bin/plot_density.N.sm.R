### get parameters
args = commandArgs(trailingOnly=TRUE)
index_matrix = args[1]
count_threshold = as.numeric(args[2])
use_user_define_thresh = args[3]
top_lim = 20

#index_matrix = 'snapshot_test_run_r2.index.matrix.txt'


index0 = read.table(index_matrix, header = F)

index = apply(index0, 1, function(x) paste(x[5:dim(index0)[2]], collapse='_'))
index_count = table(index)
print(summary(as.matrix(index_count)))
print(length(index_count))

#write.table(cbind(index_count[order(index_count)]), paste(index_matrix, '.index_count.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')
padj_thresh = 0.1

count_threshold = 2:5000#seq(10,200, by=10)
count_ISn_prod = c()
index_count_OD = index_count

for (count_threshold_i in count_threshold){
	#print(count_threshold_i)
	index_count = index_count_OD
	IS_count_pre = 0
	prod = count_threshold_i * sum(index_count>=count_threshold_i)
	count_ISn_prod = rbind(count_ISn_prod, c(count_threshold_i, sum(index_count_OD>=count_threshold_i)))#, NB_count_thresh, NB_count_thresh*sum(index_count_OD>=NB_count_thresh)))
}

### 
count_ISn_prod_1z = (count_ISn_prod[,1] - mean(count_ISn_prod[,1])) / sd(count_ISn_prod[,1]) + count_ISn_prod[,1]
count_ISn_prod_2z = (count_ISn_prod[,2] - mean(count_ISn_prod[,2])) / sd(count_ISn_prod[,2]) + count_ISn_prod[,2]

### smooth
count_ISn_prod_12z = count_ISn_prod_1z * count_ISn_prod_2z
count_ISn_prod_12z_loess = rep(0, length(count_ISn_prod_12z))
which_max_vec = c()
for (fi in seq(0.1, 0.5, by=0.01)){
	xloessi = lowess(1:length(count_ISn_prod_12z), count_ISn_prod_12z, f=fi)
	which_max_vec = c(which_max_vec, which.max(xloessi$y))
	count_ISn_prod_12z_loess = count_ISn_prod_12z_loess + xloessi$y
}


count_threshold_max = count_threshold[round(mean(which_max_vec))]
print('Smooth count thresh:')
print(count_threshold_max)
write.table(as.matrix(count_threshold_max), paste(index_matrix, '.NB_count_thresh.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')

colnames(count_ISn_prod) = c('count_threshold_i', 'ISn')
write.table(as.matrix(cbind(count_ISn_prod_1z, count_ISn_prod_2z)), paste(index_matrix, '.count_thresh.mat.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')


