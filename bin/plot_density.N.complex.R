### get parameters
args = commandArgs(trailingOnly=TRUE)
index_matrix = args[1]
count_threshold = as.numeric(args[2])
use_user_define_thresh = args[3]
top_lim = 20

if (use_user_define_thresh == 'F'){
index0 = read.table(index_matrix, header = F)

index = apply(index0, 1, function(x) paste(x[5:dim(index0)[2]], collapse='_'))
index_count = table(index)
print(summary(as.matrix(index_count)))
print(length(index_count))
write.table(cbind(index_count[order(index_count)]), paste(index_matrix, '.index_count.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')


padj_thresh = 0.1
count_threshold = 2:200
count_ISn_prod = c()
index_count_OD = index_count


for (count_threshold_i in count_threshold){
	print(count_threshold_i)
	index_count = index_count_OD
	if (length(index_count)>top_lim){
		IS_count_pre = 0
		index_count = index_count[index_count>=count_threshold_i]
		for (iter_i in 1:100){
			print(c(iter_i, length(index_count_OD)-IS_count_pre))
			### select data without top X quantile
			top = 1
			mean_95 = mean(index_count)
			sd_95 = sd(index_count)

			### get Norm p-val
			pvec = pnorm(index_count, mean=mean_95, sd=sd_95, lower.tail = FALSE)
			padjvec = p.adjust(pvec, method='fdr')
			### FDR pval
			###
			index_count = index_count[padjvec>padj_thresh]
			if (sum(padjvec>=padj_thresh)==IS_count_pre){break}
			IS_count_pre = sum(padjvec>=padj_thresh)
		}

		print(summary(as.matrix(padjvec)))
		### get count thresh
		pvec_final = pnorm(index_count_OD, mean=mean_95, sd=sd_95, lower.tail = FALSE)
		padjvec_final = p.adjust(pvec_final, method='fdr')
		counts_pfdr = cbind(index_count_OD, padjvec_final)
		print('sum(padjvec<0.01)>0')
		NB_count_thresh = min(counts_pfdr[padjvec_final<padj_thresh,1])
	} else{
		NB_count_thresh = 3
	}
	if (!is.finite(NB_count_thresh)){NB_count_thresh=10000}
	count_ISn_prod = rbind(count_ISn_prod, c(count_threshold_i, sum(index_count_OD>=NB_count_thresh), NB_count_thresh, NB_count_thresh*sum(index_count_OD>=NB_count_thresh)))
}

### ignore last N
count_ISn_prod_used = count_ISn_prod[count_ISn_prod[,4]!=count_ISn_prod[dim(count_ISn_prod)[1],4],]
### get z prod
count_ISn_prod_used_zprod = (count_ISn_prod[,2])/sd(count_ISn_prod_used[,2]) * (count_ISn_prod[,3])/sd(count_ISn_prod_used[,3])
count_ISn_prod_used_zprod[count_ISn_prod[,4]==count_ISn_prod[dim(count_ISn_prod)[1],4]] = 0
count_ISn_prod = cbind(count_ISn_prod, count_ISn_prod_used_zprod)

colnames(count_ISn_prod) = c('count_min_thresh', 'ISn', 'count_min_thresh_est', 'prod', 'zprod')
write.table(as.matrix(count_ISn_prod), paste(index_matrix, '.count_thresh.mat.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')


###### plot density
png('density_index_count.png', width = 8, height = 8, units = 'in', res = 300)
plot(density((index_count_OD)), main='density plot of index counts (log2 scale)', log='x')
abline(v=(count_threshold), col='red', lwd=1.5, lty=2)
abline(v=(count_ISn_prod[which.max(count_ISn_prod[,5]),3]), col='blue', lwd=1.5, lty=2)
abline(v=quantile(index_count_OD, top), col='black', lwd=1.5, lty=2)
dev.off()

print('initial IS number & NB_count_thresh')
print(c(sum(index_count_OD>=count_ISn_prod[which.max(count_ISn_prod[,5]),3]), count_ISn_prod[which.max(count_ISn_prod[,5]),3]))
write.table(as.matrix(count_ISn_prod[which.max(count_ISn_prod[,5]),3]), paste(index_matrix, '.NB_count_thresh.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')
} else{
	write.table(as.matrix(count_threshold), paste(index_matrix, '.NB_count_thresh.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')
}
