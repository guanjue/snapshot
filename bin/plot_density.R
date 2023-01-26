### get parameters
args = commandArgs(trailingOnly=TRUE)
index_matrix = args[1]
count_threshold = as.numeric(args[2])
top_lim = 20

index0 = read.table(index_matrix, header = F)

index = apply(index0, 1, function(x) paste(x[5:dim(index0)[2]], collapse='_'))
index_count = table(index)
print(summary(as.matrix(index_count)))
print(length(index_count))


if (length(index_count)>top_lim){
	### select data without top X quantile
	top = 0.95
	top = 0.90
	mean_95 = mean(index_count[index_count<=quantile(index_count, top)])
	var_95 = var(index_count[index_count<=quantile(index_count, top)])
	sd_95 = sd(index_count[index_count<=quantile(index_count, top)])
	print(mean_95)
	print(var_95)
	print(summary(as.matrix(index_count[index_count<quantile(index_count, top)])))

	### get NB size and prob
	prob=mean_95/var_95
	if (prob < 0.001){
		prob = 0.001
	} else if (prob > 0.999){
		prob = 0.999
	}
	size = prob*mean_95 / (1-prob)

	### get NB p-val
	pvec = pnbinom(index_count, size=size, prob=prob, lower.tail = FALSE)
	print(summary(as.matrix(pvec)))
	### FDR pval
	FDR_thresh = 1e-2
	#FDR_thresh = 1e-2
	padjvec = p.adjust(pvec, method='fdr')
	print(summary(as.matrix(padjvec)))
	### get NB count thresh
	counts_pfdr = cbind(index_count, padjvec)
	if (sum(padjvec<FDR_thresh)>0){
		print('sum(padjvec<0.01)>0')
		NB_count_thresh = min(counts_pfdr[padjvec<FDR_thresh,1])
	} else {
		print('user provide')
		NB_count_thresh = count_threshold
		print('NB model fail, use user provide count_threshold')
	}
	#if (sum(index_count>NB_count_thresh)<top_lim){
	#	print('use top N lim')
	#	NB_count_thresh = index_count[order(-index_count)][top_lim]
	#} else if (sum(index_count>NB_count_thresh)>1000){
	if (sum(index_count>NB_count_thresh)>1000){
		print('use top 1000 lim')
		NB_count_thresh = index_count[order(-index_count)][1000]
	}
} else{
	NB_count_thresh = 3
}

### make sure the threshold is greater than the dimention
if (NB_count_thresh<(dim(index0)[2]-4)){
	NB_count_thresh = (dim(index0)[2]-4+1)
}


###### plot density
png('density_index_count.png', width = 8, height = 8, units = 'in', res = 300)
plot(density(log2(index_count)), main='density plot of index counts (log2 scale)')
abline(v=log2(count_threshold), col='red', lwd=1.5, lty=2)
abline(v=log2(NB_count_thresh), col='blue', lwd=1.5, lty=2)
dev.off()

print('initial IS number: ')
print(sum(index_count>NB_count_thresh))
write.table(as.matrix(NB_count_thresh), paste(index_matrix, '.NB_count_thresh.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')

print('Dif FDR threshold Index-Sets count:')
print(summary(as.matrix(padjvec)))

FDR_thresh_vec = c(0.5, 0.4, 0.3, 0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5)
IS_n = c()
for (fdri in FDR_thresh_vec){
	IS_n = c(IS_n, sum(as.numeric(padjvec)<fdri))
	print(c(fdri, sum(as.numeric(padjvec)<fdri)) )
}

pdf('IS_num_vs_FDRthresh.pdf', width=4, height=4)
plot(FDR_thresh_vec, IS_n, log='x', xlab='FDR adjusted p-value threshold', ylab='Index-Sets Number', cex.axis=1, ylim=c(0, max(IS_n)))
lines(FDR_thresh_vec, IS_n)
abline(v=0.01, lty=2)
dev.off()


MetaIS_num_vs_FDR = read.table('MetaIS_num_vs_FDR.txt', header=F)
pdf('MetaIS_num_vs_FDRthresh.pdf', width=4, height=4)
plot(MetaIS_num_vs_FDR[,2], MetaIS_num_vs_FDR[,1], log='x', xlab='FDR adjusted p-value threshold', ylab='Meta-Index-Sets Number', cex.axis=1, ylim=c(0, max(MetaIS_num_vs_FDR[,1])))
lines(MetaIS_num_vs_FDR[,2], MetaIS_num_vs_FDR[,1])
abline(v=0.01, lty=2)
dev.off()




