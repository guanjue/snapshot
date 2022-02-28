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
padjvec = p.adjust(pvec, method='bonferroni')
print(summary(as.matrix(padjvec)))
### get NB count thresh
counts_pfdr = cbind(index_count, padjvec)
if (sum(padjvec<1e-2)>0){
	print('sum(padjvec<0.01)>0')
	NB_count_thresh = min(counts_pfdr[padjvec<1e-2,1])
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

