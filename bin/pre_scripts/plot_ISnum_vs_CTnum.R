setwd('/gpfs/group/yzz2/default/legacy/group/projects/vision/snapshot18_reproduce_0_16lim_all_NB88')

d=read.table('atac_20cell.index.matrix.txt', header=F, sep='\t')
dis = d[,-c(1:4)]
set.seed(2019)
#dis = dis[sample(dim(dis)[1], 100000),]



IS_count_vec = c(2)
IS_min_vec = c(min(table(dis[,1])))
for (i in c(2:dim(dis)[2])){
IS_mat = dis[,1:i]
IS_table = apply(IS_mat, 1, function(x) paste(x, collapse='_'))
IS_count = length(unique(IS_table))
IS_min = median(table(IS_table))
print(IS_count)
IS_count_vec = c(IS_count_vec, IS_count)
IS_min_vec = c(IS_min_vec, IS_min)
}

### get background
num_ccREs = dim(dis)[1]
set.seed(2019)
IS_count_mat_shuf = c()

for (j in c(1:5)){
print(j)
dis_shuffle = c()
for (i in c(1:dim(dis)[2])){
dis_shuffle = cbind(dis_shuffle, dis[sample(num_ccREs),i])
}
IS_count_vec_shuf = c(2)
IS_min_vec_shuf = c(min(table(dis[,1])))
for (i in c(2:dim(dis)[2])){
IS_mat = dis_shuffle[,1:i]
IS_table = apply(IS_mat, 1, function(x) paste(x, collapse='_'))
IS_count = length(unique(IS_table))
IS_min = median(table(IS_table))
print(IS_count)
IS_count_vec_shuf = c(IS_count_vec_shuf, IS_count)
IS_min_vec_shuf = c(IS_min_vec_shuf, IS_min)
}
IS_count_mat_shuf = cbind(IS_count_mat_shuf, IS_count_vec_shuf)
IS_min_mat_shuf = cbind(IS_min_mat_shuf, IS_min_vec_shuf)
}




pdf('IS_n_vs_CT_n.pdf', width=4.5, height=4.5)
par(mfrow=c(1,1))
plot(1:length(IS_count_vec), IS_count_vec, log='', ylim=c(1, max(cbind(IS_count_vec, IS_count_mat_shuf))), xlab='', ylab='')
lines(1:length(IS_count_vec), IS_count_vec, lwd=1.5)
points(1:length(IS_count_vec), rowMeans(IS_count_mat_shuf), col='gray')
lines(1:length(IS_count_vec), rowMeans(IS_count_mat_shuf), lwd=1.5, col='gray')
CI_h = rowMeans(IS_count_mat_shuf) + apply(IS_count_mat_shuf, 1, sd)/sqrt(dim(IS_count_mat_shuf)[2]) * 1.960
CI_l = rowMeans(IS_count_mat_shuf) - apply(IS_count_mat_shuf, 1, sd)/sqrt(dim(IS_count_mat_shuf)[2]) * 1.960
#lines(1:length(IS_count_vec), CI_h, lwd=1, col='gray')
#lines(1:length(IS_count_vec), CI_l, lwd=1, col='gray')
#plot(1:length(IS_count_vec), IS_min_vec, log='')
#lines(1:length(IS_count_vec), IS_min_vec, lwd=1.5)
#for (i in 1:50){
#points(1:length(IS_count_vec), IS_min_mat_shuf[,i], col='gray')
#lines(1:length(IS_count_vec), IS_min_mat_shuf[,i], lwd=1.5, col='gray')
#}
dev.off()

IS_table = apply(IS_mat, 1, function(x) paste(x, collapse='_'))
IS_ccRE_count = table(IS_table)

library(ggplot2)
pdf('ccRE_count_hist.10k.pdf', width=3.5, height=3.)
#mydata_hist = hist((IS_ccRE_count), breaks=50, xlab='', ylab='', main='')
#plot(mydata_hist$count, log="y", type='h', lwd=10, lend=2)
#box()
IS_ccRE_count = as.data.frame(IS_ccRE_count)
ggplot(as.data.frame(IS_ccRE_count), aes(x = Freq)) + geom_histogram() + scale_x_log10() + theme_bw() + theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()





