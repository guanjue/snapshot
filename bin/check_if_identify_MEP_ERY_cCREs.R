library(data.table)
library(mclust)
library(pheatmap)
library(dplyr)


### input files
#cut -f1,2,3,4 snapshot_test_run_merge.index.matrix.txt > snapshot_test_run_merge.bedinfo.bed
bed_file = 'snapshot_test_run_merge.bedinfo.bed'
input_metaIS_mat_file = 'snapshot_test_run_merge.metaISid.mat.txt'
peak_list_file = 'peak_signal_state_list.merge.txt'

### read data
bed_info = as.data.frame(fread(bed_file, header=F))[1:83701,]
#bed_info = as.data.frame(fread(bed_file, header=F))
dmerge_sp0 = as.data.frame(fread(input_metaIS_mat_file, header=F))
ct_info = as.data.frame(fread(peak_list_file, header=F))[,1]

### write meta-IS bed files
output_dir = "hg38_outputs/hg38_chrAll_analysis_merge_MetaIS_bed/"
dir.create(output_dir)
for (i in 1:max(dmerge_sp0[,16])){
	bed_info_i = bed_info[dmerge_sp0[,16]==i,]
	write.table(bed_info_i, paste0(output_dir, 'MetaIS_', i, '.bed'), quote=F, col.names=F, row.names=F, sep='\t')
}


### get dmerge_sp0 merge mat
meta_IS_vec = dmerge_sp0[,16]
dmerge_sp0_sigmat = dmerge_sp0[,-c(1,2,16)]
meta_IS_coluster_mean_mat0 = bind_rows(lapply(split(dmerge_sp0_sigmat, meta_IS_vec), colMeans))
colnames(meta_IS_coluster_mean_mat0) = ct_info
###
pdf(paste0(input_metaIS_mat_file, '.pheatmap.pdf'))
breaksList = seq(0, (max(meta_IS_coluster_mean_mat0)), by = 0.01)
my_colorbar=colorRampPalette(c('white', 'red'))(n = length(breaksList))
pheatmap(meta_IS_coluster_mean_mat0, cluster_col=F, cluster_row=F,cex=1.5, color=my_colorbar, breaks = breaksList)
dev.off()

### get cCRE per metaIS
cCRE_num_per_metaIS = c()
for (i in 1:max(meta_IS_vec)){
	cCRE_num_per_metaIS = c(cCRE_num_per_metaIS, sum(meta_IS_vec==i))
}
pdf(paste0(input_metaIS_mat_file, '.cCRE_num_per_metaIS.barplot.pdf'), height=3)
barplot(cCRE_num_per_metaIS, log='', ylim=c(0, 3000))
box()
dev.off()

pdf(paste0(input_metaIS_mat_file, '.cCRE_num_per_metaIS.barplot.all.pdf'), height=3)
barplot(cCRE_num_per_metaIS, log='y')
box()
dev.off()

### kmeans
### read raw sig mat
#draw_sig = as.data.frame(fread(input_sig_file_raw, header=F))[,-c(1:4)]
#draw_sig_scale = scale(draw_sig)
#draw_sig_ref = mean(as.matrix(draw_sig))
#draw_sig_scale = apply(draw_sig, 2, function(x) x/mean(x)*draw_sig_ref )
#draw_sig_scale = dmerge_sp0[,-c(1,2,16)]
draw_sig_scale = dmerge_sp0[,-c(1,2,16)]
draw_sig_scale[dmerge_sp0[,-c(1,2,16)]>16] = 16


### kmeans
set.seed(2019)
KM_coluster_mean_mat0_all = c()
km_num = 100
output_dir = "hg38_outputs/hg38_chrAll_analysis_merge_other_cluster_KM/"
dir.create(output_dir)
for (ki in 1:km_num){
	print(ki)
draw_sig_sigmat_km = kmeans(draw_sig_scale, max(meta_IS_vec))$cluster
### get KM mean mat
KM_coluster_mean_mat0 = bind_rows(lapply(split(as.data.frame(draw_sig_scale), draw_sig_sigmat_km), colMeans))
colnames(KM_coluster_mean_mat0) = ct_info
#rownames(KM_coluster_mean_mat0) = paste0('km:',ki,':',1: max(meta_IS_vec))
rownames(KM_coluster_mean_mat0) = (1: max(meta_IS_vec))
###
KM_coluster_mean_mat0_all = rbind(KM_coluster_mean_mat0_all, KM_coluster_mean_mat0)
###
pdf(paste0(output_dir, input_metaIS_mat_file, '.', ki, '.pheatmap.KM_and_MetaIS.pdf'))
#breaksList = seq(-max(abs(KM_coluster_mean_mat0)), max(abs(KM_coluster_mean_mat0)), by = 0.01)
breaksList = seq(-max(abs(KM_coluster_mean_mat0)), max(abs(KM_coluster_mean_mat0)), by = 0.01)
my_colorbar=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList))
pheatmap(KM_coluster_mean_mat0, cluster_col=F, cluster_row=T,cex=1.5, color=my_colorbar, breaks = breaksList)
dev.off()
}






pdf(paste0(output_dir, input_metaIS_mat_file, '.', 'all', '.pheatmap.KM_and_MetaIS.pdf'), height=15)
#breaksList = seq(-max(abs(KM_coluster_mean_mat0)), max(abs(KM_coluster_mean_mat0)), by = 0.01)
breaksList = seq(0, max(abs((max((KM_coluster_mean_mat0_all))))), by = 0.01)
my_colorbar=colorRampPalette(c('white', 'red'))(n = length(breaksList))
pheatmap((KM_coluster_mean_mat0_all), cluster_col=F, cluster_row=T,cex=1.5, color=my_colorbar, breaks = breaksList, show_rownames=F)
dev.off()

library(lsa)

KM_dist_to_sp = matrix(0, nrow=dim(KM_coluster_mean_mat0_all)[1], ncol=dim(meta_IS_coluster_mean_mat0)[1])

for (i in 1:dim(KM_coluster_mean_mat0_all)[1]){
	for (j in 1:dim(meta_IS_coluster_mean_mat0)[1]){
		#KM_dist_to_sp[i,j] = dist(rbind(KM_coluster_mean_mat0_all[i,], meta_IS_coluster_mean_mat0[j,]))
		KM_dist_to_sp[i,j] = 1 - cosine(as.numeric(KM_coluster_mean_mat0_all[i,]), as.numeric(meta_IS_coluster_mean_mat0[j,]))
	}
}

KM_dist_to_sp_which_min = apply(KM_dist_to_sp,1,which.min)

### which min KM
which_min_mat_KM = c()
for (j in 1:km_num){
	KM_dist_to_sp_which_min_j = KM_dist_to_sp_which_min[(1:max(meta_IS_vec))+max(meta_IS_vec)*(j-1)]
	which_min_mat_KM_j = c()
	for (i in 1:1:max(meta_IS_vec)){
		which_min_mat_KM_j = c(which_min_mat_KM_j, sum(KM_dist_to_sp_which_min_j==i)!=0)*1
		#which_min_mat_KM_j = c(which_min_mat_KM_j, sum(KM_dist_to_sp_which_min_j==i))*1
	}
	which_min_mat_KM = rbind(which_min_mat_KM, which_min_mat_KM_j)
}
which_min_mat_KM_sum = apply(which_min_mat_KM,2,sum)


### Hclust
set.seed(2019)
draw_sig_scale = dmerge_sp0[,-c(1,2,16)]
draw_sig_scale[dmerge_sp0[,-c(1,2,16)]>16] = 16

Hclust_coluster_mean_mat0_all = c()
km_num = 30
output_dir = "hg38_outputs/hg38_chrAll_analysis_merge_other_cluster_hclust/"
dir.create(output_dir)
for (ki in 1:km_num){
	print(ki)
used_rows = sample(dim(draw_sig_scale)[1], 20000)
draw_sig_sigmat_km = cutree(hclust(dist(draw_sig_scale[used_rows,])), k=max(meta_IS_vec))
### get KM mean mat
KM_coluster_mean_mat0 = bind_rows(lapply(split(as.data.frame(draw_sig_scale[used_rows,]), draw_sig_sigmat_km), colMeans))
colnames(KM_coluster_mean_mat0) = ct_info
#rownames(KM_coluster_mean_mat0) = paste0('km:',ki,':',1: max(meta_IS_vec))
rownames(KM_coluster_mean_mat0) = (1: max(meta_IS_vec))
###
Hclust_coluster_mean_mat0_all = rbind(Hclust_coluster_mean_mat0_all, KM_coluster_mean_mat0)
###
pdf(paste0(output_dir, input_metaIS_mat_file, '.', ki, '.pheatmap.KM_and_MetaIS.pdf'))
#breaksList = seq(-max(abs(KM_coluster_mean_mat0)), max(abs(KM_coluster_mean_mat0)), by = 0.01)
breaksList = seq(-max(abs(KM_coluster_mean_mat0)), max(abs(KM_coluster_mean_mat0)), by = 0.01)
my_colorbar=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList))
pheatmap(KM_coluster_mean_mat0, cluster_col=F, cluster_row=T,cex=1.5, color=my_colorbar, breaks = breaksList)
dev.off()
}

pdf(paste0(output_dir, input_metaIS_mat_file, '.', 'all', '.pheatmap.KM_and_MetaIS.pdf'), height=15)
#breaksList = seq(-max(abs(KM_coluster_mean_mat0)), max(abs(KM_coluster_mean_mat0)), by = 0.01)
breaksList = seq(0, max(abs((max((Hclust_coluster_mean_mat0_all))))), by = 0.01)
my_colorbar=colorRampPalette(c('white', 'red'))(n = length(breaksList))
pheatmap((Hclust_coluster_mean_mat0_all), cluster_col=F, cluster_row=T,cex=1.5, color=my_colorbar, breaks = breaksList, show_rownames=F)
dev.off()


### check which min 
Hclust_dist_to_sp = matrix(0, nrow=dim(Hclust_coluster_mean_mat0_all)[1], ncol=dim(meta_IS_coluster_mean_mat0)[1])

for (i in 1:dim(Hclust_coluster_mean_mat0_all)[1]){
	for (j in 1:dim(meta_IS_coluster_mean_mat0)[1]){
		Hclust_dist_to_sp[i,j] = dist(rbind(Hclust_coluster_mean_mat0_all[i,], meta_IS_coluster_mean_mat0[j,]))
	}
}

Hclust_dist_to_sp_which_min = apply(Hclust_dist_to_sp,1,which.min)

### which min Hclust
which_min_mat_Hclust = c()
for (j in 1:km_num){
	Hclust_dist_to_sp_which_min_j = Hclust_dist_to_sp_which_min[(1:max(meta_IS_vec))+max(meta_IS_vec)*(j-1)]
which_min_mat_Hclust_j = c()
for (i in 1:1:max(meta_IS_vec)){
	which_min_mat_Hclust_j = c(which_min_mat_Hclust_j, sum(Hclust_dist_to_sp_which_min_j==i)!=0)*1
}
which_min_mat_Hclust = rbind(which_min_mat_Hclust, which_min_mat_Hclust_j)
}
which_min_mat_Hclust_sum = apply(which_min_mat_Hclust,2,sum)


### get which.min mat
which_min_mat = c()
for (i in 1:max(meta_IS_vec)){
which_min_mat = rbind(which_min_mat, c(sum(KM_dist_to_sp_which_min==i), sum(Hclust_dist_to_sp_which_min ==i)))
}

rownames(which_min_mat) = 1:max(meta_IS_vec)
colnames(which_min_mat) = c('Kmeans', 'Hclust')

pdf('KM_Hclust_which_min_dist.pdf', width=3)
which_min_mat_p = which_min_mat/km_num
which_min_mat_p = cbind(sqrt(which_min_mat[,1]*which_min_mat[,2]))
#which_min_mat_p[which_min_mat_p>1] = 1
breaksList = seq(0, max(which_min_mat_p), by = 0.01)
my_colorbar=colorRampPalette(c('white', 'red'))(n = length(breaksList))
pheatmap(which_min_mat_p, cluster_col=F, cluster_row=F,cex=1.5, color=my_colorbar, breaks = breaksList, show_rownames=T)
dev.off()

pdf('KM_which_min_dist.barplot.100.pdf', height=3)
barplot(which_min_mat_KM_sum, log='')
box()
dev.off()


pdf('Hclust_which_min_dist.barplot.pdf', height=3)
barplot(which_min_mat_Hclust_sum, log='')
box()
dev.off()





