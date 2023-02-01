library(data.table)
library(mclust)
library(pheatmap)
library(dplyr)


file_list = c('snapshot_test_run_merge.index_binary_mat.173thesh.txt', 'snapshot_test_run_merge.index_binary_mat.500thesh.txt', 'snapshot_test_run_merge.index_binary_mat.1000thesh.txt')
for (input_binary_map_mat_file in file_list){
###
d = read.table(input_binary_map_mat_file, header=F)

dsi = as.matrix(d[,-1])
dsi[dim(dsi)[1],] = rep('0.5', dim(dsi)[2])
dsi = apply(dsi, 2, as.numeric)

colnames(dsi) = read.table('peak_signal_state_list.merge.txt', header=F)[,1]
###
pdf(paste0(input_binary_map_mat_file, '.pheatmap.pdf'), width=4)
breaksList = seq(0, (max(dsi)), by = 0.01)
my_colorbar=colorRampPalette(c('white', 'black'))(n = length(breaksList))
pheatmap(dsi, cluster_col=F, cluster_row=F,cex=1.5, color=my_colorbar, breaks = breaksList, show_colnames=T, legend=F)
dev.off()
}



