library(ggplot2)
### get parameters
args = commandArgs(trailingOnly=TRUE)
index_matrix_signal_inputfile = args[1]
input_list = args[2]
outfile = args[3]
master_pk = args[4]

### read index set matrix
read_color = function(x){
	rgb_color_int = as.numeric(unlist(strsplit(x, ',')))
	rgb_color = rgb(rgb_color_int[1],rgb_color_int[2],rgb_color_int[3],max=255)
	return(rgb_color)
}


#################################################### 
############ read input files
####################################################
### read signal matrix file
print('read signal matrix file')
signal_matrix_od = as.matrix(read.table(index_matrix_signal_inputfile, header=FALSE))
peak_bed_od = read.table(master_pk, header=FALSE, sep='\t')

### extract signal matrix without info
signal_matrix = signal_matrix_od[ , c(3:dim(signal_matrix_od)[2]) ]
signal_matrix_bed_info = signal_matrix_od[ , c(1,2)]
signal_matrix_bed_info = t( apply(signal_matrix_bed_info, 1, function(x) c(unlist(strsplit(x[1], '_')),x[2])) )
### convert to numeric matrix
class(signal_matrix) = 'numeric'
###### read colnames file
colname_file = read.table(input_list, header=F)
colname = colname_file[,1]
colnames(signal_matrix) = colname

### index set
print('get uniq index set id')
index_set_id = signal_matrix_od[,2]
index_set_id_uniq = unique(index_set_id)
### sort index
index_set_id_uniq_sort = sort(index_set_id_uniq)


for (k in c(1:length(index_set_id_uniq))){
	print(paste('index set', toString(index_set_id_uniq[k])))
	signal_matrix_tmp = signal_matrix[index_set_id==index_set_id_uniq[k],]
	signal_matrix_bed_info_tmp = signal_matrix_bed_info[index_set_id==index_set_id_uniq[k],]
	signal_table_df = c()
	for (i in c(1:dim(signal_matrix_tmp)[2])){
		signal_table_id = rep(colname[i],dim(signal_matrix_tmp)[1])
		signal_table_tmp = data.frame(signal_table_id, signal_matrix_tmp[,i])
		signal_table_df = rbind(signal_table_df, signal_table_tmp)
	}
	### save bed files
	bed_pk_IS_k = peak_bed_od[as.numeric(signal_matrix_bed_info_tmp[,1]),]
	signal_matrix_bed_info_tmp_bed = cbind(bed_pk_IS_k, signal_matrix_bed_info_tmp[,2])
	write.table(signal_matrix_bed_info_tmp_bed, file = paste(toString(k-1), '.', toString(index_set_id_uniq[k]), '.index_set.bed', sep=''), quote = FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
	### save figure
	colnames(signal_table_df) = c('celltype', 'signal')
	### keep the x-axis order
	signal_table_df$celltype = factor(signal_table_df$celltype, levels=unique(signal_table_df$celltype))
	p = ggplot(signal_table_df, aes(factor(celltype), signal), )
	p + geom_violin(scale = 'width', aes(fill = 'red')) + scale_x_discrete(labels=colname) +  stat_summary(fun.y=mean, geom='point', shape=23, size=4) + theme(legend.position="none", axis.text.x = element_text(angle=270)) + xlab(toString(index_set_id_uniq[k]))
	ggsave(paste(toString(k-1), '.', toString(index_set_id_uniq[k]), '.', outfile, sep=''), width = dim(signal_matrix)[2]+1, height = dim(signal_matrix)[2]+1)
}

#time Rscript ~/group/software/CD_viewer/bin/plot_ct_indexset_violin.R atac_20cell.bed.signal.matrix.txt.indexed.sort.txt signal_list.txt violin.pdf
