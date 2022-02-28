### get parameters
args = commandArgs(trailingOnly=TRUE)
index_matrix_ideas_state_inputfile = args[1]
signal_input_list = args[2]
ideas_state_color = args[3]
cREs_IDEASpro_outfile = args[4]

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
signal_matrix_od = as.matrix(read.table(index_matrix_ideas_state_inputfile, header=FALSE))
### extract signal matrix without info
signal_matrix = signal_matrix_od[ , c(3:dim(signal_matrix_od)[2]) ]
### convert to numeric matrix
class(signal_matrix) = 'numeric'
###### read colnames file
colname_file = read.table(signal_input_list, header=F)
colname = colname_file[,2]
colnames(signal_matrix) = colname

### index set
print('get uniq index set id')
index_set_id = signal_matrix_od[,2]
index_set_id_uniq = unique(index_set_id)
### sort index
index_set_id_uniq_sort = sort(index_set_id_uniq)

### get unique elements from the matrix
print('get uniq ideas label id')
ideas_state_matrix_flatten = as.vector(signal_matrix)
ideas_state_matrix_uniq = unique(ideas_state_matrix_flatten)
### sort index
ideas_state_matrix_uniq_sort = sort(ideas_state_matrix_uniq)

### set heatmap colors
print('set heatmap colors')
rgb_col_num = read.table(ideas_state_color,header=F)
rgb_col_num = rgb_col_num[,3]
rgb_col_num = rgb_col_num[c((length(rgb_col_num)-1):1,length(rgb_col_num))]
rgb_col_num = as.matrix(rgb_col_num)
#print(rgb_col_num)
rgb_col=apply(rgb_col_num,1,function(x) read_color(x))
my_colorbar=colorRampPalette(rgb_col)(n = dim(rgb_col_num)[1])


for (k in c(1:length(index_set_id_uniq_sort))){
	print(paste('index set', toString(index_set_id_uniq_sort[k])))
	signal_matrix_tmp = signal_matrix[index_set_id==index_set_id_uniq_sort[k],]
	###### extract counts matrix
	counts_matrix = c()
	for (i in c(1: dim(signal_matrix)[2]) ){
		### extract ith cell type data
		ideas_state_matrix_table_tmp = as.matrix(signal_matrix_tmp[,i])
		table_tmp = c()
		for (j in c( 1: length(ideas_state_matrix_uniq_sort)) ){
			### count the number of cREs have jth IDEAS state
			table_tmp[j] = sum(ideas_state_matrix_table_tmp==(j-1))
		}

		### vector to matrix
		counts_matrix = rbind(counts_matrix, table_tmp)
	}

	### transpose matrix
	counts_matrix_t = t( counts_matrix )
	### add colnames
	colnames(counts_matrix_t) = colnames(signal_matrix)

	### save figure
	png(paste(toString(k-1), '.', toString(index_set_id_uniq_sort[k]), '.', cREs_IDEASpro_outfile, sep=''), dim(signal_matrix)[2]*100+5, dim(signal_matrix)[2]*100+5)
	barplot(counts_matrix_t, col=my_colorbar)
	dev.off()
}

# Rscript plot_ct_IDEASpro_bar.R atac_20cell.bed.ideas.matrix.txt.freq.indexed.sort.txt ideas_list.txt ideas_range_color.txt bar.pdf

