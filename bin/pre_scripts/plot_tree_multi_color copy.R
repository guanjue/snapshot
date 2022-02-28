library(networkD3)
library(igraph)

####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
signal_matrix_file = args[1]
cd_tree = args[2]
signal_range_color_file = args[3]
signal_input_list = args[4]
signal_matrix_start_col = args[5]

#################################################### 
############ read input files
####################################################
### read signal matrix file
signal_matrix_od = as.matrix(read.table(signal_matrix_file, header=FALSE))
### extract signal matrix without info
signal_matrix = signal_matrix_od[ , signal_matrix_start_col:dim(signal_matrix_od)[2] ]
### get index_set name
index_set_name = signal_matrix_od[,1]
### convert to numeric matrix
class(signal_matrix) = 'numeric'

###### read colnames file
colname_file = read.table(signal_input_list, header=F)
colname = colname_file[,2]

### read cell development tree file
tree = read.table(cd_tree, header = F, sep=',')
tree.df = as.data.frame(tree)
colnames(tree.df) = c('Node.1', 'Node.2')

### get color list

signal_matrix_color = signal_matrix
###### read color file
signal_range_color = as.matrix(read.table(signal_range_color_file, header=F))

### replace NA by background signal in background color table
signal_matrix[is.na(signal_matrix)] = as.numeric(signal_range_color[dim(signal_range_color)[1],1])

### add colnames
for (i in seq(1,dim(signal_range_color)[1])){
	### read signal range and color range
	sig_range_high = signal_range_color[i,1]
	class(sig_range_high) = 'numeric'
	sig_range_low = signal_range_color[i,2]
	class(sig_range_low) = 'numeric'
	sig_range_high_col = unlist(strsplit(signal_range_color[i,3], split = ','))
	class(sig_range_high_col) = 'numeric'
	sig_range_high_col = sig_range_high_col / 255
	sig_range_low_col = unlist(strsplit(signal_range_color[i,4], split = ','))
	class(sig_range_low_col) = 'numeric'
	sig_range_low_col = sig_range_low_col / 255
	### identify signal within this signal range of the color
	signal_color_cell = as.logical((signal_matrix<=sig_range_high) * (signal_matrix>sig_range_low))
	col_range = sig_range_high - sig_range_low
	###### get r channel score
	r = sig_range_high_col[1] - ((sig_range_high-signal_matrix) / col_range) * (sig_range_high_col[1] - sig_range_low_col[1])
	### set upper & lower limit for r
	r[r>1] = 1
	r[r<0] = 0
	###### get g channel score
	g = sig_range_high_col[2] - ((sig_range_high-signal_matrix) / col_range) * (sig_range_high_col[2] - sig_range_low_col[2])
	### set upper & lower limit for g
	g[g>1] = 1
	g[g<0] = 0
	###### get b channel score
	b = sig_range_high_col[3] - ((sig_range_high-signal_matrix) / col_range) * (sig_range_high_col[3] - sig_range_low_col[3])
	### set upper & lower limit for b
	b[b>1] = 1
	b[b<0] = 0
	###### creat rgb color matrix
	color_matrix_tmp = NULL
	for (i in seq(1,dim(r)[2])){
		#print(i)
		### convert each r,g,b vector to a rgb matrix column
		cor_col_tmp = rgb( r[,i], g[,i], b[,i] )
		### cbind column to a rgb matrix
		color_matrix_tmp = cbind(color_matrix_tmp, cor_col_tmp)
	}
	###### replace color matrix cell values
	signal_matrix_color[ signal_color_cell ] = color_matrix_tmp[signal_color_cell]
}

### plot trees
for (i in seq(1,dim(signal_matrix_color)[1])){
	### get color vector from color matrix
	value_col = signal_matrix_color[i,]

	### get tree
	tree.igraph = graph.data.frame(tree.df, directed=TRUE)
	tree_names = V(tree.igraph)$name
	V(tree.igraph)$name = rep('', length(tree_names))
	### sort colnames by tree nodes id
	match_id = match(tree_names, colname)
	V(tree.igraph)$color = value_col[match_id]
	V(tree.igraph)$size = 25

	png(paste(toString(i-1), '.', signal_input_list, index_set_name[i], '.tree.png', sep = ''), width = 800, height = 800)
	plot(tree.igraph, layout = layout_as_tree(tree.igraph, root=c(1)))
	dev.off()
}

