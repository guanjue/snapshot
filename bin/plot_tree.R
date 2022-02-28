library(networkD3)
library(igraph)

####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
signal_matrix_file = args[1]
cd_tree = args[2]
input_list = args[3]
signal_matrix_start_col = args[4]
high_color = args[5]
low_color = args[6]
log2 = args[7]
smallnum = as.numeric(args[8])

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
colname_file = read.table(input_list, header=F)
colname = colname_file[,1]

### read cell development tree file
tree = read.table(cd_tree, header = F, sep=',')
tree.df = as.data.frame(tree)
colnames(tree.df) = c('Node.1', 'Node.2')

### get color list
color_number = 300
color_range = max(signal_matrix) - min(signal_matrix)

my_colorbar=colorRampPalette(c(low_color, high_color))(n = color_number+1)

### plot trees
for (i in seq(1,dim(signal_matrix)[1])){
	### convert mean signal vector to color vector
	values = signal_matrix[i,]
	### get color list
	values_id = as.integer( (values - min(signal_matrix))/color_range * color_number )
	### value to color (10X)
	value_col = my_colorbar[values_id]

	### get tree
	tree.igraph = graph.data.frame(tree.df, directed=TRUE)
	tree_names = V(tree.igraph)$name
	V(tree.igraph)$name = tree_names#rep('', length(tree_names))
	### sort colnames by tree nodes id
	match_id = match(tree_names, colname)
	V(tree.igraph)$color = value_col[match_id]
	V(tree.igraph)$size = 25

	pdf(paste(toString(i-1), '.', input_list, index_set_name[i], '.tree.pdf', sep = ''), width = 12, height = 12)
	plot(tree.igraph, layout = layout_as_tree(tree.igraph, root=c(1)))
	dev.off()
}

