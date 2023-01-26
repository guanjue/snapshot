library(networkD3)
library(igraph)

####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
signal_matrix_file = args[1]
cd_tree = args[2]
signal_range_color_file = args[3]
signal_input_list = args[4]
signal_matrix_start_col = as.numeric(args[5])

#################################################### 
############ read input files
####################################################
### read signal matrix file
signal_matrix_od = as.matrix(read.table(signal_matrix_file, header=FALSE, comment.char='~'))
### extract signal matrix without info
signal_matrix = signal_matrix_od[ , signal_matrix_start_col:dim(signal_matrix_od)[2] ]
### get index_set name
index_set_name = signal_matrix_od[,1]

###### read colnames file
colname_file = read.table(signal_input_list, header=F)
colname = colname_file[,1]

### read cell development tree file
tree = read.table(cd_tree, header = F, sep=',')
tree.df = as.data.frame(tree)
colnames(tree.df) = c('Node.1', 'Node.2')

### get color list
signal_matrix_color = signal_matrix#matrix("#FFFFFF", nrow=dim(signal_matrix)[1], ncol=dim(signal_matrix)[2])


### plot trees
for (i in seq(1,dim(signal_matrix_color)[1])){
	### get color vector from color matrix
	value_col = signal_matrix_color[i,]
	### get tree
	tree.igraph = graph.data.frame(tree.df, directed=TRUE)
	tree_names = V(tree.igraph)$name
	V(tree.igraph)$name = tree_names#rep('', length(tree_names))
	### sort colnames by tree nodes id
	match_id = match(tree_names, colname)
	V(tree.igraph)$color = value_col[match_id]
	V(tree.igraph)$size = 25
	pdf(paste(toString(i-1), '.', signal_input_list, index_set_name[i], '.tree.pdf', sep = ''), width = 12, height = 12)
	plot(tree.igraph, layout = layout_as_tree(tree.igraph, root=c(1)))
	dev.off()
}

