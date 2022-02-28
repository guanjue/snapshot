####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
signal_matrix_file = args[1]
output_filename = args[2]
signal_input_list = args[3]
signal_range_color_file = args[4]
signal_matrix_start_col = args[5]
heatmap_boarder_col = args[6]

####################################################
### use rect to plot heatmaps
color_heatmap = function(color_matrix, outputname, format, border_color){
	format(outputname, width = 1000, height = 1000) ### output name
	par(mar=c(5,0.5,0.5,0.5)) ### set heatmap margins
	colbin_len = 10 ### column bin size
	rowbin_len = 10 ### row bin size
	### row reverse
	color_matrix = color_matrix[nrow(color_matrix):1,]
	### plot areas
	plot(c(0, dim(color_matrix)[2]*colbin_len), c(0, dim(color_matrix)[1]*rowbin_len), xaxt = "n", yaxt = "n", xaxs="i", yaxs="i", type = "n", xlab = "", ylab = "",main = "")
	### add color matrix colname as heatmap colname
	axis(1, c(1 : dim(color_matrix)[2])*colbin_len-0.5*colbin_len, colnames(color_matrix), las = 2, col.axis = "black", tick=FALSE)
	### use for loop to add rectangle with different color
	for (coln in c(1 : dim(color_matrix)[2])){ ### loop columns
		for (rown in c(1 : dim(color_matrix)[1])){ ### loop rows
			### add rectangle
			rect( (coln-1)*colbin_len, (rown-1)*rowbin_len, coln*colbin_len, rown*rowbin_len, col = color_matrix[rown, coln], border=border_color, lwd = 0 )
		}
	}
	dev.off()
}

#################################################### 
############ read input files
#################################################### 
###### read signal matrix 
heatmap_save_type = png

### read signal matrix file
signal_matrix_od = as.matrix(read.table(signal_matrix_file, header=FALSE, sep='\t'))

### extract signal matrix without info
signal_matrix = signal_matrix_od[ , signal_matrix_start_col:dim(signal_matrix_od)[2] ]
### convert to numeric matrix
class(signal_matrix) = 'numeric'

###### read colnames file
colname_file = read.table(signal_input_list, header=F, sep='\t')
### add colnames
colnames(signal_matrix) = colname_file[,1]

signal_matrix_color = matrix("#FFFFFF", nrow=dim(signal_matrix)[1], ncol=dim(signal_matrix)[2])
###### read color file
signal_range_color = read.table(signal_range_color_file, header=F, sep='\t')

### add colnames
for (i in seq(1,dim(signal_range_color)[1])){
	color_rgb = unlist(strsplit(toString(signal_range_color[i,2]), split = ','))
	color_rgb = as.numeric(color_rgb)
	cor_col_i = rgb( color_rgb[1]/255, color_rgb[2]/255, color_rgb[3]/255 )
	signal_matrix_color[signal_matrix==signal_range_color[i,1]] = cor_col_i
}

####################################################
###### plot heatmap
colnames(signal_matrix_color) = colname_file[,1]
print(head(signal_matrix_color))
color_heatmap(signal_matrix_color, output_filename, heatmap_save_type, heatmap_boarder_col)




