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
signal_matrix_od = as.matrix(read.table(signal_matrix_file, header=FALSE))
print(head(signal_matrix_od))

### extract signal matrix without info
signal_matrix = signal_matrix_od[ , signal_matrix_start_col:dim(signal_matrix_od)[2] ]
### convert to numeric matrix
class(signal_matrix) = 'numeric'

###### read colnames file
colname_file = read.table(signal_input_list, header=F)
### add colnames
colnames(signal_matrix) = colname_file[,2]

signal_matrix_color = signal_matrix
###### read color file
signal_range_color = as.matrix(read.table(signal_range_color_file, header=F))
print(head(signal_range_color))

### replace NA by background signal in background color table
print(head(signal_matrix))
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

####################################################
###### plot heatmap
color_heatmap(signal_matrix_color, output_filename, heatmap_save_type, heatmap_boarder_col)





