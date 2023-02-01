####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
signal_matrix_file = args[1]
output_filename = args[2]
peak_signal_state_list = args[3]
heatmap_boarder_col = args[4]

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
			#print(color_matrix[rown, coln])
			#print(as.character(color_matrix[rown, coln]))
			rect( (coln-1)*colbin_len, (rown-1)*rowbin_len, coln*colbin_len, rown*rowbin_len, col = as.character(color_matrix[rown, coln]), border=border_color, lwd = 0 )
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
signal_matrix_od = (read.table(signal_matrix_file, header=FALSE, sep='\t', comment.char='~'))

### colnames 
colname_file = read.table(peak_signal_state_list, header=F, sep='\t')
ct_list = apply(colname_file, 1, function(x) as.character(x[1]) )

###
signal_matrix_color = signal_matrix_od[,-1]
colnames(signal_matrix_color) = ct_list
#print(head(signal_matrix_color))
####################################################
###### plot heatmap
color_heatmap(signal_matrix_color, output_filename, heatmap_save_type, heatmap_boarder_col)




