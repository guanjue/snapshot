####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
signal_matrix_file = args[1]
output_filename = args[2]
signal_input_list = args[3]
signal_matrix_start_col = as.numeric(args[4])

signal_high_color = args[5]
signal_low_color = args[6]
heatmap_boarder_col = args[7]

log2 = args[8]
smallnum = as.numeric(args[9])


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

######
### read signal matrix file
signal_matrix_od = as.matrix(read.table(signal_matrix_file, header=FALSE))
### extract signal matrix without info
signal_matrix = signal_matrix_od[ , signal_matrix_start_col:dim(signal_matrix_od)[2] ]
if (signal_matrix[dim(signal_matrix)[1],1]=='X'){
	signal_matrix[dim(signal_matrix)[1],] = rep(1,dim(signal_matrix)[2])
}
### convert to numeric matrix
class(signal_matrix) = 'numeric'

###### read colnames file
colname_file = read.table(signal_input_list, header=F)
### add colnames
colnames(signal_matrix) = colname_file[,1]

### log2 transform
if (log2=='T'){
	signal_matrix = log2(signal_matrix+smallnum)
}

### filter range
filter_range = max(signal_matrix) - min(signal_matrix)
### subtract min(filter matrix) (Entropy smaller -> less white filter)
filter_percent = (signal_matrix - min(signal_matrix) ) / filter_range


####################################################
###### convert signal matrix to color matrix (With filter)
###### get r channel score
signal_high_color_rgb = colorRamp(signal_high_color)(1)
signal_low_color_rgb = colorRamp(signal_low_color)(1)

### signal to color
rh = signal_high_color_rgb[1]/255
rl = signal_low_color_rgb[1]/255
### add filter: od_color * filter + bg_color * (1 - filter)
r = rh * filter_percent + rl * (1-filter_percent)
### set upper & lower limit for r
r[r>1] = 1
r[r<0] = 0
###### get g channel score
### signal to color
gh = signal_high_color_rgb[2]/255
gl = signal_low_color_rgb[2]/255
### add filter: od_color * (1-filter) + bg_color * filter
g = gh * filter_percent + gl * (1-filter_percent)
### set upper & lower limit for g
g[g>1] = 1
g[g<0] = 0
###### get b channel score
### signal to color
bh = signal_high_color_rgb[3]/255
bl = signal_low_color_rgb[3]/255
### add filter: od_color * (1-filter) + bg_color * filter
b = bh * filter_percent + bl * (1-filter_percent)
### set upper & lower limit for b
b[b>1] = 1
b[b<0] = 0
###### creat rgb color matrix
signal_matrix_color = NULL
for (i in seq(1,dim(r)[2])){
	#print(i)
	### convert each r,g,b vector to a rgb matrix column
	cor_col_tmp = rgb( r[,i], g[,i], b[,i] )
	### cbind column to a rgb matrix
	signal_matrix_color = cbind(signal_matrix_color, cor_col_tmp)
}
###### replace color matrix cell values
colnames(signal_matrix_color) = colname_file[,1]

####################################################
###### plot heatmap
color_heatmap(signal_matrix_color, output_filename, heatmap_save_type, heatmap_boarder_col)





