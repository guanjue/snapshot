library(data.table)


####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
data_mat_file = args[1]
peak_signal_state_list = args[2]


### read mat
d = as.data.frame(fread(data_mat_file))

### colnames 
colname_file = read.table(peak_signal_state_list, header=F, sep='\t')
ct_list = apply(colname_file, 1, function(x) as.character(x[1]) )

### add colnames
colnames(d) = c('#chr','start','end','cCREsID', ct_list)

### write output
write.table(d, data_mat_file, quote=F, sep='\t', col.names=T, row.names=F)
