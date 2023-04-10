import os
import numpy as np
from subprocess import call
from collections import Counter
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.decomposition import FastICA

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return data0

################################################################################################
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()


def QTnorm(tar, ref):
	return(np.sort(ref)[(tar.argsort()).argsort()])

################################################################################################
### get convert bedtools window output to matrix of pk and intersect function label info
def function_label_info(input_bedtools_window, id_col, lb_col, pk_col, function_col):
	#data_info_matrix = function_label_info(mark_bed_file+'.tmp01.txt', cCRE_pk_id_col, IDEAS_state_label_col+cCRE_bed_colnum, cCRE_start_col, IDEAS_state_start_col)
	data_info0=open(input_bedtools_window, 'r')
	### read DNA region orders
	data_info=[]
	for records in data_info0:
		tmp=[x.strip() for x in records.split('\t')]
		#print(tmp)
		### get intersect region; midpoint dist; TF peak length
		if ((int(tmp[function_col-1]) - int(tmp[pk_col-1]))>=0) and ((int(tmp[function_col]) - int(tmp[pk_col]))<=0) :
			### function Bin >= pk region
			tmp_vec = [ tmp[id_col-1], tmp[lb_col-1], int(tmp[function_col])-int(tmp[function_col-1]), (float(tmp[function_col])+float(tmp[function_col-1])-float(tmp[pk_col])-float(tmp[pk_col-1]))/2, int(tmp[function_col])-int(tmp[function_col-1]) ]
		elif ((int(tmp[function_col-1]) - int(tmp[pk_col-1]))<0) and ((int(tmp[function_col]) - int(tmp[pk_col]))>0) :
			### function Bin < pk region
			tmp_vec = [ tmp[id_col-1], tmp[lb_col-1], int(tmp[pk_col])-int(tmp[pk_col-1]), (float(tmp[function_col])+float(tmp[function_col-1])-float(tmp[pk_col])-float(tmp[pk_col-1]))/2, int(tmp[function_col])-int(tmp[function_col-1]) ]
		elif ((int(tmp[function_col-1]) - int(tmp[pk_col-1]))<0) and ((int(tmp[function_col]) - int(tmp[pk_col]))<=0) :
			### function Bin upstream < pk region upstream & function Bin downstream <= pk region downstream
			tmp_vec = [ tmp[id_col-1], tmp[lb_col-1], int(tmp[function_col])-int(tmp[pk_col-1]), (float(tmp[function_col])+float(tmp[function_col-1])-float(tmp[pk_col])-float(tmp[pk_col-1]))/2, int(tmp[function_col])-int(tmp[function_col-1]) ]
		elif ((int(tmp[function_col-1]) - int(tmp[pk_col-1]))>=0) and ((int(tmp[function_col]) - int(tmp[pk_col]))>0) :
			### function Bin upstream >= pk region upstream & function Bin downstream > pk region downstream
			tmp_vec = [ tmp[id_col-1], tmp[lb_col-1], int(tmp[pk_col])-int(tmp[function_col-1]), (float(tmp[function_col])+float(tmp[function_col-1])-float(tmp[pk_col])-float(tmp[pk_col-1]))/2, int(tmp[function_col])-int(tmp[function_col-1]) ]
		data_info.append(tmp_vec)
	data_info0.close()
	###### return output dict
	return (data_info)

################################################################################################
### merge state bed for each cell-type
def merge_ct_state(input_bed):
	#input_bed = 'MEP1.bbd.sort.OD.bed'
	bedmat = read2d_array(input_bed, str)
	unique_state = np.unique(bedmat[:,3])
	call('rm -f ' + input_bed+'.MS.bed', shell=True)
	###
	for s in unique_state:
		#print(s)
		call('cat ' + input_bed + ' | awk -F \'\t\' -v OFS=\'\t\' -v state_i='+str(s)+' \'{if ($4==state_i) print $0}\' > ' + input_bed+'.'+str(s)+'.bed', shell=True)
		call('bedtools merge -i ' + input_bed+'.'+str(s)+'.bed' + ' > ' + input_bed+'.'+str(s)+'.M.bed', shell=True)
		call('cat ' + input_bed+'.'+str(s)+'.M.bed' + ' | awk -F \'\t\' -v OFS=\'\t\' -v state_i='+str(s)+' \'{print $0,state_i}\' >> ' + input_bed+'.MS.bed', shell=True)
		call('rm ' + input_bed+'.'+str(s)+'.M.bed', shell=True)
		call('rm ' + input_bed+'.'+str(s)+'.bed', shell=True)

################################################################################################
### get peak's function labels
def get_cRE_function_state(data_info_matrix, id_col, lb_col, cover_col, middist_col, functionlen, bed_od_file, bed_od_idcol, outputname):
	### read DNA region function state info matrix
	pk_id_list = []
	data_function1={}
	data_function1_maxcover={} ### coverage size
	data_function1_middist={} ### midpoint dist
	data_function1_statelen={} ### function state len
	### initialize problem counter
	k=0
	for info in data_info_matrix:
		pk_id = info[id_col-1]
		lb_tmp = info[lb_col-1]
		### creat pk_id_list for keeping the id order in output
		pk_id_list.append(pk_id)
		if not (pk_id in data_function1):
			data_function1[pk_id] = lb_tmp
			data_function1_maxcover[pk_id] = info[cover_col-1]
			data_function1_middist[pk_id] = info[middist_col-1]
			data_function1_statelen[pk_id] = info[functionlen-1]
		elif (lb_tmp!='0') and (data_function1[pk_id]=='0'):
			### if interesect non-0-state, use non-0-state replace 0 state
			data_function1[pk_id] = lb_tmp
			data_function1_maxcover[pk_id] = info[cover_col-1]
			data_function1_middist[pk_id] = info[middist_col-1]
			data_function1_statelen[pk_id] = info[functionlen-1]
		elif info[cover_col-1] > data_function1_maxcover[pk_id] and lb_tmp!='0':
			### if multiple cover; select the highest covering state
			data_function1[pk_id] = lb_tmp
			data_function1_maxcover[pk_id] = info[cover_col-1]
			data_function1_middist[pk_id] = info[middist_col-1]
			data_function1_statelen[pk_id] = info[functionlen-1]
		elif info[cover_col-1] == data_function1_maxcover[pk_id]: ### if 2 states cover the same region with same length
			if np.absolute(info[middist_col-1]) < np.absolute(data_function1_middist[pk_id]) and lb_tmp!='0': 
				### if cover the same; check mid point distance
				data_function1[pk_id] = lb_tmp
				data_function1_maxcover[pk_id] = info[cover_col-1]
				data_function1_middist[pk_id] = info[middist_col-1]
				data_function1_statelen[pk_id] = info[functionlen-1]
			elif info[middist_col-1] == data_function1_middist[pk_id]: ### if 2 states cover the same region with same length; with same midpoint dist
				if info[functionlen-1] < data_function1_statelen[pk_id] and lb_tmp!='0':
					### if cover same & mid point distance same; check state len 
					data_function1[pk_id] = lb_tmp
					data_function1_maxcover[pk_id] = info[cover_col-1]
					data_function1_middist[pk_id] = info[middist_col-1]
					data_function1_statelen[pk_id] = info[functionlen-1]
				else: ### if 2 states cover the same region with same length; with same midpoint dist; with same state length ...give attention!
					k=k+1
					print('problem!')
					print(k)
	### read original bed file to get the pk id list
	bed_od_file=open(bed_od_file,'r')
	bed_od_id_list = []
	for records in bed_od_file:
		bed_od_id_list.append(records.split()[bed_od_idcol-1])
	bed_od_file.close()
	### write function label output
	result=open(outputname,'w')
	for pkid in bed_od_id_list:
		if pkid in data_function1:
			tmp=data_function1[pkid]
			result.write(tmp+'\n')
		else:
			tmp=records
			result.write('NA'+'\n')
	result.close()


################################################################################################
### get index/signal matrix
def get_mark_matrix(peak_bed, peak_info_column, mark_list, output_file, method, sort_sigbed, script_folder, signal_type='peak', genome_size_file=None):
	#######
	### sort input bed files
	sort_bed_file = peak_bed + '.sort.bed'
	call('cp ' + sort_bed_file + ' ' + output_file, shell=True)
	### For state bigbed file, get 200bp bin
	if method == 'window':
		print('make 200bp windows')
		call('bedtools makewindows -g '+genome_size_file+' -w 200 > genome.200bp.bed', shell=True)
		call('sort -k1,1 -k2,2n genome.200bp.bed > genome.200bp.sort.bed', shell=True)
		call('rm genome.200bp.bed', shell=True)
	##############################
	### generate index mark matrix
	mark_list_vec = open(mark_list, 'r')
	for mark_bed in mark_list_vec:	
		tmp = [x.strip() for x in mark_bed.split('\t')]
		#print(tmp)
		#######
		### read bianry label file list
		if signal_type=='signal':
			mark_bed_bw_file = tmp[2]
		elif signal_type=='peak':
			mark_bed_bw_file = tmp[1]
			### sort bianry label bed files
			if sort_sigbed == 'T':
				call('sort -k1,1 -k2,2n ' + mark_bed_bw_file + ' > ' + mark_bed_bw_file+'.sort.bed', shell=True)
			else:
				call('cp ' + mark_bed_bw_file + ' ' + mark_bed_bw_file+'.sort.bed', shell=True)
		elif signal_type=='state':
			mark_bed_bw_file = tmp[3]
			### sort bianry label bed files
			call('bigBedToBed ' + mark_bed_bw_file + ' ' + mark_bed_bw_file+'.sort.OD.bed', shell=True)
			call('bedtools map -a genome.200bp.sort.bed -b '+ mark_bed_bw_file+'.sort.OD.bed' + ' -null 0 -c 4 -o max > ' + mark_bed_bw_file + '.sort.bed', shell=True)
			call('rm ' + mark_bed_bw_file+'.sort.OD.bed', shell=True)
		#######
		### use bedtools to generate the index/signal matrix
		if method == 'intersect':
			### used bedtools intersect to get the binary label of each peak
			call('bedtools intersect -c -a ' + sort_bed_file + ' -b ' + mark_bed_bw_file+'.sort.bed' + ' > ' + mark_bed_bw_file+'.tmp01.txt', shell=True)
			call('cat ' + mark_bed_bw_file+'.tmp01.txt' + ' | awk -F \'\t\' -v OFS=\'\t\' \'{if ($5<=1) print $0; else print $1, $2, $3, $4, 1}\' > ' + mark_bed_bw_file+'.tmp01b.txt', shell=True)
			call('mv ' + mark_bed_bw_file+'.tmp01b.txt' + ' ' + mark_bed_bw_file+'.tmp01.txt', shell=True)
			call('rm '+mark_bed_bw_file+'.sort.bed', shell=True)
		elif method == 'map':
			### used bedtools map to get the average signal of each peak
			#call('bedtools map -c ' + str(signal_col) + ' -null 0 -o mean -a ' + sort_bed_file + ' -b ' + mark_bed_file+'.sort.bed' + ' > ' + mark_bed_bw_file+'.tmp01.txt', shell=True)
			call('bigWigAverageOverBed ' + mark_bed_bw_file + ' ' + sort_bed_file + ' ' + mark_bed_bw_file+'.tab', shell=True)
			call('sort -k1,1n'+ ' ' + mark_bed_bw_file+'.tab' + ' > ' + mark_bed_bw_file+'.sort.tab', shell=True)
			call('paste ' + sort_bed_file + ' ' + mark_bed_bw_file+'.sort.tab' + ' | awk -F \'\t\' -v OFS=\'\t\' \'{print $1, $2, $3, $4, $9}\' > ' + mark_bed_bw_file+'.tmp01.txt', shell=True)
			call('rm'+ ' ' +mark_bed_bw_file+'.tab'+ ' ' +mark_bed_bw_file+'.sort.tab', shell=True)
		elif method == 'window':
			### used bedtools map to get the average signal of each peak
			merge_ct_state(mark_bed_bw_file+'.sort.bed')
			call('bedtools window -a ' + sort_bed_file + ' -b ' + mark_bed_bw_file+'.sort.bed.MS.bed' + ' -w 0 > ' + mark_bed_bw_file+'.tmp01.txt', shell=True)
			### convert bedtools window output to matrix of pk and intersect function label info (intersect region; midpoint dist; TF peak length)
			data_info_matrix = function_label_info(mark_bed_bw_file+'.tmp01.txt', 4, 8, 2, 6)
			### get peak's function labels based on intersect region; midpoint dist; TF peak length
			get_cRE_function_state(data_info_matrix, 1, 2, 3, 4, 5, sort_bed_file, 4, mark_bed_bw_file+'.tmp01.txt')
			call('rm '+mark_bed_bw_file+'.sort.bed', shell=True)
			call('rm '+mark_bed_bw_file+'.sort.bed.MS.bed', shell=True)
		### cut the map number column
		call('cut -f'+ str(peak_info_column) +" -d$'\t' " + mark_bed_bw_file+'.tmp01.txt' + ' | awk -F \'\t\' -v OFS=\'\t\' \'{if ($1=="NA") print 0; else print $1}\' > ' + mark_bed_bw_file+'.tmp02.txt', shell=True)
		### cbind to matrix
		call('paste ' + output_file + ' ' + mark_bed_bw_file+'.tmp02.txt' + ' > ' + output_file+'.tmp.txt' + ' && mv ' + output_file+'.tmp.txt ' + output_file, shell=True)
		### remove tmp files
		call('rm ' + mark_bed_bw_file+'.tmp01.txt' + ' ' + mark_bed_bw_file+'.tmp02.txt', shell=True)
	if method == 'window':
		call('rm genome.200bp.sort.bed', shell=True)
	mark_list_vec.close()

################################################################################################
### get merged peaks
def merge_pk(peak_signal_state_list, outputname, script_folder):
	import os.path
	import os
	### 
	if os.path.isfile('all_pk.bed'):
		call('rm all_pk.bed', shell=True)
	### read filenames in peak list
	for file_info in open(peak_signal_state_list, 'r'):
		filename = file_info.split('\t')[1]
		call('cat ' + filename + ' >> all_pk.bed', shell=True)
	### sort merge_pk
	call('sort -k1,1 -k2,2n all_pk.bed > all_pk.sort.bed', shell=True)
	### merge peak
	outputfile_name = outputname + '.sort.bed'
	call('bedtools merge -i all_pk.sort.bed > ' + outputfile_name, shell=True)
	### add pk id
	call('cat ' + outputfile_name + ' | awk -F \'\t\' -v OFS=\'\t\' \'{print $1, $2, $3, NR}\' > ' + outputfile_name + '.tmp.txt', shell=True)
	call('mv ' + outputfile_name + '.tmp.txt ' + outputfile_name, shell=True)
	### rm tmp files
	call('rm all_pk.bed all_pk.sort.bed', shell=True)
	print('cCRE number:')
	call('wc -l ' + outputfile_name, shell=True)
	### return filename
	return(outputfile_name)

################################################################################################
### vector most frequent element
def frequent(element_vector, state_color_mat):
	#print(Counter(element_vector))
	most_freq_element = Counter(element_vector).most_common(1)[0][0]
	###### return output dict
	return most_freq_element

def rgb_to_hex(rgb):
	rbg_code = '#%02x%02x%02x' % (rgb[0], rgb[1], rgb[2])
	return rbg_code

def frequent(element_vector, state_color_mat):
	most_freq_element = Counter(element_vector)
	rgb_merge = np.array([0,0,0])
	denominator = 0
	#print(most_freq_element)
	### sort count
	most_freq_element_sort = dict(sorted(most_freq_element.items(), key=lambda item: item[1], reverse=True))
	### sum count
	all_freq_sum = sum(most_freq_element_sort.values())
	proportion_n = 0
	k = 0
	for i in most_freq_element:
		k = k+1
		proportion_n = proportion_n + (most_freq_element[i])
		rgb_merge = rgb_merge + state_color_mat[int(i),:] * (most_freq_element[i])
		denominator = denominator + (most_freq_element[i])
		if (proportion_n / all_freq_sum) >= 0.6:
			break
		#if k >=1:
		#	break
	#print(rgb_merge)
	rgb_merge = np.round(rgb_merge / denominator, 0)
	#print(rgb_merge)
	rgb_merge = rgb_merge.astype(int)
	#print(rgb_merge)
	rgb_merge = rgb_to_hex(np.round(rgb_merge, 0))
	#print(rgb_merge)
	###### return output dict
	return rgb_merge


################################################################################################
### column based calculation
def matrix_col_cal(matrix, function, para=None):
	### get total column number
	column = matrix.shape[1]
	### for loop columns
	col_score_list = []
	for i in range(0, column):
		### extract column vector
		col_vec = matrix[:,i]
		### calculation
		if para is None:
			col_vec_score = function(col_vec)
		elif para.shape[1] == column:
			col_vec_score = function(col_vec, para[:,i])
		else:
			col_vec_score = function(col_vec, para)
		col_score_list.append(col_vec_score)
	col_score_list = np.array(col_score_list)
	###### return output dict
	return col_score_list

################################################################################################
### column based calculation
def matrix_col_cal_merge(matrix, function, state_color_mat, para=None):
	### get total column number
	column = matrix.shape[1]
	### for loop columns
	col_score_list = []
	for i in range(0, column):
		### extract column vector
		col_vec = matrix[:,i]
		### calculation
		if para is None:
			col_vec_score = function(col_vec, state_color_mat)
		elif para.shape[1] == column:
			col_vec_score = function(col_vec, para[:,i])
		else:
			col_vec_score = function(col_vec, para)
		col_score_list.append(col_vec_score)
	col_score_list = np.array(col_score_list)
	###### return output dict
	return col_score_list

################################################################################################
### use QDA to rescue peaks with rare pattern
def QDA_rescue(index_label_vector, signal_matrix, index_X, count_threshold, qda_num):
	##################
	### use QDA to reassign labels
	index_label_vector = np.array(index_label_vector)
	index_label_vector_od = index_label_vector
	change_num_array = []
	print('X_num:')
	print(np.sum(index_label_vector==index_X))
	qda_num = 1
	if qda_num > 0:
		for i in range(0,qda_num):
			#print('QDA iteration: ' + str(i))
			filter_rows = index_label_vector!=index_X
			signal_matrix_train = signal_matrix[filter_rows,:] + np.random.uniform(low=0.000001, high=0.00001, size=signal_matrix[filter_rows,:].shape)
			ica = FastICA(n_components=10)
			#signal_matrix_train_ica = ica.fit_transform(signal_matrix_train)
			#signal_matrix_ica = ica.fit_transform(signal_matrix)
			index_label_vector_train = index_label_vector[filter_rows]
			clf = QuadraticDiscriminantAnalysis(priors=np.repeat(1/np.unique(index_label_vector_train).shape[0], np.unique(index_label_vector_train).shape[0] ))
			clf.fit(signal_matrix_train, index_label_vector_train)

			### rescued index_vector
			index_label_vector_pre = index_label_vector
			index_label_vector = clf.predict(signal_matrix)
			
			### print posterior probabilities of classification
			index_label_vector_post = clf.predict_proba(signal_matrix)
			index_label_vector_post_max = np.max(index_label_vector_post, axis=1)
			print(index_label_vector_post.shape)
			print(index_label_vector_post_max.shape)

			print('X_num: after rescue')
			print(np.sum(index_label_vector==index_X))

			print(np.sum(index_label_vector!=index_label_vector_od))
			index_label_vector[index_label_vector_post_max<0.5] = index_label_vector_od[index_label_vector_post_max<0.5]
			
			print('X_num: after rescue filter non condident ones')
			print(np.sum(index_label_vector==index_X))
			print(np.sum(index_label_vector!=index_label_vector_od))
			print('index_label_vector')
			print('index_label_vector')
			print('index_label_vector')
			print(index_label_vector)
			
			#index_label_vector[filter_rows] = index_label_vector_pre[filter_rows]

			### print the number of peak label changes
			#print('QDA changed label number: ')
			change_num = np.sum(index_label_vector_pre!=index_label_vector)
			#print(change_num)
			change_num_array.append(change_num)
			if change_num == 0:
				break
	else:
		for i in range(0,qda_num):
			change_num_array.append(0)

	### get change_num_array
	change_num_array = np.array(change_num_array)
	change_num_array = change_num_array.reshape(change_num_array.shape[0], 1)

	### generate rescued signal dict 1st
	index_set_mean_signal_matrix_dict_QDA_rescue = {}
	index_label_vector_QDA_rescue = []
	index_uniq_vec = []

	### get new index set matrix
	for index, index_signal in zip(index_label_vector, signal_matrix):
		if not (index in index_set_mean_signal_matrix_dict_QDA_rescue):
			index_set_mean_signal_matrix_dict_QDA_rescue[ index ] = ''
			index_label_vector_QDA_rescue.append(index)
			index_uniq_vec.append(index)
		else:
			index_label_vector_QDA_rescue.append(index)

	index_label_vector_QDA_rescue = np.array(index_label_vector_QDA_rescue)
	#print(index_uniq_vec)
	#print('QDA changed label number: ')
	change_num = np.sum(index_label_vector_QDA_rescue!=index_label_vector)
	#print(change_num)
	
	### filter by count_thresh
	to_indexX = {}
	for index in index_uniq_vec:
		#print(index)
		#print('OD count: '+str(np.sum(index_label_vector_od == index)))
		index_new_num = np.sum(index_label_vector_QDA_rescue == index)
		#print('QDA rescued count: '+str(index_new_num))
		if index_new_num < count_threshold:
			to_indexX[index] = ''

	### generate rescued signal dict 2nd
	index_set_mean_signal_matrix_dict_QDA_rescue = {}
	index_label_vector_QDA_rescue = []
	index_uniq_vec = []

	### get new index set matrix
	for index, index_signal in zip(index_label_vector, signal_matrix):
		if index in to_indexX:
			index = index_X
		if not (index in index_set_mean_signal_matrix_dict_QDA_rescue):
			index_set_mean_signal_matrix_dict_QDA_rescue[ index ] = [ index_signal ]
			index_label_vector_QDA_rescue.append(index)
			index_uniq_vec.append(index)
		else:
			index_set_mean_signal_matrix_dict_QDA_rescue[ index ].append(index_signal)
			index_label_vector_QDA_rescue.append(index)

	### return index_label_vector_QDA_rescue & index_set_mean_signal_matrix_dict_QDA_rescue
	return { 'index_label_vector_QDA_rescue': index_label_vector_QDA_rescue, 'index_set_mean_signal_matrix_dict_QDA_rescue':index_set_mean_signal_matrix_dict_QDA_rescue, 'change_num_array': change_num_array }

################################################################################################
### get pass count threshold index dict
def pass_count_thresh(index_matrix, count_threshold):
	##################
	###### extract index_set pass count threshold
	index_vector = []
	for index_array in index_matrix:
		index_label = ''
		for i in range(0,len(index_array)-1):
			index_label = index_label + index_array[i] + '_'
		index_label = index_label + index_array[len(index_array)-1]
		### append to index vector
		index_vector.append(index_label)
	### index_vector 2 np array
	index_vector = np.array(index_vector)
	#print(index_vector)
	### index peak counts (dict)
	index_uniq_count_dict = Counter(index_vector)
	#print(index_uniq_count_dict)
	### get index dict of index pass the count threshold
	pass_thresh_index_dict = {}
	for index in index_uniq_count_dict:
		if index_uniq_count_dict[ index ] >= count_threshold:
			pass_thresh_index_dict[ index ] = index_uniq_count_dict[ index ]
	### return pass_thresh_index_dict
	print(pass_thresh_index_dict)
	return { 'pass_thresh_index_dict': pass_thresh_index_dict, 'index_vector':index_vector }

################################################################################################
### get index_set signal matrix
def get_index_set_mean_signal_matrix(signal_matrix_file, pass_thresh_index_dict, count_threshold, index_vector, qda_num, log2_sig='F', normalization_method='F', quantile_1='F', smallnum=0.0):
	##################
	###### get index_set signal matrix
	### read signal matrix
	signal_matrix_od = read2d_array(signal_matrix_file, 'str')
	### bed info
	bed_info = signal_matrix_od[:, 3].reshape(signal_matrix_od.shape[0], 1)
	### signal matrix
	signal_matrix = signal_matrix_od[:, (5-1):].astype(float)
	### QTnorm
	signal_matrix_rowMeans = np.mean(signal_matrix, axis=1)
	#signal_matrix_QT = signal_matrix
	#for i in range(0,signal_matrix.shape[1]):
	#	signal_matrix_QT[:,i] = QTnorm(signal_matrix[:,i], signal_matrix_rowMeans)
	#signal_matrix = signal_matrix_QT


	### cap outliers
	for i in range(0,signal_matrix.shape[1]):
		signal_matrix[signal_matrix[:,i]>np.quantile(signal_matrix[:,i], 0.99999),i] = np.quantile(signal_matrix[:,i], 0.99999)
	### adjust signal
	if log2_sig == 'T':
		print('get_index_set_mean_signal_matrix: signal_matrix')
		print(signal_matrix)
		print('log_sig:')
		signal_matrix = np.log2(signal_matrix+smallnum)
		print(signal_matrix)
	else:
		signal_matrix = signal_matrix+smallnum

	##########################################################################################
	### normalization section
	### normalize by scale, scale_quantile, QTnorm
	### S3norm can only be done after Index-Set clustering (for Identify common-peak-IS and common-background-IS)
	if normalization_method == 'scale_quantile':
		quantile_1 = float(quantile_1)
		signal_matrix_quantile_sig = np.quantile(signal_matrix, quantile_1, axis=0)
		signal_matrix = signal_matrix / signal_matrix_quantile_sig
	elif normalization_method == 'QTnorm':
		### get ref average signal
		signal_matrix_rowMeans = np.mean(signal_matrix, axis=1)
		signal_matrix_QTnorm = signal_matrix
		### QTnorm
		for i in range(0,signal_matrix.shape[1]):
			signal_matrix_QTnorm[:,i] = QTnorm(signal_matrix[:,i], signal_matrix_rowMeans)
		###
		signal_matrix = signal_matrix_QTnorm
	else:
		signal_matrix_mean = np.mean(signal_matrix, axis=0)
		signal_matrix_std = np.std(signal_matrix, axis=0)
		signal_matrix = (signal_matrix -  signal_matrix_mean) / signal_matrix_std	
	##########################################################################################


	### get index set mean signal matrix
	index_set_mean_signal_matrix_dict = {}
	index_set_vector = []
	index_label_vector = []

	### get index_X
	index_vec = index_vector[0].split('_')
	index_X = ''
	index_common_1 = ''
	for i in range(0, len(index_vec)-1):
		index_X = index_X + 'X_'
		index_common_1 = index_common_1 + '1_'
	index_X = index_X + 'X'
	index_common_1 = index_common_1 + '1'

	print('Initial IS number:')
	print(len(np.unique(index_vector)))
	for index, index_signal in zip(index_vector, signal_matrix):
		### if the index_set is not in pass_thresh_index_dict, replace index by X_
		if not (index in pass_thresh_index_dict):
			index = index_X

		### append to index_label_vector for function matrix analysis
		index_label_vector.append(index)
		### get index set mean signal
		if not (index in index_set_mean_signal_matrix_dict):
			index_set_mean_signal_matrix_dict[ index ] = [ index_signal ]
			index_set_vector.append(index)
		else:
			index_set_mean_signal_matrix_dict[ index ].append(index_signal)

	######
	### QDA rescue
	index_set_mean_signal_matrix_dict_QDA_rescue_info = QDA_rescue(index_label_vector, signal_matrix, index_X, count_threshold, 1)
	index_label_vector_QDA_rescue = index_set_mean_signal_matrix_dict_QDA_rescue_info['index_label_vector_QDA_rescue']
	index_set_mean_signal_matrix_dict_QDA_rescue = index_set_mean_signal_matrix_dict_QDA_rescue_info['index_set_mean_signal_matrix_dict_QDA_rescue']
	change_num_array = index_set_mean_signal_matrix_dict_QDA_rescue_info['change_num_array']
	######
	### get index mean signal matrix
	index_set_mean_signal_matrix = []
	index_set_vector_qda = []
	print('check colmeans: ')
	for index in index_set_vector:
		if index in index_set_mean_signal_matrix_dict_QDA_rescue:
			### read signal matrix of each index_set
			index_set_signal_matrix_individual = np.array(index_set_mean_signal_matrix_dict_QDA_rescue[ index ]).astype(float)
			### get mean vector
			#index_set_signal_mean_vector_individual = np.mean(index_set_signal_matrix_individual, axis=0)
			### only use top 100 cCREs to get mean
			index_set_signal_matrix_individual_top = np.sort(index_set_signal_matrix_individual, axis=0)[-100:,:]
			#index_set_signal_mean_vector_individual = np.median(index_set_signal_matrix_individual_top, axis=0)
			index_set_signal_mean_vector_individual = np.median(index_set_signal_matrix_individual, axis=0)
			print(np.median(index_set_signal_matrix_individual_top, axis=0))
			print(np.median(index_set_signal_matrix_individual, axis=0))
			index_set_mean_signal_matrix.append(index_set_signal_mean_vector_individual)
			index_set_vector_qda.append(index)
	######
	### cbind index_set_label and signal matrix
	index_label_vector_QDA_rescue = np.array(index_label_vector_QDA_rescue).reshape(len(index_label_vector_QDA_rescue), 1)
	index_set_vector_qda = np.array(index_set_vector_qda).reshape(len(index_set_vector_qda), 1)
	index_set_mean_signal_matrix = np.array(index_set_mean_signal_matrix)
	### index_set_sort_id
	sort_id = np.argsort(index_set_vector_qda, axis=0)
	sort_id = sort_id[:,0]

	### normalize signal to noise ratio
	colmin = index_set_mean_signal_matrix.min(axis=0)
	if np.sum(index_set_vector_qda[:,0]==index_common_1)!=0:
		common_pk_sig = index_set_mean_signal_matrix[index_set_vector_qda[:,0]==index_common_1,:] - colmin
	else:
		common_pk_sig = index_set_mean_signal_matrix.max(axis=0)
	common_pk_sig_ref = np.max(common_pk_sig)
	colmin_ref = np.median(colmin)
	common_pk_sig_ref = np.max(signal_matrix)
	colmin_ref = np.min(signal_matrix)
	### make colmin_ref > 0
	common_pk_sig_ref = common_pk_sig_ref - colmin_ref
	colmin_ref = 1
	#common_pk_sig_ref = 10
	#colmin_ref = 1
	if normalization_method=='S3norm':
		print('S3norm')
		print(common_pk_sig_ref)
		print(colmin_ref)
		index_set_mean_signal_matrix = (index_set_mean_signal_matrix - colmin) / common_pk_sig * common_pk_sig_ref + colmin_ref
		signal_matrix = (signal_matrix - colmin) / common_pk_sig * common_pk_sig_ref + colmin_ref

	signal_matrix[signal_matrix<0] = 0

	
	### cbind for output
	#write2d_array(signal_matrix, 'outputname'+'.signal_matrix.renormalized.txt')
	index_set_mean_signal_matrix = np.concatenate((index_set_vector_qda, index_set_mean_signal_matrix), axis=1)[sort_id,:]
	index_signal_matrix = np.concatenate((bed_info, index_label_vector_QDA_rescue, signal_matrix), axis=1)
	index_signal_matrix = index_signal_matrix[np.argsort(index_signal_matrix[:,1], axis=0),:] ### sort index label functional state matrix

	return { 'index_set_mean_signal_matrix': index_set_mean_signal_matrix, 'sort_id': sort_id, 'index_label_vector': index_label_vector_QDA_rescue, 'index_set_vector': index_set_vector_qda, 'index_signal_matrix': index_signal_matrix, 'change_num_array': change_num_array }

################################################################################################
### index_set function matrix
def get_index_set_function_matrix(function_matrix_file, pass_thresh_index_dict, sort_id, index_label_vector, index_set_vector, state_color_file):
	##################
	###### get index_set signal matrix
	### read functional state matrix
	function_matrix_od = read2d_array(function_matrix_file, 'str')
	### bed info
	bed_info = function_matrix_od[:, 3].reshape(function_matrix_od.shape[0], 1)
	### functional state matrix
	function_matrix = function_matrix_od[:, (5-1):]
	### get index set functional state matrix
	index_set_function_matrix_dict = {}
	for index, function_vector in zip(index_label_vector, function_matrix):
		index = index[0]
		### get index set mean signal
		if not (index in index_set_function_matrix_dict):
			index_set_function_matrix_dict[ index ] = [ function_vector ]
		else:
			index_set_function_matrix_dict[ index ].append(function_vector)
	######
	state_color_mat = read2d_array(state_color_file, str)
	state_color_vec = state_color_mat[:,1]
	state_color_mat = np.array([string.split(",") for string in state_color_vec])
	state_color_mat = state_color_mat.astype(int)
	print('state_color_vec')
	print('state_color_vec')
	print('state_color_vec')
	print(state_color_vec)
	print(state_color_mat)
	print(state_color_mat.shape)

	### get index set most frequent functional state matrix
	index_set_mostfreqfun_matrix = []
	for index in index_set_vector:
		index = index[0]
		matrix = np.array(index_set_function_matrix_dict[ index ], dtype = str)
		### extract the most frequent functional state in each index set
		#matrix_col = matrix_col_cal(matrix, frequent)
		#print(matrix)
		matrix_col = matrix_col_cal_merge(matrix, frequent, state_color_mat)
		#print(matrix_col)
		### append most frequent functional state
		index_set_mostfreqfun_matrix.append( matrix_col )
	### cbind index_set_label and most frequent functional state matrix
	index_set_mostfreqfun_matrix = np.concatenate((index_set_vector, index_set_mostfreqfun_matrix), axis=1)[sort_id,:]
	index_fun_matrix = np.concatenate((bed_info, index_label_vector, function_matrix), axis=1)
	index_fun_matrix = index_fun_matrix[np.argsort(index_fun_matrix[:,1], axis=0),:] ### sort index label functional state matrix
	return { 'index_set_mostfreqfun_matrix': index_set_mostfreqfun_matrix, 'index_fun_matrix': index_fun_matrix }

################################################################################################
### index set sth some score matrix
def get_index_set(outputname, signal_matrix_file, count_threshold, log2_sig, normalization_method, quantile_1, smallnum, qda_num, function_method, state_color_file, script_folder, peak_signal_state_list):
	### read index matrix file
	index_matrix_file = outputname + '.index.matrix.txt'
	index_matrix_od = read2d_array(index_matrix_file, 'str')
	### bed info: chrom start end id
	bed_info = index_matrix_od[:, 0:5]
	### index matrix
	index_matrix = index_matrix_od[:, (5-1):]
	print(index_matrix)
	print(count_threshold)

	##################
	###### extract index_set pass count threshold
	pass_thresh_index_info = pass_count_thresh(index_matrix, count_threshold)
	pass_thresh_index_dict = pass_thresh_index_info['pass_thresh_index_dict']
	index_vector = pass_thresh_index_info['index_vector']

	##################
	###### get index_set signal matrix
	print('get index_set mean signal matrix...')
	signal_matrix_file = outputname + '.signal.matrix.txt'
	print(index_vector)
	print(index_vector.shape)
	print(len(pass_thresh_index_info))
	index_set_mean_signal_info = get_index_set_mean_signal_matrix(signal_matrix_file, pass_thresh_index_dict, count_threshold, index_vector, qda_num, log2_sig, normalization_method, quantile_1, smallnum)
	index_set_mean_signal_matrix = index_set_mean_signal_info['index_set_mean_signal_matrix']
	sort_id = index_set_mean_signal_info['sort_id']
	index_label_vector = index_set_mean_signal_info['index_label_vector']
	index_set_vector = index_set_mean_signal_info['index_set_vector']
	index_signal_matrix = index_set_mean_signal_info['index_signal_matrix']
	change_num_array = index_set_mean_signal_info['change_num_array']

	write2d_array(change_num_array, outputname+'.change_num_array.txt')
	write2d_array(index_set_mean_signal_matrix, outputname+'.meansig.txt')
	write2d_array(index_signal_matrix, outputname+'.sig.txt')
	print('get index_set mean signal matrix...DONE')




	##################
	###### get index_set function matrix
	print('get index_set most freqent functional state matrix...')
	function_matrix_file = outputname + '.function.matrix.txt'
	if function_method == 'mostfreq':
		### if functional information is most frequent state information
		print(index_label_vector)
		print(index_label_vector.shape)
		print(index_set_vector)
		print(index_set_vector.shape)
		#index_set_freqfun_info = get_index_set_function_matrix(function_matrix_file, MetaIS_id_vec_unique[:,0], sort_id, MetaIS_id_vec, MetaIS_id_vec_unique, state_color_file)
		index_set_freqfun_info = get_index_set_function_matrix(function_matrix_file, pass_thresh_index_dict, sort_id, index_label_vector, index_set_vector, state_color_file)
		index_set_mostfreqfun_matrix = index_set_freqfun_info['index_set_mostfreqfun_matrix']
		index_fun_matrix = index_set_freqfun_info['index_fun_matrix']
		write2d_array(index_set_mostfreqfun_matrix, outputname+'.indexset_fun.txt')
		write2d_array(index_fun_matrix, outputname+'.fun.txt')
	elif function_method == 'mean':
		### if functional information is numerical information
		index_set_funcion_mean_signal_info = get_index_set_mean_signal_matrix(function_matrix_file, pass_thresh_index_dict, count_threshold, index_vector)
		index_set_funcion_mean_signal_matrix = index_set_funcion_mean_signal_info['index_set_mean_signal_matrix']
		sort_id = index_set_funcion_mean_signal_info['sort_id']
		index_label_vector = index_set_funcion_mean_signal_info['index_label_vector']
		index_set_vector = index_set_funcion_mean_signal_info['index_set_vector']
		index_signal_matrix = index_set_funcion_mean_signal_info['index_signal_matrix']
		write2d_array(index_set_funcion_mean_signal_matrix, outputname+'.indexset_fun.txt')
		write2d_array(index_signal_matrix, outputname+'.fun.txt')	
	else:
		print('ERROR: get index_set most freqent functional state matrix...METHOD not found!!!')
	print('get index_set most freqent functional state matrix...DONE')

################################################################################################

################################################################################################
def snapshot(master_peak_bed, peak_signal_state_list, genome_size_file, outputname, count_threshold, siglog2, normalization_method, quantile_1, sigsmallnum, function_color_file, cd_tree, input_folder, output_folder, script_folder, qda_num, use_user_thresh, have_function_state_files, index_matrix_txt, signal_matrix_txt, function_state_matrix_txt):
	### set working directory
	os.chdir(input_folder)

	signal_high_color = 'red'
	signal_low_color = 'white'
	index_set_sig_matrix_start_col = 2
	index_sig_matrix_start_col = 3
	index_set_fun_matrix_start_col = 2
	index_fun_matrix_start_col = 3
	index_set_boarder_color = 'gray'
	index_boarder_color = 'NA'
	function_method = 'mostfreq'

	################################################################################################
	print('Check if have function state data...')
	print('have_function_state_files..: ' + have_function_state_files)
	################################################################################################	


	################################################################################################
	###### Used user provided master peaks OR Generate merged master peak -> get sorted peak bed files
	print('01: Get master peak list')
	if master_peak_bed==None:
		merge_pk_name = merge_pk(peak_signal_state_list, outputname, script_folder)
	else:
		call('cat '+master_peak_bed + ' | awk -F \'\t\' -v OFS=\'\t\' \'{print $1, $2, $3, NR}\' > ' + outputname + '.sort.bed', shell=True)
		merge_pk_name = outputname + '.sort.bed'
	print('merged peak generated: ' + merge_pk_name + '... Done')
	################################################################################################


	################################################################################################
	###### get used matrices
	print('02: Get used matrices')
	### get bed binary matrix
	peak_label_column = 5
	output_file_index = outputname + '.index.matrix.txt'
	signal_type = 'peak'
	sort_sigbed = 'T'
	method = 'intersect'
	print('get binary matrix...')
	if index_matrix_txt=='F':
		print('Generate Index matrix')
		get_mark_matrix(outputname, peak_label_column, peak_signal_state_list, output_file_index, method, sort_sigbed, script_folder, signal_type)
	else:
		print('Use user provided index matrix')
		call('tail -n+2 '+index_matrix_txt + ' > ' + output_file_index, shell=True)
		# write ct list file by using the first row of index matrix
		ct_list_file = peak_signal_state_list
		call('head -n1 '+index_matrix_txt + ' | awk -F \'\t\' -v OFS=\'\t\' \'{for(i=5;i<=NF;i++) print $i}\' > ' + ct_list_file, shell=True)

	### get signal matrix
	peak_signal_column = 5
	output_file_signal = outputname + '.signal.matrix.txt'
	signal_type = 'signal'
	sort_sigbed = 'T'
	method = 'map'
	print('get signal matrix...')
	if signal_matrix_txt=='F':
		print('Generate signal matrix')
		get_mark_matrix(outputname, peak_signal_column, peak_signal_state_list, output_file_signal, method, sort_sigbed, script_folder, signal_type)
	else:
		print('Use user provided signal matrix')
		call('tail -n+2 '+signal_matrix_txt + ' > ' + output_file_signal, shell=True)

	if have_function_state_files != 'F':
		### get function label matrix
		output_file_function = outputname + '.function.matrix.txt'
		signal_type = 'state'
		sort_sigbed = 'F'
		print('get function matrix...')
		peak_function_column = 1
		method = 'window'
		if function_state_matrix_txt == 'F':
			print('Generate functional state matrix')
			get_mark_matrix(outputname, peak_function_column, peak_signal_state_list, output_file_function, method, sort_sigbed, script_folder, signal_type, genome_size_file)
		else:
			print('Use user provided functional state matrix')
			call('tail -n+2 '+function_state_matrix_txt + ' > ' + output_file_function, shell=True)
	################################################################################################


	################################################################################################
	###### Determine Index-Set counts threshold
	###### plot index_count density plot
	print('03: IS count threshold')
	print('get NB_count_thresh')
	#call('time Rscript ' + script_folder + 'plot_density.N.sm.R ' + outputname + '.index.matrix.txt' + ' ' + str(count_threshold) + ' ' + str(use_user_thresh), shell=True)
	call('time Rscript ' + script_folder + 'plot_density.R ' + outputname + '.index.matrix.txt' + ' ' + str(count_threshold) + ' ' + str(use_user_thresh), shell=True)
	###### read NB count threshold
	#NB_count_thresh = int(read2d_array(outputname + '.index.matrix.txt.NB_count_thresh.txt', 'str')[0])
	NB_count_thresh = int(read2d_array(outputname + '.index.matrix.txt.NB_count_thresh.txt', 'str')[0])
	#NB_count_thresh = 0
	print('NB model based thresholds:')
	print(NB_count_thresh)
	###### if NB_count_thresh is greater than user provided count threshold, then use NB_count_thresh
	#if NB_count_thresh > 1000:
	#	print('replace user provide count_threshold by NB_count_thresh')
	#	count_threshold = NB_count_thresh
	#	print(count_threshold)
	#else:
	#	print('use user provide count_threshold')
	if count_threshold != 0:
		count_threshold = count_threshold
	else:
		count_threshold = NB_count_thresh

	#count_threshold = 500
	print('count_threshold: '+str(count_threshold))
	################################################################################################


	################################################################################################
	###### get index_set matrices
	print('04: get Index-Sets matrix')
	get_index_set(outputname, output_file_signal, count_threshold, siglog2, normalization_method, quantile_1, sigsmallnum, qda_num, function_method, function_color_file, script_folder, peak_signal_state_list)
	################################################################################################


	################################################################################################
	###### get Meta-ISs matrix
	print('05: get Meta-ISs matrix')
	call('time Rscript ' + script_folder + 'get_meta_IS.R ' + outputname+'.meansig.txt' + ' ' + peak_signal_state_list + ' ' + outputname+'.sig.txt' + ' ' + outputname + ' ' + merge_pk_name + ' ' + outputname+'.function.matrix.txt' + ' ' + have_function_state_files, shell=True)
	### read Meta-IS-IDs
	MetaIS_id_vec = read2d_array(outputname+'.metaISid.mat.txt', str)[:,15]
	print(MetaIS_id_vec.shape)
	MetaIS_id_vec = MetaIS_id_vec.astype(int)
	#MetaIS_id_vec = MetaIS_id_vec.reshape(MetaIS_id_vec.shape[0], 1)
	MetaIS_id_vec_unique = np.unique(MetaIS_id_vec)
	MetaIS_id_vec_unique = MetaIS_id_vec_unique.astype(int)
	#MetaIS_id_vec_unique = MetaIS_id_vec_unique.reshape(MetaIS_id_vec_unique.shape[0], 1)
	print('MetaIS_id_vec...')
	print(MetaIS_id_vec)
	print(MetaIS_id_vec.shape)
	print(np.min(MetaIS_id_vec_unique))
	print(np.max(MetaIS_id_vec_unique))
	################################################################################################


	################################################################################################
	###### get index_set matrices
	if have_function_state_files != 'F':
		print('06: get function matrix')
		index_set_function_matrix = read2d_array(outputname+'.function.matrix.txt', str)[:,4:]
		index_set_function_matrix = index_set_function_matrix.astype('int')
		###
		state_color_mat = read2d_array(function_color_file, str)
		state_color_vec = state_color_mat[:,1]
		state_color_mat = np.array([string.split(",") for string in state_color_vec])
		state_color_mat = state_color_mat.astype(int)
		### get index set most frequent functional state matrix
		META_index_set_mostfreqfun_matrix = []
		#for index in MetaIS_id_vec_unique:
		for i in range(np.min(MetaIS_id_vec_unique), np.max(MetaIS_id_vec_unique)+1):
			index = i
			print(index)
			matrix = np.array(index_set_function_matrix[MetaIS_id_vec==index,:], dtype = float)
			### extract the most frequent functional state in each index set
			matrix_col = matrix_col_cal_merge(matrix, frequent, state_color_mat)
			### append most frequent functional state
			META_index_set_mostfreqfun_matrix.append( matrix_col )
		### cbind index_set_label and most frequent functional state matrix
		META_index_set_mostfreqfun_matrix = np.concatenate((MetaIS_id_vec_unique.reshape(MetaIS_id_vec_unique.shape[0], 1), META_index_set_mostfreqfun_matrix), axis=1)
		print(META_index_set_mostfreqfun_matrix)
		write2d_array(META_index_set_mostfreqfun_matrix, outputname+'.metaIS.fun.txt')
		call('sort -k1,1n ' + outputname+'.metaIS.fun.txt' + ' > ' + outputname+'.metaIS.fun.sort.txt', shell=True)
		call('time Rscript ' + script_folder + 'plot_rect.1.R ' + outputname+'.metaIS.fun.sort.txt' + ' ' + outputname+'.Meta_indexset_fun.png' + ' ' + peak_signal_state_list + ' ' + index_set_boarder_color, shell=True)
	################################################################################################


	################################################################################################
	############ plot figures
	###### for signal
	### plot binary heatmaps 
	print('07: plot heatmaps')
	print('use plot_rect to plot signal index & index set heatmap...')
	call('cat ' + outputname+'.meansig.txt' + ' | awk -F \'\t\' -v OFS=\'\t\' \'{print $1}\' | awk -F \'_\' -v OFS=\'\t\' \'{print $0}\'  > ' + outputname+'.index_binary_mat.txt.tmp1', shell=True)
	call('cat ' + outputname+'.meansig.txt' + ' | awk -F \'\t\' -v OFS=\'\t\' \'{print $1}\' | awk -F \'_\' -v OFS=\'\t\' \'{for(i=1;i<=NF;i++) printf "%s\\t",$i; printf "\\n"}\'  > ' + outputname+'.index_binary_mat.txt.tmp2', shell=True)
	call('paste '+ outputname+'.index_binary_mat.txt.tmp1' + ' ' + outputname+'.index_binary_mat.txt.tmp2' + ' > ' + outputname+'.index_binary_mat.txt', shell=True)
	call('rm '+ outputname+'.index_binary_mat.txt.tmp1' + ' ' + outputname+'.index_binary_mat.txt.tmp2', shell=True)
	call('time Rscript ' + script_folder + 'plot_rect_sig.R ' + outputname+'.index_binary_mat.txt' + ' ' + outputname+'.index_binary_mat.png' + ' ' + peak_signal_state_list + ' ' + str(index_set_sig_matrix_start_col) + ' ' + 'black' + ' ' + 'white' + ' ' + index_set_boarder_color + ' ' + siglog2 + ' ' + str(sigsmallnum), shell=True)
	### plot heatmaps 
	print('08: plot heatmaps')
	print('use plot_rect to plot signal index & index set heatmap...')
	call('time Rscript ' + script_folder + 'plot_rect_sig.R ' + outputname+'.meansig.txt' + ' ' + outputname+'.meansig.png' + ' ' + peak_signal_state_list + ' ' + str(index_set_sig_matrix_start_col) + ' ' + signal_high_color + ' ' + signal_low_color + ' ' + index_set_boarder_color + ' ' + siglog2 + ' ' + str(sigsmallnum), shell=True)
	### plot tree
	print('09: plot tree')
	print('plot mean signal of cell differentiation tree')
	call('if [ ! -d signal_tree_IndexSet ]; then mkdir signal_tree_IndexSet; fi', shell=True)
	call('time Rscript ' + script_folder + 'plot_tree.R ' + outputname+'.meansig.txt' + ' ' + cd_tree + ' ' + peak_signal_state_list + ' ' + str(index_set_sig_matrix_start_col) + ' ' + signal_high_color + ' ' + signal_low_color + ' ' + siglog2 + ' ' + str(sigsmallnum), shell=True)
	call('mv *tree.pdf signal_tree_IndexSet/', shell=True)
	### plot violin
	print('10: plot violin')
	print('plot signal violin plot & create index set bed files ...')
	call('if [ ! -d signal_violin_IndexSet ]; then mkdir signal_violin_IndexSet; fi', shell=True)
	call('time Rscript ' + script_folder + 'plot_sig_violin.R ' + outputname+'.sig.txt' + ' ' + peak_signal_state_list + ' ' + 'violin.pdf' + ' ' + outputname + '.sort.bed' , shell=True)
	call('mv *violin.pdf signal_violin_IndexSet/', shell=True)
	call('if [ ! -d '+outputname+'_IndexSets_bed_files ]; then mkdir '+outputname+'_IndexSets_bed_files; fi', shell=True)
	call('mv *.index_set.bed '+outputname+'_IndexSets_bed_files/', shell=True)

	###### for functional state
	if have_function_state_files != 'F':
		### plot heatmaps functional
		if function_method == 'mostfreq':
			### plot functional state
			print('11: use plot_rect to plot function heatmap...')
			#call('time Rscript ' + script_folder + 'plot_rect.R ' + outputname+'.indexset_fun.txt' + ' ' + outputname+'.indexset_fun.png' + ' ' + peak_signal_state_list + ' ' + function_color_file + ' ' + str(index_set_fun_matrix_start_col) + ' ' + index_set_boarder_color, shell=True)
			call('time Rscript ' + script_folder + 'plot_rect.1.R ' + outputname+'.indexset_fun.txt' + ' ' + outputname+'.indexset_fun.png' + ' ' + peak_signal_state_list + ' ' + index_set_boarder_color, shell=True)
			### plot tree
			print('12: plot functional state of cell differentiation tree')
			call('if [ ! -d fun_tree_IndexSet ]; then mkdir fun_tree_IndexSet; fi', shell=True)
			call('time Rscript ' + script_folder + 'plot_tree_multi_color.1.R ' + outputname+'.indexset_fun.txt' + ' ' + cd_tree + ' ' + function_color_file + ' ' + peak_signal_state_list + ' ' + str(index_set_fun_matrix_start_col), shell=True)
			call('mv *tree.pdf fun_tree_IndexSet/', shell=True)
			### plot bar
			print('13: plot functional state barplot...')
			call('if [ ! -d fun_bar_IndexSet ]; then mkdir fun_bar_IndexSet; fi', shell=True)
			call('time Rscript ' + script_folder + 'plot_fun_bar.R ' + outputname+'.fun.txt' + ' ' + peak_signal_state_list + ' ' + function_color_file + ' ' + 'bar.pdf' , shell=True)
			call('mv *bar.pdf fun_bar_IndexSet/', shell=True)


	### plot tree
	print('14: plot mean signal of cell differentiation tree')
	call('if [ ! -d signal_tree_MetaIS ]; then mkdir signal_tree_MetaIS; fi', shell=True)
	call('time Rscript ' + script_folder + 'plot_tree.R ' + outputname+'.metaISid.meansigmat.txt' + ' ' + cd_tree + ' ' + peak_signal_state_list + ' ' + str(index_set_sig_matrix_start_col) + ' ' + signal_high_color + ' ' + signal_low_color + ' ' + siglog2 + ' ' + str(sigsmallnum), shell=True)
	call('mv *tree.pdf signal_tree_MetaIS/', shell=True)
	### plot violin
	print('15: plot signal violin plot & create index set bed files ...')
	call('if [ ! -d signal_violin_MetaIS ]; then mkdir signal_violin_MetaIS; fi', shell=True)
	call('time Rscript ' + script_folder + 'plot_sig_violin.metaIS.R ' + outputname+'.MetaIS.forviolin.sig.txt' + ' ' + peak_signal_state_list + ' ' + 'violin.pdf' + ' ' + outputname + '.sort.bed' , shell=True)
	call('mv *violin.pdf signal_violin_MetaIS/', shell=True)

	if have_function_state_files != 'F':
		### plot tree
		print('16: plot functional state of cell differentiation tree')
		call('if [ ! -d fun_tree_MetaIS ]; then mkdir fun_tree_MetaIS; fi', shell=True)
		call('time Rscript ' + script_folder + 'plot_tree_multi_color.1.MetaIS.R ' + outputname+'.metaIS.fun.sort.txt' + ' ' + cd_tree + ' ' + function_color_file + ' ' + peak_signal_state_list + ' ' + str(index_set_fun_matrix_start_col), shell=True)
		call('mv *tree.pdf fun_tree_MetaIS/', shell=True)
		###
		print('17: plot functional state barplot...')
		call('if [ ! -d fun_bar_MetaIS ]; then mkdir fun_bar_MetaIS; fi', shell=True)
		call('time Rscript ' + script_folder + 'plot_fun_bar.MetaIS.R ' + outputname+'.metaIS_all.fun.txt' + ' ' + peak_signal_state_list + ' ' + function_color_file + ' ' + '.metaIS.bar.pdf' , shell=True)
		call('mv *bar.pdf fun_bar_MetaIS/', shell=True)
		### 
	print('18: add colnames...')
	call('time Rscript ' + script_folder + 'add_colnames.R ' + outputname+'.index.matrix.txt' + ' ' + peak_signal_state_list, shell=True)
	call('time Rscript ' + script_folder + 'add_colnames.R ' + outputname+'.signal.matrix.txt' + ' ' + peak_signal_state_list, shell=True)
		
	if have_function_state_files != 'F':
		call('time Rscript ' + script_folder + 'add_colnames.R ' + outputname+'.function.matrix.txt' + ' ' + peak_signal_state_list, shell=True)



	################################################################################################
	###### mv all output to output folder
	call('if [ -d ' + output_folder + ' ]; then rm -r ' + output_folder + '; mkdir ' + output_folder + '; fi', shell=True)
	call('if [ ! -d ' + output_folder + ' ]; then mkdir ' + output_folder + '; fi', shell=True)
	call('mkdir ' + output_folder + '/IndexSet_output', shell=True)
	call('mkdir ' + output_folder + '/MetaIndexSet_output', shell=True)

	### mv merge peak file
	call('mv ' + outputname + '.sort.bed' + ' ' + output_folder, shell=True)
	### mv matrix
	call('mv ' + outputname + '.index.matrix.txt' + ' ' + output_folder, shell=True)
	call('mv ' + outputname + '.signal.matrix.txt' + ' ' + output_folder, shell=True)
	if have_function_state_files != 'F':
		call('mv ' + outputname + '.function.matrix.txt' + ' ' + output_folder, shell=True)
	###
	call('mv ' + outputname + '.index_binary_mat.txt' + ' ' + output_folder + '/IndexSet_output', shell=True)
	call('mv ' + outputname + '.meansig.txt' + ' ' + output_folder + '/IndexSet_output', shell=True)
	call('mv ' + outputname + '.sig.txt' + ' ' + output_folder + '/IndexSet_output', shell=True)
	call('mv ' + outputname + '.change_num_array.txt' + ' ' + output_folder + '/IndexSet_output', shell=True)
	call('mv ' + outputname + '.index.matrix.txt.NB_count_thresh.txt' + ' ' + output_folder + '/IndexSet_output', shell=True)
	call('mv ' + outputname + '.index.matrix.txt.IS_num_vs_FDRthresh.pdf' + ' ' + output_folder + '/IndexSet_output', shell=True)
	if have_function_state_files != 'F':
		call('mv ' + outputname + '.indexset_fun.txt' + ' ' + output_folder + '/IndexSet_output', shell=True)
		call('mv ' + outputname + '.fun.txt' + ' ' + output_folder + '/IndexSet_output', shell=True)

	### Meta-Index-Set matrix
	call('mv ' + outputname + '.IS_metaISid.mat.final.txt' + ' ' + output_folder, shell=True)
	call('mv ' + outputname + '.meta_cluster_cCRE_ave_merge.pdf' + ' ' + output_folder, shell=True)
	call('mv ' + outputname + '.IS_cluster.AIC.pdf' + ' ' + output_folder + '/MetaIndexSet_output', shell=True)
	call('mv ' + outputname + '.meta_cluster_ave.pdf' + ' ' + output_folder + '/MetaIndexSet_output', shell=True)
	call('mv ' + outputname + '.meta_cluster_label.pdf' + ' ' + output_folder + '/MetaIndexSet_output', shell=True)
	call('mv ' + outputname + '.metaISid.meansigmat.txt' + ' ' + output_folder + '/MetaIndexSet_output', shell=True)
	if have_function_state_files != 'F':
		call('mv ' + outputname + '.metaIS.fun.txt' + ' ' + output_folder + '/MetaIndexSet_output', shell=True)

	### index set merged figures
	call('mv *.png ' + output_folder, shell=True)
	### individual index set figures & bed files
	### Index-Set results
	call('mv '+outputname+'_IndexSets_bed_files ' + output_folder, shell=True)
	call('mv signal_tree_IndexSet ' + output_folder + '/IndexSet_output', shell=True)
	call('mv signal_violin_IndexSet ' + output_folder + '/IndexSet_output', shell=True)
	if have_function_state_files != 'F':
		call('mv fun_tree_IndexSet ' + output_folder + '/IndexSet_output', shell=True)
		call('mv fun_bar_IndexSet ' + output_folder + '/IndexSet_output', shell=True)
	### Meta-Index-Set results
	call('mv '+outputname+'_MetaISs_bed_files ' + output_folder, shell=True)
	call('mv signal_tree_MetaIS ' + output_folder + '/MetaIndexSet_output', shell=True)
	call('mv signal_violin_MetaIS ' + output_folder + '/MetaIndexSet_output', shell=True)
	if have_function_state_files != 'F':
		call('mv fun_tree_MetaIS ' + output_folder + '/MetaIndexSet_output', shell=True)
		call('mv fun_bar_MetaIS ' + output_folder + '/MetaIndexSet_output', shell=True)

	### remove tmp files
	call('rm ' + outputname + '.metaISid.mat.txt', shell=True)
	call('rm ' + outputname + '.MetaIS.forviolin.sig.txt', shell=True)
	if have_function_state_files != 'F':
		call('rm ' + outputname + '.metaIS.fun.sort.txt', shell=True)
		call('rm ' + outputname + '.metaIS_all.fun.txt', shell=True)
	################################################################################################


##############################################
# time python /Users/gzx103/Documents/zhang_lab/projects/scripts/snapshot/bin/snapshot.py -p peak_list.txt -n atac_4cell -t 5 -s signal_list.txt -l F -z F -x 0.01 -f function_list.txt -m mostfreq -c function_color_list.txt -e cd_tree.txt -i /Users/gzx103/Documents/zhang_lab/projects/scripts/snapshot/test_data/input_data/ -o /Users/gzx103/Documents/zhang_lab/projects/scripts/snapshot/test_data/output_result/ -b /Users/gzx103/Documents/zhang_lab/projects/scripts/snapshot/bin/
##############################################

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hm:p:n:t:l:z:d:x:f:c:e:i:o:s:q:u:b:a:r:g:j:")
	except getopt.GetoptError:
		print('time python3 snapshot.py error, Please check parameter settings')
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print('time python3 snapshot.py, Please check Snapshot github page (https://github.com/guanjue/snapshot) for details')
			sys.exit()
		elif opt=="-m":
			master_peak_bed=str(arg.strip())
		elif opt=="-p":
			peak_signal_state_list=str(arg.strip())
		elif opt=="-n":
			outputname=str(arg.strip())
		elif opt=="-t":
			count_threshold=int(arg.strip())
		elif opt=="-l":
			siglog2=str(arg.strip())
		elif opt=="-z":
			normalization_method=str(arg.strip())
		elif opt=="-d":
			quantile_1=str(arg.strip())
		elif opt=="-x":
			sigsmallnum=float(arg.strip())
		elif opt=="-f":
			genome_size_file=str(arg.strip())
		elif opt=="-c":
			function_color_file=str(arg.strip())
		elif opt=="-e":
			cd_tree=str(arg.strip())
		elif opt=="-i":
			input_folder=str(arg.strip())
		elif opt=="-o":
			output_folder=str(arg.strip())
		elif opt=="-s":
			script_folder=str(arg.strip())
		elif opt=="-q":
			qda_num=int(arg.strip())
		elif opt=='-u':
			use_user_thresh = str(arg.strip())
		elif opt=='-b':
			index = str(arg.strip())
		elif opt=='-a':
			index_matrix_txt = str(arg.strip())
		elif opt=='-r':
			signal_matrix_txt = str(arg.strip())
		elif opt=='-g':
			function_state_matrix_txt = str(arg.strip())
		elif opt=='-j':
			have_function_state_files = str(arg.strip())


	############ Default parameters
	###### required parameters
	try:
		print('User provide peak-bed signal-bigWig state-bigBed file list: -p '+str(peak_signal_state_list))
		print('User provide output filename: -n '+str(outputname))
		print('User provide genome chrom size file: -f '+str(genome_size_file))
		print('User provide Epigenetic State color list file: -c'+str(function_color_file))
		print('User provide pairwise cell-type tree file: -e '+str(cd_tree))
		print('User provide input folder: -i '+str(input_folder))
		print('User provide output folder: -o '+str(output_folder))
		print('User provide source code folder: -s '+str(script_folder))
	except NameError:
		print('Missing required parameter(s): time python3 snapshot.py -p peak_signal_state_list[required] -n outputname[required] -f genome_chrom_size_file[required] -c function_color_list[required] -e cell_development_tree[required] -i input_folder[required] -o output_folder[required] -b script_folder[required]')	
		return()

	###
	###### optional parameters
	try:
		print('User provide master_peak_bed: -m '+str(master_peak_bed))
	except NameError:
		print('Default: Generating Pool & Merged master_peak_bed file')
		master_peak_bed = None

	try:
		print('User provide minimum number of bins per Index-set: -t '+str(count_threshold))
	except NameError:
		print('Default: The minimum number of bins per Index-set is set to 0, which means let Snapshot determine this number automatically')
		count_threshold = 0

	try:
		print('add small number to avoid 0: -x '+str(sigsmallnum))
	except NameError:
		print('Default: add small number 0.01')
		sigsmallnum = 0.01

	try:
		print('User decide if use log scale for signals: -l '+str(siglog2))
	except NameError:
		print('Default: use linear scale for signals')
		siglog2 = 'F'

	try:
		print('User decide if convert signal to Z score: -z '+str(normalization_method))
	except NameError:
		print('Default: use raw signals')
		normalization_method = 'F'

	try:
		print('User decide if quantile / median for scale: -d '+str(quantile_1))
	except NameError:
		print('Default: use default quantile / median for scale if normalization_method=scale_quantile')
		quantile_1 = 0.5

	try:
		print('User provide QDA round number: -q '+str(qda_num))
	except NameError:
		print('Default: 1 round of QDA rescue')
		qda_num = 1

	try:
		print('Use user given thresh for index-set number lower lim: -u '+str(use_user_thresh))
	except NameError:
		print('Default: Not use user given thresh for index-set number lower lim')
		use_user_thresh = 'F'

	try:
		print('Use user provide info for have_function_state_files: -b '+str(have_function_state_files))
	except NameError:
		print('Default: Not use user provide info for have_function_state_files')
		have_function_state_files = 'T'

	try:
		print('Use user provide index_matrix.txt: -a '+str(index_matrix_txt))
	except NameError:
		print('Default: Not use user provide index_matrix.txt')
		index_matrix_txt = 'F'

	try:
		print('Use user provide signal_matrix.txt: -a '+str(signal_matrix_txt))
	except NameError:
		print('Default: Not use user provide signal_matrix.txt')
		signal_matrix_txt = 'F'

	try:
		print('Use user provide function_state_matrix.txt: -a '+str(function_state_matrix_txt))
	except NameError:
		print('Default: Not use user provide function_state_matrix.txt')
		function_state_matrix_txt = 'F'




	print('Starting Snapshot pipeline .......')
	snapshot(master_peak_bed, peak_signal_state_list, genome_size_file, outputname, count_threshold, siglog2, normalization_method, quantile_1, sigsmallnum, function_color_file, cd_tree, input_folder, output_folder, script_folder, qda_num, use_user_thresh, have_function_state_files, index_matrix_txt, signal_matrix_txt, function_state_matrix_txt)

if __name__=="__main__":
	main(sys.argv[1:])



