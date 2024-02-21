import sys
import os
import re
import pandas as pd
from itertools import compress
from collections import Counter
import argparse
import gzip

# set working directory
git_dir = os.path.expanduser("~/git/DA-analysis")
python_dir = git_dir + "/python"
os.chdir(python_dir)
sys.path.append(python_dir)

### dynamically build CLI
parser = argparse.ArgumentParser()
## build the CLI
grid_file = git_dir + '/sh/grids/preprocessing/inner-dex-Pliner-2018-scATAC-chunks.txt'
grid = pd.read_csv(grid_file, sep='\t')
for arg_name in list(grid):
	param_name = '--' + arg_name
	param_dtype = str(grid[arg_name].dtype)
	# convert to pandas
	param_type = {'object': str,
				  'int64': int,
				  'float64': float,
				  'bool': bool
				  }[param_dtype]
	parser.add_argument(param_name, type=param_type)

# parse all arguments
args = parser.parse_args()
print(args)

chunk_dir = args.chunk_dir
chunk_file = args.chunk_file
chunk_file_index_size = 1250000

srr_code, read, chunk_idx = chunk_file.split('_')

base_dir = '/scratch/st-bkkwon-1/DA-analysis'

metadata_dir = '{}/metadata/'.format(base_dir)
raw_data_dir = '{}/raw_data/matched_bulk/Pliner_2018_scATAC/'.format(base_dir)
processed_data_dir = '{}/processed_data/matched_bulk/Pliner_2018_scATAC/chunks/'.format(base_dir)

if srr_code == 'SRR6652584':
	scATAC_meta_file = '{}/GSM2970930_sciATAC_HSMM1_indextable.txt.gz'.format(metadata_dir)
elif srr_code == 'SRR6652585':
	scATAC_meta_file = '{}/GSM2970931_sciATAC_HSMM2_indextable.txt.gz'.format(metadata_dir)

scATAC_meta = pd.read_csv(
	scATAC_meta_file, compression='gzip', index_col=0,
	dtype = {'barcode': str, 'timepoint': str, 'barcode_index': 'uint16'}
	)
scATAC_meta = scATAC_meta[(scATAC_meta['timepoint'] == '0 hours') | (scATAC_meta['timepoint'] == '72 hours')]

out_file = '{}/{}.dex.fastq'.format(processed_data_dir, chunk_file)
write_line = False
metadata_idxs = pd.DataFrame()
barcode_idx = 0
new_file = True
filter_index_df = True

scATAC_index_file = '{}Pliner_scATAC_{}_index_conversion_filtered.txt.gz'.format(metadata_dir, srr_code)
scATAC_index = pd.read_csv(
	scATAC_index_file, compression='gzip',
	index_col=0,
	dtype = {'ReadNumber': 'uint32', 'barcode_index': 'uint16'}
	)
scATAC_index = scATAC_index.sort_values(by=['ReadNumber']).reset_index(drop=True)

fastq_file = '{}/{}'.format(chunk_dir, chunk_file)
output_text = ''

with open(fastq_file, 'r') as fastq_in:
	line_counter, line_idx = 0, 0
	while True:
		line = fastq_in.readline()
		if line == '':
			break
		line_idx += 1
		line_counter += 1
		print(line_idx, line_counter)
		# if line_idx > 10000:
		# 	break
		line = line.strip()
		if len(line) > 0:
			if (line[0]=='@' and line_counter == 1):
				# print(line)
				temp_line = line.split(' ')[1]
				temp_idx, temp_read = temp_line.split('/')
				assert(temp_read == read)
				temp_idx = int(temp_idx)
				if filter_index_df:
					filter_index_df = False
					start_idx = temp_idx
					end_idx = start_idx + chunk_file_index_size
					scATAC_index = scATAC_index[(scATAC_index['ReadNumber'] >= start_idx) & (scATAC_index['ReadNumber'] <= end_idx)]
				if (any(scATAC_index['ReadNumber']==temp_idx)):
					meta_line = scATAC_index[scATAC_index['ReadNumber']==temp_idx]
					barcode_index = meta_line['barcode_index'].item()
					barcode = scATAC_meta[scATAC_meta['barcode_index']==barcode_index]['barcode'].item()
					write_line = True
					line = "@{}:".format(barcode)+line[1:]
			if write_line == True:
				output_text += line + '\n'
			if line_counter == 4:
				line_counter = 0
				write_line = False

fastq_out = open(out_file, "w")
fastq_out.write(output_text)
fastq_out.close()
fastq_in.close()