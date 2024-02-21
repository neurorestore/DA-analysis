setwd("~/git/DA-analysis")

library(argparse)
library(stringr)
library(tidyverse)
library(magrittr)
library(Seurat)
library(hdf5r)
library(reshape2)
library(Matrix)

parser = ArgumentParser(prog = 'inner-get-peaks-from-bed.R')
grid = read.delim("sh/grids/preprocessing/get-peaks-from-bed.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

source('R/functions/detect_system.R')
# args$qval = 0.01
macs2_options = paste0(
	'--nomodel --shift 37 --ext 73 --qval ', format(args$qval, scientific=TRUE), ' -B --SPMR --call-summits'
)
message(macs2_options)

setwd(args$data_dir)

if (!file.exists(args$combined_bed_file)){
	files = list.files(args$data_dir)
	bed_files = files[endsWith(files, '.bed')]
	if (length(bed_files) == 0){
		bed_files = files[endsWith(files, '.bed.gz')]
	}

	flag = system2(command="cat", 
			args=c(paste(bed_files, collapse = ' '),
					">", args$combined_bed_file
				)
	)
}

flag = system2(command=args$macs2_path, 
	args=c("callpeak", 
			"-t", args$combined_bed_file, 
			"-f", "BED",
			"-g", args$gsize,
			macs2_options,
			"-n", args$prefix
			)
	)
if (flag != 0) {
	stop("'MACS' call failed")
}