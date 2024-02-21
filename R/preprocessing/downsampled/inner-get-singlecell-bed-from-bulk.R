setwd("~/git/DA-analysis")

library(argparse)
library(stringr)
library(tidyverse)
library(magrittr)
library(Seurat)
library(hdf5r)
library(reshape2)
library(Matrix)
library(stringi)
library(GenomicRanges)
library(plyr)

parser = ArgumentParser(prog = 'inner-get-singlecell-bed-from-bulk.R')
grid = read.delim("sh/grids/preprocessing/get-singlecell-bed-from-bulk.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

set.seed(args$seed)

bulk_files = list.files(args$bulk_snap_dir, full.names=T)
bulk_files = bulk_files[endsWith(bulk_files, '.snap')]
replicates_cells = split(1:args$n_cells, cut(seq_along(1:args$n_cells), length(bulk_files), labels = FALSE))

combined_fragments_df = data.frame()

for (i in 1:length(bulk_files)){
	snap_input_file = bulk_files[i]
	sample_name = gsub('-sorted.snap', '', basename(snap_input_file))
	print(paste0('Loading ', sample_name))
	h5_file = H5File$new(snap_input_file, mode = "r")
	barcodes = readDataSet(h5_file[['BD']][['name']])

	fragments_df = data.frame(
		'chrom' = readDataSet(h5_file[['FM']][['fragChrom']]),
		'start' = readDataSet(h5_file[['FM']][['fragStart']]),
		'len' = readDataSet(h5_file[['FM']][['fragLen']])
		) %>% 
		mutate(end = start + len - 1) %>%
		relocate(end, .after=start)
	
	number_of_cells = replicates_cells[[i]]
	number_of_fragments = floor(length(number_of_cells) * args$fragments_proportion * nrow(fragments_df))
	fragments_idx = sample(1:nrow(fragments_df), number_of_fragments, replace=T)

	new_fragments_df = fragments_df[fragments_idx,]
	new_fragments_df$barcode = sample(paste0(sample_name, '_Cell_', number_of_cells), number_of_fragments, replace=T)
	new_fragments_df %<>%
		select(-len)

	rownames(new_fragments_df) = NULL

	combined_fragments_df = rbind(combined_fragments_df, new_fragments_df)
}

colnames(combined_fragments_df) = c('chrom', 'chromStart', 'chromEnd', 'name')
write.table(combined_fragments_df, gzfile(args$output_bed_filename), row.names=F, sep='\t')