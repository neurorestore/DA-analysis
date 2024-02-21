setwd("~/git/DA-analysis")

library(argparse)
library(reshape2)
library(tidyverse)
library(magrittr)
library(Libra)
library(Matrix)
library(Seurat)
library(edgeR)
library(limma)
library(DESeq2)
library(digest)
library(parallel)
library(SummarizedExperiment)
library(qsmooth)

parser = ArgumentParser(prog = 'inner-run-da-by-comparison-pseudo-compar.R')
grid = read.delim("sh/grids/analysis/run-da-by-comparison-pseudo-compar.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

source('R/analysis/differential_accessibility/da_misc_functions.R')
## make sure seed is constant
hash_seed = substr(digest(args$output_filename, algo='xxhash32'), 1, 3)
set.seed(strtoi(hash_seed, 16))


temp = readRDS(args$input_filename)
mat = temp$mat
metadata = temp$metadata
rownames(metadata) = metadata$barcode
mat = mat[,colnames(mat) %in% rownames(metadata)]
metadata = metadata[colnames(mat),]
rm(temp)

if (args$experiment=='bulk_ATAC') {
	metadata$replicate = metadata$barcode
}

metadata = preprocess_meta(metadata, args$input_filename)

paper_da_type = args$paper_da_type
da_types = da_family_args[[args$da_family]][[args$da_method]]

output_df = data.frame()

label_of_interest = args$label
temp_meta = metadata %>%
	filter(
		(label==label_of_interest),
		!is.na(replicate))
replicates = unique(temp_meta$replicate)

if (length(replicates) >= 4){
	new_labels_1 = sample(replicates, floor(length(replicates)/2))	
	temp_meta %<>%
		mutate(
			label_ori = label,
			label = ifelse(
				replicate %in% new_labels_1,
				'pseudogroup1',
				'pseudogroup2'
			)
		)

	compar1 = 'pseudogroup1'
	compar2 = 'pseudogroup2'

	rownames(temp_meta)=temp_meta$barcode
	temp_mat = mat[,rownames(temp_meta)]
	temp_mat = temp_mat[rowSums(temp_mat) > 0,]

	if (args$binarization == 't'){
		temp_mat@x[temp_mat@x > 0] = 1
	}

	gc_content_df = readRDS('peaks/GC_content/Luecken_peaks_gc_content.rds')
	rownames(gc_content_df) = gc_content_df$gene

	for (curr_cell_type in cell_types) {
		print(curr_cell_type)
		temp_meta_2 = temp_meta %>%
			filter(cell_type == curr_cell_type)
		temp_mat_2 = temp_mat[,temp_meta_2$barcode]
		temp_mat_2 = temp_mat_2[rowSums(temp_mat_2) > 0,]
		output_df = run_da_wrapper(temp_mat_2, temp_meta_2, args, gc_content_df = gc_content_df, da_types = da_types)
	}

	if (nrow(output_df) > 0){
		output_df %<>% mutate(paper_da_type=paper_da_type)
		saveRDS(output_df, args$output_filename)
	}
}
