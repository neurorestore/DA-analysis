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
library(parallel)
library(qsmooth)
library(SummarizedExperiment)

parser = ArgumentParser(prog = 'inner-run-da-by-comparison.R')
grid = read.delim("sh/grids/analysis/run-da-by-comparison.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

source('R/analysis/differential_accessibility/da_misc_functions.R')

temp = readRDS(args$input_filename)
mat = temp$mat
metadata = temp$metadata
rownames(metadata) = metadata$barcode
mat = mat[,colnames(mat) %in% rownames(metadata)]
metadata = metadata[colnames(mat),]
rm(temp)

metadata = preprocess_meta(metadata, args$input_filename)

paper_da_type = args$paper_da_type
if (!args$alternative_implementations) {
	da_types = da_family_args[[args$da_family]][[args$da_method]]
} 

output_df = data.frame()
compar1 = args$compar1
compar2 = args$compar2

temp_meta = metadata %>%
	filter(
		(label==compar1|label==compar2),
		!is.na(replicate))
	
if (paper_da_type=='cell type'){
	temp_meta %<>%
		mutate(cell_type=paste0(compar1, '_vs_', compar2))
}

rownames(temp_meta)=temp_meta$barcode
temp_mat = mat[,rownames(temp_meta)]
temp_mat = temp_mat[rowSums(temp_mat) > 0,]


# filter experiment
if (args$min_cell_filter != 0 & args$da_family == 'pseudobulk'){
	temp_mat = temp_mat[rowSums(temp_mat >= 1) >= args$min_cell_filter,]
}

if (args$min_cell_percent_filter != 0 & args$da_family == 'pseudobulk'){
	temp_mat = temp_mat[rowSums(temp_mat >= 1) >= floor(ncol(temp_mat)*(as.numeric(args$min_cell_percent_filter)/100)),]	
}

# binarization experiment
if (args$binarization == 't'){
	temp_mat@x[temp_mat@x > 0] = 1
}

# normalization experiment
if (args$normalization == 'qsmooth'){
	peak_matrix_check = grepl('peak_matrix', args$input_filename)
	if (!paper_check | !peak_matrix_check){
		stop('Unable to get GC content for this job')
		quit()
	}
	gc_content_df = readRDS('peaks/GC_content/peaks_gc_content.rds')[[args$paper]]
	rownames(gc_content_df) = gc_content_df$gene
}

# additional experiments
if (!is.na(args$experiment)) {
	
	# matched bulk -- ArchR background experiment
	if (args$experiment == 'ArchR_background') {
		rownames(temp_meta) = temp_meta$barcode
		ArchR_bg_res = readRDS(paste0('ArchR/bg_comparisons/', args$paper, '-ArchR_bg_comparisons.rds'))
		ArchR_bg_list = ArchR_bg_res[[paste0(compar1, '_vs_', compar2)]]
		filtered_barcodes = c(ArchR_bg_list[[1]]$barcode, ArchR_bg_list[[2]]$barcode)
		temp_meta %<>% filter(barcode %in% filtered_barcodes)
	}

	# matched bulk -- spurious peaks
	if (args$experiment == 'spurious_peaks') {
		other_data_dir = 'final_data/matched_bulk/'
		other_temp = readRDS(paste0(other_data_dir, args$paper, '_spurious_peak_matrix.rds'))
		other_mat = other_temp$mat[,colnames(mat)]
		temp_mat = rbind(temp_mat, other_mat)
		rm(other_temp)
	}
}


cell_types = data.frame(table(temp_meta$cell_type)) %>%
	filter(Freq > 3) %>%
	pull(Var1) %>%
	as.character()

output_df = data.frame()
for (curr_cell_type in cell_types) {
	print(curr_cell_type)
	temp_meta_2 = temp_meta %>%
		filter(cell_type == curr_cell_type)
	temp_mat_2 = temp_mat[,temp_meta_2$barcode]
	temp_mat_2 = temp_mat_2[rowSums(temp_mat_2) > 0,]
	if (!args$alternative_implementations) {
		output_df = run_da_wrapper(output_df, temp_mat_2, temp_meta_2, args, gc_content_df = gc_content_df, da_types = da_types)	
	} else {
		output_df = run_alt_DA(mat1, mat2, meta1, meta2, test = args$da_method)
	}
}
	
if (nrow(output_df) > 0){
	output_df %<>% mutate(paper_da_type=paper_da_type)
	saveRDS(output_df, args$output_filename)	
}