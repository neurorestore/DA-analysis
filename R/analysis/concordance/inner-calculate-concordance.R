setwd("~/git/DA-analysis")

library(tidyverse)
library(magrittr)
options(stringsAsFactors = FALSE)
library(argparse)
library(Matrix)
library(stringr)
library(sparseMatrixStats)
library(fgsea)

parser = ArgumentParser(prog = 'inner-calculate-concordance.R')
grid = read.delim("sh/grids/analysis/calculate-concordance.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

source("R/analysis/concordance/concordance_misc_functions.R")

bulk_da = readRDS(args$bulk_filename)
sc_da = readRDS(args$sc_filename)

bulk_data = readRDS(args$bulk_data_filename)
bulk_mat = bulk_data$mat
bulk_meta = bulk_data$meta
sc_data = readRDS(args$bulk_data_filename)
sc_mat = sc_data$mat
sc_meta = sc_data$meta

peak_sizes = readRDS(args$peak_size_filenames)
overlapping_windows = readRDS(args$overlapping_windows_filename)

cell_type_paper_da_type = bulk_da %>% 
	select(cell_type, paper_da_type) %>%
	distinct()

output_df = data.frame()

for (i in 1:nrow(cell_type_paper_da_type)){
	curr_cell_type = cell_type_paper_da_type$cell_type[i]
	curr_paper_da_type = cell_type_paper_da_type$paper_da_type[i]
	print(paste0('curr cell type: ', curr_cell_type, '; curr paper da type:', curr_paper_da_type))
	bulk_da_2 %<>% filter(cell_type == curr_cell_type, paper_da_type == curr_paper_da_type)
	sc_da_2 %<>% filter(cell_type == curr_cell_type, paper_da_type == curr_paper_da_type)

	# add gene statistics
	labels = unlist(str_split(curr_cell_type, '_vs_')[[1]])
	temp_bulk_meta = bulk_meta %>% filter(label %in% c(labels[1], labels[2]))
	temp_bulk_mat = bulk_mat[,temp_bulk_meta$barcode]
	temp_bulk_bin_mat = temp_bulk_mat
	bulk_gene_stats = data.frame(
		'mean'=rowMeans(temp_bulk_mat),
		'percentage_open'=rowMeans(temp_bulk_bin_mat)
	) %>%
	left_join(peak_sizes) %>%
	left_join(overlapping_windows)
	bulk_da_2 %<>% left_join(bulk_gene_stats, by = 'gene') 

	temp_sc_meta = sc_meta %>% filter(label %in% c(labels[1], labels[2]))
	temp_sc_mat = sc_mat[,temp_sc_meta$barcode]
	temp_sc_bin_mat = temp_sc_mat
	sc_gene_stats = data.frame(
		'mean'=rowMeans(temp_sc_mat),
		'percentage_open'=rowMeans(temp_sc_bin_mat)
	) %>%
	left_join(peak_sizes) %>%
	left_join(overlapping_windows)
	sc_da_2 %<>% left_join(sc_gene_stats, by = 'gene') 

	## experiments
	# multiome, multiome go -- non_overlapping windows
	# multiome, multiome go -- lowly expressed genes
	# matched bulk, multiome -- logFC filtering
	experimental_res = experimental_design(bulk_da_2, sc_da_2, args$experiment)

	for (curr_filter in names(experimental_res)) {
		bulk_da_3 = experimental_res[[curr_filter]]$bulk_da
		sc_da_3 = experimental_res[[curr_filter]]$sc_da
		output_df = calculate_concordance(output_df, bulk_da_3, sc_da_3, k=args$k, filter = curr_filter)
	}
}

if (nrow(output_df) > 0) {
	saveRDS(output_df, args$output_filename)
}