setwd("~/git/DA-analysis")

library(tidyverse)
library(magrittr)
options(stringsAsFactors = FALSE)
library(argparse)
library(tidyverse)
library(magrittr)
library(Matrix)
library(DropletUtils)

parser = ArgumentParser(prog = 'inner-get-downsample-matrix.R')
grid = read.delim("sh/grids/preprocessing/get-downsample-matrix.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

source("R/functions/detect_system.R")
input = readRDS(args$input_filename)

downsample_proportions = c(100, 200, 500, 1000, 2000, 5000, 10000)
final_output = list()

comparison_grid = readRDS('metadata/comparisons.rds') %>%
	filter(paper == args$paper)

mean_cell_count = readRDS('metadata/matched_bulk_comparisons_counts_per_cell.rds')
for (downsample_proportion in downsample_proportions) {
	downsample_out_list = list()

	for (i in 1:nrow(comparison_grid)){
		compar1 = unname(comparison_grid$compar1[i])
		compar2 = unname(comparison_grid$compar2[i])

		temp_mat = input$mat
		temp_meta = input$metadata 
		temp_meta %<>%
			filter(
				barcode %in% colnames(temp_mat),
				label %in% c(compar1, compar2)
			)
		temp_mat = temp_mat[,temp_meta$barcode]

		temp_mat = temp_mat[,temp_meta$barcode]
		temp_mat = temp_mat[rowSums(temp_mat) > 0,]

		downsample_proportion_factor = mean_cell_count %>%
			filter(compar1 == !!compar1, compar2 == !!compar2) %>%
			pull(mean_counts_per_cell)
		downsample_proportion_factor = downsample_proportion/downsample_proportion_factor
		temp_mat = downsampleMatrix(temp_mat, downsample_proportion_factor)
		temp_mat = temp_mat[rowSums(temp_mat) > 0,]
		downsample_out_list[[paste0(compar1, '_vs_', compar2)]] = list('mat'=temp_mat,
			'metadata'=temp_meta)
	}
	final_output[[paste0('downsample_proportion=', downsample_proportion)]] = downsample_out_list
}

saveRDS(final_output, args$output_filename)