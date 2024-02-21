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

# do downsampling here
downsample_proportions = c(10, 30, 50, 100, 300, 500, 1000)
final_output = list()

comparison_grid = readRDS('metadata/comparisons.rds') %>%
	filter(paper == args$paper)

cell_count_df = readRDS('metadata/cell_count_df.rds') %>%
	mutate(
		cell_type = as.character(cell_type),
		label = as.character(label),
		replicate = as.character(replicate)
	)

cell_types_of_interest = cell_count_df %>%
	filter(
		paper %in% grid$paper,
		cell_type %in% c(comparison_grid$compar1, comparison_grid$compar2)
	) %>%
	group_by(paper, cell_type) %>%
	summarise(cell_count = sum(cell_count)) %>%
	ungroup() %>%
	filter(cell_count > 1000) %>%
	pull(cell_type)	

comparison_grid %<>%
	filter(
		(compar1 %in% cell_types_of_interest & compar2 %in% cell_types_of_interest)
	)

for (downsample_proportion in downsample_proportions) {
	downsample_out_list = list()
	for (i in 1:nrow(comparison_grid)){
		compar1 = unname(comparison_grid$compar1[i])
		compar2 = unname(comparison_grid$compar2[i])

		temp_mat = input$mat
		temp_meta = input$metadata 

		if (args$paper == 'Luecken_2021') {
			temp_meta$label = gsub('\\/', '_', temp_meta$label)	
		}

		temp_meta %<>%
			filter(
				barcode %in% colnames(temp_mat),
				label %in% c(compar1, compar2)
			)
		temp_mat = temp_mat[,temp_meta$barcode]
		sampled_barcodes = temp_meta %>%
			group_by(
				label
			) %>%
			slice_sample(n = downsample_proportion) %>%
			pull(barcode)
		temp_meta %<>%
			filter(barcode %in% sampled_barcodes)
		temp_mat = temp_mat[,temp_meta$barcode]
		temp_mat = temp_mat[rowSums(temp_mat) > 0,]
		downsample_out_list[[paste0(compar1, '_vs_', compar2)]] = list('mat'=temp_mat,
			'metadata'=temp_meta)
	}
	final_output[[paste0('downsample_proportion=', downsample_proportion)]] = downsample_out_list
}

saveRDS(final_output, args$output_filename)