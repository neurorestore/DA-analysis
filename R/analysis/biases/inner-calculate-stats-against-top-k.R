setwd("~/git/DA-analysis")

library(tidyverse)
library(magrittr)
options(stringsAsFactors = FALSE)
library(argparse)
library(Matrix)
library(stringr)
library(sparseMatrixStats)

parser = ArgumentParser(prog = 'inner-inner-calculate-stats-against-top-k.R')
grid = read.delim("sh/grids/analysis/inner-calculate-stats-against-top-k.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)
setwd(base_dir)

da_res = readRDS(args$concordance_filename)
tmp = readRDS(args$data_filename)
mat = tmp$mat
metadata = tmp$metadata
rm(tmp)
gc()

output_df = data.frame()
grid_to_run = da_res %>% 
	select(cell_type, paper_da_type) %>%
	distinct()

for (i in 1:nrow(grid_to_run)){
	print(i)
	curr_cell_type = inner_grid_to_run$cell_type[i]
	curr_paper_da_type = inner_grid_to_run$paper_da_type[i]

	print(paste0('curr cell type: ', curr_cell_type, '; curr paper da type: ', curr_paper_da_type))
	da_res_2 = da_res %>% filter(cell_type == curr_cell_type, paper_da_type == curr_paper_da_type)

	output_df = stats_against_top_k_helper_function(output_df, da_res_2, mat, metadata, curr_cell_type, curr_paper_da_type, k = args$k, args)
}

if (nrow(output_df) > 0) {
	saveRDS(output_df, args$output_filename)
}