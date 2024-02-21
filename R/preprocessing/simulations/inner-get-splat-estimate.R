setwd("~/git/DA-analysis")
options(stringsAsFactors = FALSE)
library(argparse)
library(tidyverse)
library(magrittr)
library(Matrix)

parser = ArgumentParser(prog = 'inner-get-splat-estimate.R')
grid = read.delim("sh/grids/preprocessing/get-splat-estimate.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

source('R/preprocessing/simulations/splat-estimate-functions.R')

counts = readRDS(args$input_filename)$mat
cell_filter = sample(colnames(counts), floor(ncol(counts)*(as.numeric(args$random_sample_proportion))))
counts = counts[,cell_filter]
counts = counts[rowSums(counts) > 0,]
params = splatEstimate.matrix(counts)

saveRDS(params, args$output_filename)