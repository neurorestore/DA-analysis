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
library(SummarizedExperiment)
library(qsmooth)

parser = ArgumentParser(prog = 'inner-run-da-by-comparison-simulated.R')
grid = read.delim("sh/grids/analysis/run-da-by-comparison-simulated.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

source('R/analysis/differential_accessibility/da_misc_functions.R')

sc = readRDS(args$input_filename)
methods_list = da_family_args[[args$da_family]]

## Libra take counts
if (args$binarization == 't'){
	temp_mat = GetAssayData(sc, slot='counts')
	temp_mat@x[temp_mat@x > 0] = 1
	sc %<>% SetAssayData(temp_mat, slot='counts')
}

temp_mat = GetAssayData(sc, slot='counts')
temp_meta = sc@meta.data

## only for downsampled data
if (args$normalization == 'qsmooth'){
	seed = gsub('.*-seed=', 'seed=', basename(args$input_filename))
	seed = gsub('\\.rds', '', seed)
	gc_content_df = readRDS('peaks/GC_content/downsampled_simulated_peaks_gc_content.rds')[[seed]]
	rownames(gc_content_df) = gc_content_df$gene
} else {
	gc_content_df = NULL
}

output_df = run_da_wrapper(temp_mat, temp_meta, args, gc_content_df = gc_content_df, da_types = methods_list)
if (nrow(output_df) > 0){
	saveRDS(output_df, args$output_filename)
}
