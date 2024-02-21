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
library(peakRAM)

parser = ArgumentParser(prog = 'inner-run-scalability-by-comparison.R')
grid = read.delim("sh/grids/analysis/run-scalability-by-comparison.txt")
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
da_types = da_family_args[[args$da_family]][[args$da_method]]

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

	output_row = peakRAM(
		run_da_wrapper(data.frame(), temp_mat_2, temp_meta_2, args, gc_content_df = gc_content_df, da_types = da_types)
	)
	output_df %<>% rbind(output_row)
}
	
if (nrow(output_df) > 0){
	saveRDS(output_df, args$output_filename)	
}