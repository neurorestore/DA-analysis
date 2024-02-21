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
library(doParallel)
library(parallel)
library(foreach)

parser = ArgumentParser(prog = 'inner-get-fragments-bed-from-snap.R')
grid = read.delim("sh/grids/preprocessing/get-fragments-bed-from-snap.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

if (!dir.exists(args$output_dir)){
	dir.create(args$output_dir)
}

buffer_size = 500

temp = readRDS(args$metadata_file)
mat = temp$mat
metadata = temp$metadata
rownames(metadata) = metadata$barcode
mat = mat[,colnames(mat) %in% rownames(metadata)]
metadata = metadata[colnames(mat),]
rm(temp)
rm(mat)

if (args$cellranger_barcodes == 1){
	barcodes_df = read.csv('metadata/cellranger_arc_barcodes_mapping.csv',row.names=1)
}

barcodes_mapped_threshold = 0.8

snap_input_file = args$snap_input_file
h5_file = H5File$new(snap_input_file, mode = "r")
barcodes_idx_df = data.frame('ori_barcode'=readDataSet(h5_file[['BD']][['name']])) %>%
	mutate(barcode=ori_barcode)
sample = gsub('\\.snap', '', basename(snap_input_file))

if (args$paper == 'Boukhaled_2022'){
	barcodes_idx_df$barcode = as.character(unlist(sapply(barcodes_idx_df$barcode, function(x){substr(x, 9, nchar(x))})))
}

if (args$paper == 'Squair_2022'){
	barcodes_idx_df$barcode = gsub('-1', '', barcodes_idx_df$barcode)
	barcodes_idx_df$barcode = paste0(sample, '-', barcodes_idx_df$barcode)
}

if (args$paper == 'Luecken_2021'){
	barcodes_idx_df$barcode = gsub('-1', '', barcodes_idx_df$barcode)
	sample = gsub('site', 's', sample)
	sample = gsub('_', '', sample)
	sample = gsub('donor', 'd', sample)
	if (!grepl('10', sample)){
		sample = gsub('0', '', sample)
	}
	barcodes_idx_df$barcode = paste0(barcodes_idx_df$barcode, '-', sample)
}

barcodes_mapped = TRUE
sample_barcode = FALSE

if (args$by_sample == 1){
	temp_meta = metadata %>%
		filter(replicate==!!sample)
} else {
	temp_meta = metadata
}

if (nrow(temp_meta) > 0){
	if (any(grepl('_', barcodes_idx_df$barcode))){
		ori_barcodes = data.frame(do.call(rbind, lapply(barcodes_idx_df$barcode, function(x){
		stri_reverse(stri_split_fixed(stri_reverse(x),"_",n = 2)[[1]])[2:1]
		})))
		colnames(ori_barcodes) = c('sample', 'barcode')
		ori_barcodes %<>%
			mutate(full_barcode = paste0(sample, '_', barcode))
		sample_barcode = TRUE
		if (args$cellranger_barcodes == 1){
			barcodes_mapped = FALSE
			true_barcode_name = 'atac_barcode'
			for (barcode_name in colnames(barcodes_df)){
					if (mean(ori_barcodes$barcode %in% barcodes_df[,barcode_name]) > barcodes_mapped_threshold){
							barcodes_mapped = TRUE
							true_barcode_name = barcode_name
					}
			}
			if (barcodes_mapped == TRUE){
				temp_barcodes_df = barcodes_df[,c('rna_barcode', true_barcode_name)]
				colnames(temp_barcodes_df)[2] = 'atac_barcode'
				colnames(ori_barcodes)[2] = 'atac_barcode'
				temp_barcodes_df %<>% 
					inner_join(ori_barcodes, by='atac_barcode') %>%
					filter(!is.na(rna_barcode)) %>%
					mutate(rna_barcode=paste0(sample, '_', rna_barcode)) %>%
					select(-sample, -atac_barcode) %>%
					dplyr::rename(atac_barcode=full_barcode)

				if (!'rna_barcode' %in% colnames(temp_meta)){
					temp_meta %<>% mutate(rna_barcode = barcode)
				}
				temp_meta %<>% left_join(temp_barcodes_df)
			}
		}
	} else{
		if (args$cellranger_barcodes == 1){
			barcodes_mapped = FALSE
			true_barcode_name = 'atac_barcode'
			for (barcode_name in colnames(barcodes_df)){
					if (mean(barcodes_idx_df$barcode %in% barcodes_df[,barcode_name]) > 0.8){
							barcodes_mapped = TRUE
							true_barcode_name = barcode_name
					}
			}
			if (barcodes_mapped == TRUE){
				temp_barcodes_df = barcodes_df[,c('rna_barcode', true_barcode_name)]
				colnames(temp_barcodes_df)[2] = 'atac_barcode'
				if (!'rna_barcode' %in% colnames(temp_meta)){
					temp_meta %<>% mutate(rna_barcode = barcode)
				}
				if (grepl('_', temp_meta$rna_barcode[1])){
					temp_barcodes_df %<>% mutate(rna_barcode = paste0(sample, '_', rna_barcode))
				}

				temp_meta %<>% left_join(temp_barcodes_df) %>%
					mutate(atac_barcode = paste0(replicate, '_', atac_barcode))

				barcodes_idx_df$barcode = paste0(sample, '_', barcodes_idx_df$barcode)
			}
		}
	}

	if (barcodes_mapped == TRUE){

		barcodes_output_filename = paste0(args$output_dir, '/', sample, '_barcodes.txt')

		if (args$paper == 'Buenrostro_2018' & grepl('bulk_ATAC', sample)){
			barcodes_output_filename = gsub('Corces_2016', 'Buenrostro_2018', barcodes_output_filename)
		} 

		if (args$cellranger_barcodes == 1){
			selected_barcodes = barcodes_idx_df %>%
				filter(barcode %in% temp_meta$atac_barcode) %>%
				pull(ori_barcode)
		} else {
			selected_barcodes = barcodes_idx_df %>%
				filter(barcode %in% temp_meta$barcode) %>%
				pull(ori_barcode)
		}
		write.table(selected_barcodes, file = barcodes_output_filename, append = FALSE, quote = FALSE, sep = "\t",
			eol = "\n", na = "NA", dec = ".", row.names = FALSE,
			col.names = FALSE, qmethod = c("escape", "double"),
			fileEncoding = "")
		flag = system2(command=args$snaptools_path, 
			args=c("dump-fragment", 
				"--snap-file", snap_input_file, 
				"--output-file", args$bed_output_filename, 
				"--barcode-file", barcodes_output_filename,
				"--buffer-size", buffer_size
				)
		)
	}
}