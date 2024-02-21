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

parser = ArgumentParser(prog = 'inner-combine-pmat-from-snap.R')
grid = read.delim("sh/grids/preprocessing/combine-pmat-from-snap.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

filenames = list.files(args$snap_dir)
filenames = filenames[filenames != args$snap_dir]
filenames = paste0(args$snap_dir, filenames)
filenames = filenames[endsWith(filenames, '.snap')]
metadata = readRDS(args$metadata_filename)$metadata

# check if barcodes are cellranger-based
if (args$cellranger_barcodes == 1){
	barcodes_df = read.csv('metadata/cellranger_arc_barcodes_mapping.csv',row.names=1)
}

full_count_df = data.frame()
barcodes_mapped_threshold = 0.8

for (snap_input_file in filenames){
	print(snap_input_file)
	h5_file = H5File$new(snap_input_file, mode = "r")
	barcodes = readDataSet(h5_file[['BD']][['name']])
	sample = gsub('\\.snap', '', basename(snap_input_file))
	
	## custom preprocessing steps
	if (args$paper == 'Squair_2022'){
		barcodes = gsub('-1', '', barcodes)
		barcodes = paste0(sample, '-', barcodes)
	}

	if (args$paper == 'Luecken_2021'){
		barcodes = gsub('-1', '', barcodes)
		sample = gsub('site', 's', sample)
		sample = gsub('_', '', sample)
		sample = gsub('donor', 'd', sample)
		if (!grepl('10', sample)){
			sample = gsub('0', '', sample)
		}
		barcodes = paste0(barcodes, '-', sample)
	}

	if (args$file_prefix == 'Pliner_2018_scATAC'){
		barcodes = paste0(sample, '_', barcodes)
	}

	barcodes_mapped = TRUE
	sample_barcode = FALSE

	## if data is per-sample or combined-sample
	if (args$by_sample == 1){
		temp_meta = metadata %>%
			filter(replicate==!!sample)
	} else {
		temp_meta = metadata
	}

	if (any(grepl('_', barcodes))){
		ori_barcodes = data.frame(do.call(rbind, lapply(barcodes, function(x){
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
					if (mean(barcodes %in% barcodes_df[,barcode_name]) > 0.8){
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

				barcodes = paste0(sample, '_', barcodes)
			}
		}
	}

	## only process if majority of barcodes can be mapped
	if (barcodes_mapped == TRUE){
		## get peaks
		peaks_df = read.table(args$peak_filename) %>% 
			set_colnames(c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'FC', 'log_q', 'log_p', 'summit_pos')) %>%
			mutate(
				chrom = sub('.', '', chrom)
			)
		peaks_df$chrom = gsub("'","",peaks_df$chrom)
		peaks_df %<>%
			mutate(
				peak_id = paste0(chrom, '_', chromStart, '_', chromEnd)
			)

		peaks_granges = GenomicRanges::GRanges(peaks_df$chrom, IRanges(peaks_df$chromStart, peaks_df$chromEnd))

		# filter out overlapping granges from qval=0.01
		if (args$qval != 0.01){
			clean_peaks_filename = gsub('-qval=.*', '_peaks.narrowPeak', args$peak_filename)
			clean_peaks_df = read.table(clean_peaks_filename) %>% 
				set_colnames(c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'FC', 'log_q', 'log_p', 'summit_pos')) %>% 
				mutate(
					chrom = sub('.', '', chrom)
				)
			clean_peaks_df$chrom = gsub("'","",clean_peaks_df$chrom)
			clean_peaks_df %<>%
				mutate(
					peak_id = paste0(chrom, '_', chromStart, '_', chromEnd)
				)

			clean_peaks_granges = GenomicRanges::GRanges(clean_peaks_df$chrom, IRanges(clean_peaks_df$chromStart, clean_peaks_df$chromEnd))
			peaks_granges =  subsetByOverlaps(peaks_granges, clean_peaks_granges, invert = TRUE)
			rm(clean_peaks_df, clean_peaks_granges)
			gc()
		}
		
		barcodes_idx = data.frame(
			'barcode'=barcodes,
			'start_idx'=readDataSet(h5_file[['FM']][['barcodePos']]),
			'len'=readDataSet(h5_file[['FM']][['barcodeLen']])
			)
		
		fragments_df = data.frame(
			'chrom' = readDataSet(h5_file[['FM']][['fragChrom']]),
			'start' = readDataSet(h5_file[['FM']][['fragStart']]),
			'len' = readDataSet(h5_file[['FM']][['fragLen']])
			) %>% 
			mutate(end = start + len - 1) %>%
			relocate(end, .after=start)
		
		fragments_df$barcode = unlist(sapply(1:nrow(barcodes_idx), function(x){
			return (rep(barcodes_idx$barcode[x], barcodes_idx$len[x]))
		})) 

		if (args$cellranger_barcodes == 1){
			fragments_df %<>%
				filter(barcode %in% temp_meta$atac_barcode)
		} else {
			fragments_df %<>%
				filter(barcode %in% temp_meta$barcode)
		}

		fragments_granges = GRanges(fragments_df$chrom, 
						   IRanges(fragments_df$start, fragments_df$end), 
						   barcode=fragments_df$barcode
		)
		
		ov = findOverlaps(fragments_granges, peaks_granges)
		ind = data.frame('idy'=fragments_df$barcode[queryHits(ov)], 'idx'=peaks_df$peak_id[subjectHits(ov)])
		count_df = plyr::count(ind, vars = c("idx", "idy"))
		colnames(count_df) = c('gene', 'barcode', 'freq')
		if (args$cellranger_barcodes==1) {
			temp_meta %<>%
				filter(atac_barcode %in% count_df$barcode)
			count_df %<>%
				filter(barcode %in% temp_meta$atac_barcode) %>%
				dplyr::rename(atac_barcode = barcode) %>%
				left_join(temp_meta[,c('atac_barcode', 'barcode')]) %>%
				select(gene, barcode, freq)

		} else {
			temp_meta %<>%
				filter(barcode %in% count_df$barcode)
			count_df %<>%
				filter(barcode %in% temp_meta$barcode)
		}

		full_count_df = rbind(full_count_df, count_df)
	}
}

## finally create list of metadata and matrix
combined_mat = xtabs(freq~gene+barcode, full_count_df, sparse=T)
metadata = metadata[metadata$barcode %in% colnames(combined_mat),]
combined_mat = combined_mat[,metadata$barcode]
saveRDS(list('mat'=combined_mat, 'metadata'=metadata), args$output_filename)