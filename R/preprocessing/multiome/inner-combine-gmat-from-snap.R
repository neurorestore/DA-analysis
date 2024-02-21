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

parser = ArgumentParser(prog = 'inner-combine-gmat-from-snap.R')
grid = read.delim("sh/grids/preprocessing/combine-gmat-from-snap.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

filenames = list.files(args$snap_dir)
filenames = filenames[filenames != args$snap_dir]
filenames = paste0(args$snap_dir, filenames)
filenames = filenames[endsWith(filenames, '.snap')]
rna_tmp = readRDS(args$rna_data_file)
metadata = rna_tmp$metadata

if (args$cellranger_barcodes == 1){
	barcodes_df = read.csv('metadata/cellranger_arc_barcodes_mapping.csv',row.names=1)
}

combined_mat = NULL
barcodes_mapped_threshold = 0.8

for (snap_input_file in filenames){
	print(snap_input_file)
	h5_file = H5File$new(snap_input_file, mode = "r")
	barcodes = readDataSet(h5_file[['BD']][['name']])
	sample = gsub('\\.snap', '', basename(snap_input_file))

	if (args$paper == 'Boukhaled_2022'){
		barcodes = as.character(unlist(sapply(barcodes, function(x){substr(x, 9, nchar(x))})))
	}

	if (args$paper == 'Squair_2022'){
		barcodes = gsub('-1', '', barcodes)
		barcodes = paste0(sample, '-', barcodes)
	}

	barcodes_mapped = TRUE
	sample_barcode = FALSE

	temp_meta = metadata %>%
		filter(replicate==!!sample)

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

	## only process if barcodes are mapped properly
	if (barcodes_mapped == TRUE){
		
		gtf = read.delim(args$gtf_filename, header = FALSE, sep="\t",stringsAsFactors=FALSE, skip=1, fill=TRUE)
		colnames(gtf) = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame",
						  "attributes")
		gene_ids =  gsub(".*gene_id (.*?);.*", "\\1", gtf$attributes)
		transcript_ids = gsub(".*transcript_id (.*?);.*", "\\1", gtf$attributes)
		index=gene_ids!=""
		gene_granges = GRanges(seqnames=gtf$seqname[index],
							   ranges=IRanges(gtf$start[index],gtf$end[index]),
							   strand=gtf$strand[index],
							   tx_id=transcript_ids[index],
							   gene_id=gene_ids[index])
		gene_granges$name = gene_ids
		
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
		
		ov = findOverlaps(fragments_granges, gene_granges)
		
		ind = data.frame('idy'=fragments_df$barcode[queryHits(ov)], 'idx'=gene_ids[subjectHits(ov)])
		count_df = plyr::count(ind, vars = c("idx", "idy"))
		gmat = xtabs(freq~idx+idy, count_df, sparse=T)

		if (args$cellranger_barcodes==1) {
			temp_meta %<>%
				filter(atac_barcode %in% colnames(gmat))
			gmat = gmat[, temp_meta$atac_barcode]
			colnames(gmat) = temp_meta$barcode
		} else {
			temp_meta %<>%
				filter(barcode %in% colnames(gmat))	
			gmat = gmat[, temp_meta$barcode]
		}

		if (is.null(combined_mat)){
			combined_mat = gmat
		} else {
			if (any(!rownames(gmat) %in% rownames(combined_mat))){
				cbind_mat = gmat[rownames(gmat) %in% rownames(combined_mat),]
				ori_rownames = rownames(cbind_mat)
				cbind_mat_missing_genes = rownames(combined_mat)[!rownames(combined_mat) %in% rownames(gmat)]
				cbind_mat = rbind(cbind_mat, 
					matrix(0, nrow=length(cbind_mat_missing_genes), ncol=ncol(gmat)))
				rownames(cbind_mat) = c(ori_rownames, cbind_mat_missing_genes)
				cbind_mat = cbind_mat[rownames(combined_mat),]
				combined_mat = cbind(combined_mat, cbind_mat)

				rbind_mat_missing_genes = rownames(gmat)[!rownames(gmat) %in% rownames(combined_mat)]
				ori_colnames = colnames(combined_mat)[!colnames(combined_mat) %in% colnames(gmat)]
				rbind_mat = cbind(matrix(0, nrow=length(rbind_mat_missing_genes), ncol=length(ori_colnames)), 
					gmat[rbind_mat_missing_genes,])
				colnames(rbind_mat) = c(ori_colnames, colnames(gmat))
				rbind_mat = rbind_mat[,colnames(combined_mat)]
				rownames(rbind_mat) = rbind_mat_missing_genes
				combined_mat = rbind(combined_mat, rbind_mat)
			} else{
				gmat = gmat[rownames(combined_mat),]
				combined_mat = cbind(combined_mat, gmat)
			}
			print(dim(combined_mat))
		}
	}
}

## finally create list of metadata and matrix
metadata = metadata[metadata$barcode %in% colnames(combined_mat),]
combined_mat = combined_mat[,metadata$barcode]
saveRDS(list('mat'=combined_mat, 'metadata'=metadata), args$output_filename)