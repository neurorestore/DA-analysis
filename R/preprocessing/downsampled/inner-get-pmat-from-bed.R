setwd("~/git/DA-analysis")

library(argparse)
library(stringr)
library(tidyverse)
library(magrittr)
library(Seurat)
library(reshape2)
library(Matrix)
library(stringi)
library(GenomicRanges)
library(plyr)

parser = ArgumentParser(prog = 'inner-get-pmat-from-bed.R')
grid = read.delim("sh/grids/preprocessing/get-pmat-from-bed.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

peaks_df = read.table(args$peak_filename)
colnames(peaks_df) = c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'FC', 'log_q', 'log_p', 'summit_pos')
peaks_df %<>%
	mutate(
		chrom = sub('\\.', '', chrom),
		chrom = paste0('chr', chrom)
	)
# somehow mutate doesn't work
peaks_df$chrom = gsub("b'","",peaks_df$chrom)
peaks_df$chrom = gsub("'","",peaks_df$chrom)
peaks_df %<>%
	mutate(
		peak_id = paste0(chrom, '_', chromStart, '_', chromEnd)
	)

peaks_granges = GenomicRanges::GRanges(peaks_df$chrom, IRanges(peaks_df$chromStart, peaks_df$chromEnd))

fragments_df = read.table(gzfile(args$bed_filename), header=T)
fragments_df %<>%
	mutate(barcode = name)

fragments_granges = GRanges(fragments_df$chrom, 
				   IRanges(fragments_df$chromStart, fragments_df$chromEnd), 
				   barcode=fragments_df$barcode
)

ov = findOverlaps(fragments_granges, peaks_granges)

ind = data.frame('idy'=fragments_df$barcode[queryHits(ov)], 'idx'=peaks_df$peak_id[subjectHits(ov)])
count_df = plyr::count(ind, vars = c("idx", "idy"))
colnames(count_df) = c('gene', 'barcode', 'freq')

mat = xtabs(freq~gene+barcode, count_df, sparse=T)
metadata = data.frame('barcode'=colnames(mat)) %>%
	mutate(
		'replicate'=gsub('_Cell.*', '', barcode),
		'cell_type'=paste0('cell_1')
		)
labels = unique(metadata$replicate)
stim_labels = sample(labels, floor(length(labels)/2))

metadata %<>% 
	mutate(
		label = ifelse(replicate %in% stim_labels, 'stim', 'unst')
	)
rownames(metadata) = metadata$barcode
mat = mat[,metadata$barcode]
sc = CreateSeuratObject(mat, meta.data=metadata)
saveRDS(sc, args$output_filename)