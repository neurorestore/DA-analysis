setwd("~/git/DA-analysis")

library(tidyverse)
library(magrittr)
options(stringsAsFactors = FALSE)
library(argparse)
library(tidyverse)
library(magrittr)
library(Matrix)
library(dplyr)
library(GenomicRanges)

parser = ArgumentParser(prog = 'inner-get-enhancer-promoter-matrix.R')
grid = read.delim("sh/grids/preprocessing/get-enhancer-promoter-matrix.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name), type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

input = readRDS(args$input_filename)

annotation_files = c(
	'enhancer' = 'GRCh38-ELS.bed', 
	'promoter' = 'GRCh38-PLS.bed'
)

enhancer_df = read.delim(paste0(args$metadata_dir, 'GRCh38-ELS.bed'), header = FALSE) %>% set_colnames(c('chrom', 'chromStart', 'chromEnd', 'name', 'details1', 'details2'))
enhancer_granges = GenomicRanges::GRanges(enhancer_df$chrom, IRanges(enhancer_df$chromStart, enhancer_df$chromEnd))

promoter_df = read.delim(paste0(args$metadata_dir, 'GRCh38-PLS.bed'), header = FALSE) %>% set_colnames(c('chrom', 'chromStart', 'chromEnd', 'name', 'details1', 'details2'))
promoter_granges = GenomicRanges::GRanges(promoter_df$chrom, IRanges(promoter_df$chromStart, promoter_df$chromEnd))

peaks_df = data.frame(do.call(rbind, lapply(rownames(input$mat), function(x){
	chrom_details = strsplit(x, '_')[[1]]
	chrom_pos = chrom_details[(length(chrom_details)-1):length(chrom_details)]
	chrom = paste0(chrom_details[1:(length(chrom_details)-2)], collapse='_')
	c(chrom, chrom_pos) 
}))) %>%
	set_colnames(c('chrom', 'chromStart', 'chromEnd')) %>%
	type_convert()
peaks_granges = GenomicRanges::GRanges(peaks_df$chrom, IRanges(peaks_df$chromStart, peaks_df$chromEnd))

enhancer_hits =  countOverlaps(peaks_granges, enhancer_granges)
promoter_hits =  countOverlaps(peaks_granges, promoter_granges)

if (args$annotation == 'enhancer') {
	hits = promoter_hits == 0 & enhancer_hits != 0
} else if (args$annotation == 'promoter') {
	hits = enhancer_hits == 0 & promoter_hits != 0
}

input$mat = input$mat[hits,]

saveRDS(input, args$output_filename)