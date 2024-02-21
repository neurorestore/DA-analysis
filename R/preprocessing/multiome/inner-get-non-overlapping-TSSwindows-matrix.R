setwd("~/git/DA-analysis")

library(tidyverse)
library(magrittr)
options(stringsAsFactors = FALSE)
library(argparse)
library(tidyverse)
library(magrittr)
library(Matrix)
library(dplyr)

parser = ArgumentParser(prog = 'inner-get-non-overlapping-TSSwindows-matrix.R')
grid = read.delim("sh/grids/preprocessing/get-non-overlapping-TSSwindows-matrix.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

source("R/functions/detect_system.R")
input = readRDS(args$input_filename)

gtf = read.delim(args$gtf_filename, header = FALSE, sep="\t",stringsAsFactors=FALSE, skip=1, fill=TRUE)
colnames(gtf) = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
gtf %<>%
	mutate(
		gene = gsub(".*gene_id (.*?);.*", "\\1", attributes),
	) %>%
	filter(gene %in% rownames(input$mat)) %>%
	dplyr::rename(chr=seqname) %>%
	dplyr::select(gene, chr, strand, start, end) %>%
	arrange(chr, strand, start, end)

curr_chr = ''
curr_strand = ''
filtered_genes_df = data.frame()
for (i in 1:nrow(gtf)){
	print(i)
	tmp_row = gtf[i,]
	chr = tmp_row$chr
	strand = tmp_row$strand
	start = tmp_row$start
	end = tmp_row$end

	tmp_row$idx = i

	if (chr != curr_chr) {
		prev_end = -1
		curr_chr = chr
	}

	if (strand != curr_strand) {
		prev_end = -1
		curr_strand = strand
	}

	if (start > prev_end){
		print('Adding new entry')
		filtered_genes_df %<>% rbind(., tmp_row)
		prev_add = TRUE
	}
	if (start <= prev_end & prev_add) {
		print('Remove previous entry')
		filtered_genes_df = filtered_genes_df[-nrow(filtered_genes_df),]
		prev_add = FALSE
	}
	prev_end = end
}

input$mat = input$mat[filtered_genes_df$gene,]
saveRDS(input, args$output_filename)