setwd("~/git/DA-analysis")
library(ArchR)
library(tidyverse)
library(magrittr)
library(argparse)

parser = ArgumentParser(prog = 'inner-get-ArchR-bg.R')
grid = read.delim("sh/grids/preprocessing/get-ArchR-bg.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

source('R/preprocessing/ArchR_functions.R')
metadata = readRDS(args$metadata_file)$metadata

proj = loadArchRProject(args$archR_dir)

archR_meta = data.frame(getCellColData(proj))
archR_meta$barcode = gsub('.*\\#', '', rownames(archR_meta))

comparisons_meta = archR_meta %>%
	left_join(
		metadata %>% select(barcode, label)
	)

comparison_grid = readRDS('metadata/comparisons.rds') %>%
	filter(paper == args$paper)

output_comparison_list = list()

for (i in 1:nrow(comparison_grid)){
	compar1 = comparison_grid$compar1[i]
	compar2 = comparison_grid$compar2[i]

	comparisons_meta0 = comparisons_meta %>%
		filter(label %in% c(compar1, compar2))

	label_counts = comparisons_meta0 %>%
		dplyr::select(label) %>%
		table()
	foreground_label = names(which.min(label_counts))
	background_label = names(label_counts)[names(label_counts) != foreground_label]

	match_bgd = matchBiasCellGroups(
		input = comparisons_meta0, 
		groups = comparisons_meta0$label,
		useGroups = foreground_label,
		bgdGroups = background_label,
		bias = c("TSSEnrichment", "log10(nFrags)"),
		k = 10,
		n = max(unname(label_counts))
	)

	curr_comparison = paste0(compar1, '_vs_', compar2)
	foreground_meta = comparisons_meta0[match_bgd@listData[[1]][[1]]$cells,] %>% as_tibble()
	background_meta = comparisons_meta0[match_bgd@listData[[1]][[1]]$bgd,] %>% as_tibble()
	out_list = list(
		'foreground_meta' = foreground_meta, 
		'background_meta' = background_meta
	)
 	output_comparison_list[[curr_comparison]] = out_list
}

saveRDS(output_comparison_list, args$output_filename)