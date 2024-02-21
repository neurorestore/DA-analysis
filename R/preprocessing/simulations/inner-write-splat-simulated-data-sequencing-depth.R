setwd("~/git/DA-analysis")
options(stringsAsFactors = F)
library(argparse)

parser = ArgumentParser(prog = 'inner-write-splat-simulated-data-sequencing-depth.R')
grid = read.delim("sh/grids/preprocessing/write-splat-simulated-data-sequencing-depth.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

library(Seurat)
library(splatterBatch)
library(scater)
library(tidyverse)
library(magrittr)
library(Matrix)

# calculate group probabilities
group_probs = 1 / args$n_reps
# assign groups
unst = sample(seq(args$n_reps), args$n_reps/2)

params_dir = args$params_dir
params_files = list.files(params_dir, full.names=T)
params_files = params_files[!grepl('cell_filter', params_files)]
  
# split the vector by chunk number by specifying 
replicates_cells = split(1:args$n_cells, cut(seq_along(1:args$n_cells), length(params_files), labels = FALSE))

for (i in 1:length(params_files)){
	params_file = params_files[i]
	params = readRDS(params_file)

	# get sequencing depth from filename
	sequencing_depth = gsub('.*-downsample_proportion=', '', basename(params_file))
	sequencing_depth = gsub('-splat_estimate_params.rds', '', sequencing_depth)
	sequencing_depth = as.numeric(sequencing_depth)

	replicate_cells = replicates_cells[[i]]
	replicate_size = length(replicate_cells)

	# generate simulated cells
	sim = splatterBatch::splatSimulate(
		params = params,
		seed = args$simulation_iter,
		batchCells = replicate_size,
		de.prob = args$de_prob,
		de.facLoc = args$de_loc,
		verbose = T
	) %>% logNormCounts() %>% as.Seurat()

	sim %<>%
		RenameCells(
			new.names = paste0('Cell', replicate_cells)
		)

	sim@meta.data %<>%
		dplyr::mutate(cell_type = paste0('cell_', 1)) %>%
		mutate(
			label = ifelse(i %in% unst, 'unst', 'stim'),
			replicate = paste0('Depth=', sequencing_depth)
			) %>%
		set_rownames(colnames(GetAssayData(sim)))
	sim@meta.data$Cell = rownames(sim@meta.data)

	if (i==1){
		combined_sim = sim
	} else {
		combined_sim %<>% merge(., sim)
	}
}

saveRDS(combined_sim, args$output_filename)