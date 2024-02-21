setwd("~/git/DA-analysis")
options(stringsAsFactors = F)
library(argparse)

parser = ArgumentParser(prog = 'inner-write-splat-simulated-data.R')
grid = read.delim("sh/grids/preprocessing/write-splat-simulated-data.txt")
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

params = readRDS(args$params_file)
# calculate group probabilities
group_probs = 1 / args$n_reps
# assign groups
unst = sample(seq(args$n_reps), args$n_reps/2)

# generate simulated cells
sim = splatterBatch::splatSimulateGroups(
	params = params,
	seed = args$simulation_iter,
	batchCells = args$n_cells,
	de.prob = args$de_prob,
	de.facLoc = args$de_loc,
	group.prob = rep(group_probs, args$n_reps), verbose = T
) %>% logNormCounts() %>% as.Seurat()

# adjust metadata for default input to Seurat
sim@meta.data %<>%
	dplyr::mutate(cell_type = paste0('cell_', 1)) %>%
	dplyr::rename(label = Group) %>%
	mutate(replicate = gsub("Group", "Replicate ", label)) %>%
	mutate(label = as.numeric(gsub("Group", "", label))) %>%
	mutate(label = ifelse(label %in% unst, 'unst', 'stim')) %>%
	set_rownames(colnames(GetAssayData(sim)))

saveRDS(sim, args$output_filename)