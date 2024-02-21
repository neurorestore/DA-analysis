setwd("~/git/DA-analysis")

library(tidyverse)
library(magrittr)
options(stringsAsFactors = FALSE)
library(Matrix)
library(stringr)


outer_data_dirs = 'data/da_analysis_simulated'
data_types = c('peaks', 'da_analysis_simulated/peaks_with_sequencing_depth', 'downsampled')
output_dir = 'data/summaries/false_discoveries/'

for (data_type ind data_types) {
	filenames = list.files(file.path(outer_data_dir, data_type), full.names=T)
	sig_summary_df = data.frame()

	for (filename in filenames) {
		temp_df = readRDS(filename)
		temp_params = temp_df %>% 
			dplyr::select(filename, seq_type, pseudobulk, min_cell_filter, pseudo_repl, binarization, da_family, da_method, latent_vars, normalization, percent_filter) %>% 
			distinct()
		
		seurat_filename = gsub('-pseudobulk=.*', '.rds', filename)
		sc = readRDS(seurat_filename)
		temp_mat = GetAssayData(sc, slot='counts')
		bin_temp_mat = temp_mat
		bin_temp_mat@x[bin_temp_mat > 0] =1

		count_meta = data.frame(
			'mean'=rowMeans(temp_mat),
			'percentage_open'=rowMeans(bin_temp_mat)
		) %>%
		left_join(peak_sizes) %>%
		set_rownames(NULL)

		for (da_type in unique(temp_df$da_type)){
			curr_temp_df = temp_df %>%
				filter(da_type == !!da_type) 
			number_of_regions = nrow(curr_temp_df)
			
			curr_temp_df %<>%
				filter(!is.na(p_val_adj), !is.na(avg_logFC)) %>%
				left_join(peak_sizes, by=c('gene'='peaks')) %>%
				left_join(count_meta, by='gene')

			sig_genes_by_pval = sum(curr_temp_df$p_val < 0.05)
			sig_genes_by_pval_adj = sum(curr_temp_df$p_val_adj < 0.05)
			number_of_regions = nrow(curr_temp_df)

			curr_temp_params = temp_params %>%
				mutate(
					da_type = da_type,
					pval_cutoff=0.05,
					l2fc_cutoff=0,
					number_of_regions=number_of_regions,
					number_of_da_regions=sig_genes_by_pval_adj,
					number_of_nominal_da_regions=sig_genes_by_pval,
					percentile = NA,
					percetile_type = NA
				) 
			sig_summary_df %<>% rbind(curr_temp_params)
			
			# get decile-level false discoveries
			var_names = c('mean', 'percentage_open', 'peak_size')
			for (var_name in var_names) {
				curr_temp_df %<>% mutate(percentile = ntile(get(var_name), 10))
				for (i in 1:10) {

					curr_temp_df1 = curr_temp_df %>%
						filter(percentile == i)

					sig_genes_by_pval = sum(curr_temp_df1$p_val < 0.05)
					sig_genes_by_pval_adj = sum(curr_temp_df1$p_val_adj < 0.05)
					number_of_regions = nrow(curr_temp_df1)

					curr_temp_params = temp_params %>%
						mutate(
							da_type = da_type,
							pval_cutoff=0.05,
							l2fc_cutoff=0,
							number_of_regions=number_of_regions,
							number_of_da_regions=sig_genes_by_pval_adj,
							number_of_nominal_da_regions=sig_genes_by_pval,
							percentile = paste0('percentile_', i),
							percetile_type = 'var_name'
						) 
					sig_summary_df %<>% rbind(curr_temp_params)
				}
			}
		}
	}
	output_filename = paste0(output_dir, 'false_discoveries_', data_type, '.rds')
	saveRDS(sig_summary_df, output_filename)	
}
