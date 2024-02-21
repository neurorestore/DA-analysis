library(tidyverse)
library(magrittr)
library(Matrix)
library(sparseMatrixStats)

output_dir = 'da_analysis_false_discoveries/'
data_dir = paste0(output_dir, 'by_comparison/')

compar_files = list.files(data_dir, full.names=T)
tmp = readRDS('/multiome/Luecken_2021_peak_matrix.rds')
metadata = tmp$metadata
mat = tmp$mat

peak_sizes = readRDS('peaks/sc_peaks/Luecken_2021/final_peaks.rds')

curr_base_file = ''
sig_summary_df = data.frame()

for (compar_file in compar_files){
	temp_df = readRDS(compar_file)
	temp_params = temp_df %>% 
		dplyr::select(filename, seq_type, pseudobulk, min_cell_filter, pseudo_repl, binarization, da_family, da_method, label, latent_vars, normalization, percent_filter) %>% 
		distinct()
	
	temp_meta = metadata %>%
		filter(label == temp_params$label)
	temp_mat = mat[,temp_meta$barcode]
	temp_mat = temp_mat[rowSums(temp_mat) > 0,]
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
						percetile_type = var_name
					) 
				sig_summary_df %<>% rbind(curr_temp_params)
			}
		}
	}
}
saveRDS(sig_summary_df, 'data/summaries/false_discoveries/false_discoveries_luecken.rds')