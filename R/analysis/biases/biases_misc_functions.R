top_k_summary_helper = function(tmp_da_res, k){
	
	if (nrow(tmp_da_res) > 0){
		
		top_k_genes_by_mean = tmp_da_res %>%
			dplyr::slice(1:k) %>%
			summarise(
				mean_expr = mean(mean),
				sd_expr = mean(sd),
				percentage_open = mean(percentage_open),
				peak_width = mean(peak_width)
			) %>%
			mutate(summary_type = 'mean')
		return (top_k_genes_by_mean)
	}
}

stats_against_top_k_helper_function = function(output_df, da_res, mat, metadata, curr_cell_type, curr_paper_da_type, k = 1000, args = NULL) {
	
	da_types = da_res %>%
		filter(cell_type == curr_cell_type) %>%
		pull(de_type) %>%
		unique()

	for (da_type in da_types){
		if (curr_paper_da_type != 'cell type'){
			cell_type_and_label = str_split(curr_cell_type, '-')[[1]]
			actual_cell_type = cell_type_and_label[1]
			labels = str_split(cell_type_and_label[2], '_vs_')[[1]]
			tmp_cell_filter = tmp_metadata %>%
				filter(
					cell_type == actual_cell_type,
					label %in% labels, 
					barcode %in% colnames(tmp_mat),
					!is.na(replicate)
				) %>%
				pull(barcode)
		} else {
			labels = str_split(curr_cell_type, '_vs_')[[1]]
			if (args$paper == 'Luecken_2021'){
				labels = gsub('\\/', '_', labels)
			}
			tmp_cell_filter = tmp_metadata %>%
				filter(
					label %in% labels, 
					barcode %in% colnames(tmp_mat),
					!is.na(replicate)
				) %>%
				pull(barcode)
		}

		tmp_mat_2 = tmp_mat[,tmp_cell_filter]
		tmp_mat_2 = tmp_mat_2[rowSums(tmp_mat_2) > 0,]
		tmp_mat_2_bin = tmp_mat_2
		tmp_mat_2_bin@x[tmp_mat_2_bin@x > 0] = 1 

		tmp_gene_features = data.frame(
			'gene'=rownames(tmp_mat_2),
			'mean'=rowMeans(tmp_mat_2),
			'sd'=rowSds(tmp_mat_2),
			'percentage_open'=rowMeans(tmp_mat_2_bin)
		)

		if (args$feature == 'peak') {
			tmp_gene_features$peak_width = unlist(lapply(tmp_gene_features$gene, function(x){
				gene = gsub('-', '_', x)
				gene_coords = str_split(gene, '_')[[1]]
				gene_size = as.numeric(gene_coords[length(gene_coords)]) -
					as.numeric(gene_coords[length(gene_coords)-1])
			}))
		} else {
			tmp_gene_features$peak_width = 0
		}

		tmp_da_res = da_res %>%
			filter(de_type == da_type) %>% unique()
		tmp_da_res %<>%
			left_join(tmp_gene_features, by='gene') %>%
			filter(!is.na(mean)) %>%
			arrange(p_val, abs(test_statistic))

		out_row = top_k_summary_helper(tmp_da_res, k)
		output_df = rbind(output_df, out_row)
	}
	return (output_df)
}