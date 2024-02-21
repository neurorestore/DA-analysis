experimental_design = function(bulk_da, sc_da, experiement) {
	output_list = list('none'=list('bulk_da'=bulk_da, 'sc_da'=sc_da))

	# multiome, multiome go -- lowly expressed genes
	if (experiement == 'non_overlapping_windows') {
		sc_da %<>% filter(!overlapping_windows)
		bulk_da %<>% filter(gene %in% sc_da$gene)
		output_list = list('non_overlapping_windows' = list('bulk_da'=bulk_da, 'sc_da'=sc_da))
	}

	# multiome, multiome go -- lowly expressed genes
	if (experiement == 'lowly_expressed_genes') {
		bulk_da %<>% 
			mutate(percentile = ntile(mean, 3)) %>%
			filter(percentile == 1)
		sc_da %<>% filter(gene %in% bulk_da$gene)
		output_list = list('lowly_expressed_genes' = list('bulk_da'=bulk_da, 'sc_da'=sc_da))
	}

	# matched bulk, multiome -- logFC filtering
	if (experiment == 'logFC_filtering') {
		output_list = list()
		quantile_filters = quantile(abs(temp_sc_da$avg_logFC), seq(0.1, 0.9, 0.1))

		for (i in 1:length(quantile_filters)) {
			quantile_filter = unname(quantile_filters[i])
			curr_filter = paste0('logFC_filter_', names(quantile_filters)[i])					
			temp_sc_da = temp_sc_da %>%
				mutate(abs_avg_logFC = abs(avg_logFC)) %>%
				filter(abs_avg_logFC > quantile_filter) %>%
				dplyr::select(-abs_avg_logFC)
			temp_bulk_da = bulk_da %>% filter(gene %in% temp_sc_da_gene)

			output_list[[curr_filter]] = list('bulk_da'=temp_bulk_da, 'sc_da'=temp_sc_da)
			
			curr_filter = paste0(curr_filter, '_pseudo')
			pseudo_size = sum(abs(temp_sc_da$avg_logFC) > quantile_filter)

			temp_sc_da_2 = temp_sc_da %>%
				sample_n(pseudo_size)
			temp_bulk_da = bulk_da %>% filter(gene %in% temp_sc_da_gene)

			output_list[[curr_filter]] = list('bulk_da'=temp_bulk_da, 'sc_da'=temp_sc_da)
		}
	}

	return (output_list)
}

calculate_concordance = function(output_df, bulk_da, sc_da, k=500, filter='NA'){
	bulk_size = nrow(bulk_da)
	sc_size = nrow(sc_da)

	duplicated_bulk_genes = bulk_da$gene[duplicated(bulk_da$gene)]
	duplicated_sc_genes = sc_da$gene[duplicated(sc_da$gene)]

	sc_da %<>% filter(!is.na(p_val), !is.na(p_val_adj), !gene %in% duplicated_sc_genes) %<>% unique()
	bulk_da %<>% filter(!is.na(p_val), !is.na(p_val_adj), !gene %in% duplicated_bulk_genes) %<>% unique()
	
	# replace p=0 with minimum p-value
	sc_da_min = min(sc_da$p_val_adj[sc_da$p_val_adj > 0])
	bulk_da_min = min(bulk_da$p_val_adj[bulk_da$p_val_adj > 0])
	sc_da %<>% 
		mutate(p_val_adj = ifelse(p_val_adj <= sc_da_min, sc_da_min, p_val_adj))
	bulk_da %<>% 
		mutate(p_val_adj = ifelse(p_val_adj <= bulk_da_min, bulk_da_min, p_val_adj))
	## repeat for raw p-values
	sc_da_min = min(sc_da$p_val[sc_da$p_val > 0])
	bulk_da_min = min(bulk_da$p_val[bulk_da$p_val > 0])
	sc_da %<>% 
		mutate(p_val = ifelse(p_val <= sc_da_min, sc_da_min, p_val))
	bulk_da %<>%
		mutate(p_val = ifelse(p_val <= bulk_da_min, bulk_da_min, p_val))
	
	# filter to genes detected in both single-cell and bulk data
	genes = intersect(bulk_da$gene, sc_da$gene)
	sc_da %<>% filter(gene %in% genes) %>% arrange(gene)
	bulk_da %<>% filter(gene %in% genes) %>% arrange(gene)
	
	k = as.integer(k)
	
	vec1 = bulk_da %>%
		arrange(p_val, desc(abs(test_statistic))) %>%
		pull(gene) %>%
		head(k)
	vec2 = sc_da %>%
		arrange(p_val, desc(abs(test_statistic))) %>%
		pull(gene) %>%
		head(k)
	
	concordance_curve = map_dbl(seq_len(k), ~ {
		v1 = vec1[seq_len(.)]
		v2 = vec2[seq_len(.)]
		length(intersect(v1, v2))
	})
	denom = k * (k + 1) / 2
	aucc = sum(concordance_curve) / denom
	out_row = data.frame(c('aucc'=aucc, 'overlap'=length(genes), 'bulk_features'=bulk_size, 'sc_features'=sc_size), 'filter'= filter)
	output_df %<>% rbind(out_row)
	return (output_df)
}

convert_genes_to_go = function(temp_res, go_list, n_permutations=1e4, min_size=10, max_size=1000){
	temp_res %<>%
		mutate(gene = gsub('\\..*', '', gene))

	add_col_grid = temp_res %>%
		dplyr::select(cell_type, de_family, de_method, de_type, paper_da_type) %>%
		distinct()

	gsea = data.frame()

	for (i in 1:nrow(add_col_grid)){
		add_col_row = add_col_grid[i,]
		ranks = temp_res %>%
			filter(
				cell_type == add_col_row$cell_type,
				de_family == add_col_row$de_family,
				de_method == add_col_row$de_method,
				de_type == add_col_row$de_type,
				paper_da_type == add_col_row$paper_da_type
			) %>%
			drop_na(test_statistic) %$%
			setNames(abs(test_statistic), gene) %>% 
			sort(decreasing = TRUE)
		## replace infinite values
		ranks[is.infinite(ranks)] = max(ranks[!is.infinite(ranks)])
		
		temp_gsea = fgsea(pathways = go_list,
					 stats = ranks,
					 nproc = 1,
					 nperm = n_permutations,
					 minSize = min_size,
					 maxSize = max_size) %>%
			dplyr::select(-leadingEdge) %>%
			# flag cell type and comparison 
			dplyr::rename(gene=pathway) %>%
			mutate(
				avg_logFC = 0,
				test_statistic = -log(pval, 10),
				p_val = pval,
				p_val_adj = padj
			) %>%
			cbind(add_col_row) %>%
			data.frame()

		if (nrow(temp_gsea) > 0){
			# temp_gsea = temp_gsea[, colnames(temp_res)]
			temp_gsea %<>%
				select(cell_type, gene, avg_logFC, test_statistic, p_val, p_val_adj, de_family, de_method, de_type, paper_da_type)
			gsea = rbind(gsea, temp_gsea)
		}
	}
	# saveRDS(gsea, save_filename)
	return (gsea)
}