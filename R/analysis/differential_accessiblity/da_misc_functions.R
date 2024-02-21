da_family_args = list(
	'singlecell'=list(
		'wilcox'='wilcox', 
		't'='t',
		'LR'='LR',
		'negbinom'='negbinom',
		'permutation'='permutation',
		'fisher'='exact',
		'binomial'='binomial',
		'LR_peaks'='LR_peaks'
		),
	'pseudobulk'=list(
		'edgeR'=c('LRT', 'QLF'), 
		'DESeq2'=c('LRT','Wald'), 
		'limma'=c('trend','voom')
		),
	'others'=list(
		'snapatac'='exact'
	)
)

preprocess_meta = function(metadata, input_filename) {
	if (grepl('Luecken_2021', input_filename)){
		metadata$label = gsub('\\/', '_', metadata$label)
	}

	if (grepl('Argelaguet_2022', input_filename)){
		metadata %<>% 
			mutate(
				cell_type=gsub('\\.', '_', celltype.predicted),
				replicate=gsub('\\.', '_', sample)
			) %>%
			filter(stage == 'E8.5')
		if (paper_da_type=='cell type'){
			metadata %<>% mutate(
				label=cell_type
			)
		} else if (paper_da_type=='condition'){
			metadata %<>% mutate(
				label=gsub('\\.', '_', stage)
			)
		} else if (paper_da_type=='cell type + condition'){
			metadata %<>% mutate(
				label=cell_type,
				cell_type=gsub('\\.', '_', stage)
			)
		}
	}

	if (grepl('Squair_2022', input_filename)){
		metadata %<>% 
			dplyr::rename(label_ori=label) %>%
			mutate(
				cell_type=gsub(' ', '_', layer2),
				cell_type=gsub(',', '', cell_type),
				replicate=library_id
			)
		if (paper_da_type=='cell type'){
			metadata %<>% mutate(
				label=cell_type
			)
		} else if (paper_da_type=='condition'){
			metadata %<>% mutate(
				label=label_ori
			)
		} else if (paper_da_type=='cell type + condition'){
			metadata %<>% mutate(
				label=cell_type,
				cell_type=label_ori
			)
		}
	}

	if (grepl('Boukhaled_2022', input_filename)){
		metadata %<>% 
			mutate(
				cell_type=cell_annotation
			) %>%
			filter(disease == 'Healthy')
		if (paper_da_type=='cell type'){
			metadata %<>% mutate(
				label=cell_type
			)
		} else if (paper_da_type=='condition'){
			metadata %<>% mutate(
				label=disease
			)
		} else if (paper_da_type=='cell type + condition'){
			metadata %<>% mutate(
				label=cell_type,
				cell_type=disease
			)
		}
	}
}

run_da_wrapper = function(output_df, temp_mat_2, temp_meta_2, args, gc_content_df = NULL, da_types = NULL) {
	if (args$normalization == 'qsmooth'){
		gc_content = gc_content_df[rownames(temp_mat_2),]$gc
		names(gc_content) = rownames(temp_mat_2)

		## manual snippet because the code written is stupidly bugging
		nGroups = 50
		round = TRUE
		gc_groups = Hmisc::cut2(gc_content, g = nGroups)
		gcBinNormCounts = matrix(NA, nrow=nrow(temp_mat_2), ncol=ncol(temp_mat_2),  dimnames=list(rownames(temp_mat_2), colnames(temp_mat_2)))
		for(ii in seq_len(nlevels(gc_groups))){
			id = which(gc_groups==levels(gc_groups)[ii])
			countBin = temp_mat_2[id,]
			qs = qsmooth::qsmooth(SummarizedExperiment(list(counts=countBin)), group_factor=factor(temp_meta_2$replicate))
			normCountBin = qsmooth::qsmoothData(qs)
			if(round) normCountBin = base::round(normCountBin)
			normCountBin[normCountBin<0] = 0
			gcBinNormCounts[id,] = normCountBin
		}
		gcBinNormCounts = as(gcBinNormCounts, "dgCMatrix")
		temp_mat_2 = gcBinNormCounts
		rm(gcBinNormCounts)
		gc()
	}

	tryCatch({
		for (da_type in da_types){
			if (args$experiment=='bulk_ATAC'){
				temp_res = run_de(temp_mat_2, meta = temp_meta_2,
							de_family = args$da_family,
							de_method = args$da_method,
							de_type = da_type,
							min_cells=1,
							n_threads=8)
			} else {
				if (args$normalization == 'f' | args$normalization == 'qsmooth'){
					temp_res = run_de(temp_mat_2, meta = temp_meta_2,
						de_family = args$da_family,
						de_method = args$da_method,
						de_type = da_type,
						n_threads=8)
				} else {
					temp_res = run_de(temp_mat_2, meta = temp_meta_2,
						de_family = args$da_family,
						de_method = args$da_method,
						de_type = da_type,
						n_threads=8,
						normalization=args$normalization)
				}
			}
			if (paper_da_type!='cell type'){
				temp_res %<>% mutate(cell_type=paste0(cell_type, '-', compar1, '_vs_', compar2))
			}
			output_df = rbind(output_df, temp_res)
		},
		error = function(e){
			print(e)
			print(paste0('Skipping ', curr_cell_type))
		},
		finally = {
			message(paste0(curr_cell_type, ' done.'))
		})
	}
	return (output_df)
}