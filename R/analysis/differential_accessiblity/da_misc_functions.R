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
		}
	},
	error = function(e){
		print(e)
		print(paste0('Skipping ', curr_cell_type))
	},
	finally = {
		message(paste0(curr_cell_type, ' done.'))
	})
	return (output_df)
}

## alternative implementations from ArchR and glmGamPoi
sparseMatWilcoxon = function(mat1, mat2, meta1 = NULL, meta2 = NULL){
	n1 = ncol(mat1)
    n2 = ncol(mat2)
    ## no background matching ##
    # stopifnot(n1==n2)

	df = wilcoxauc(cbind(mat1,mat2), c(rep("Top", ncol(mat1)),rep("Bot", ncol(mat2))))
	df = df[which(df$group=="Top"),]

	#Sparse Row Sums
	m1 = Matrix::rowSums(mat1, na.rm=TRUE)
	m2 = Matrix::rowSums(mat2, na.rm=TRUE)
	offset = 1 #quantile(c(mat1@x,mat2@x), 0.99) * 10^-4
	log2FC = log2((m1 + offset) / (m2 + offset))
	log2Mean = log2(((m1 + offset) + (m2 + offset)) / 2)

	out = data.frame(
		log2Mean = log2Mean,
		log2FC = log2FC,
		fdr = df$padj, 
		pval = df$pval, 
		mean1 = Matrix::rowMeans(mat1, na.rm=TRUE), 
		mean2 = Matrix::rowMeans(mat2, na.rm=TRUE), 
		n = ncol(mat1),
		auc = df$auc
	)

	return(out)

}

sparseMatTTest = function(mat1, mat2, meta1 = NULL, meta2 = NULL, m0 = 0){
    
    n1 = ncol(mat1)
    n2 = ncol(mat2)
    ## no background matching ##
    # stopifnot(n1==n2)
    n = n1 + n2
    
    #Sparse Row Means
    m1 = Matrix::rowMeans(mat1, na.rm=TRUE)
    m2 = Matrix::rowMeans(mat2, na.rm=TRUE)
    
    #Sparse Row Variances
    v1 = computeSparseRowVariances(mat1@i + 1, mat1@x, m1, n1)
    v2 = computeSparseRowVariances(mat2@i + 1, mat2@x, m2, n2)
    
    #Calculate T Statistic
    se = sqrt( (1/n1 + 1/n2) * ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2) )
    tstat = (m1-m2-m0)/se
    pvalue = 2*pt(-abs(tstat), n - 2)
    fdr = p.adjust(pvalue, method = "fdr")
    
    #Sparse Row Sums
    m1 = Matrix::rowSums(mat1, na.rm=TRUE)
    m2 = Matrix::rowSums(mat2, na.rm=TRUE)
    offset = 1 #quantile(c(mat1@x,mat2@x), 0.99) * 10^-4
    log2FC = log2((m1 + offset)/(m2 + offset))
    log2Mean = log2(((m1+offset) + (m2+offset)) / 2)

    out = data.frame(
      log2Mean = log2Mean,
      log2FC = log2FC,
      fdr = fdr, 
      pval = pvalue, 
      mean1 = m1 / n1, 
      mean2 = m2 / n2, 
      var1 = v2,
      var2 = v2,
      n = n1
    )
    return(out)
}

glmGamPoi_de = function(mat1, mat2, meta1 = NULL, meta2 = NULL) {
	full_mat = as.matrix(cbind(mat1, mat2))
	full_meta = rbind(meta1, meta2)
	full_meta %<>% 
		mutate(
			label = gsub(' ', '_', label),
			label = gsub('-', '_', label),
			label = gsub('\\+', '_', label)
		)

	fit = glm_gp(full_mat, col_data = full_meta, design = ~ label, on_disk = FALSE)
	de_res = test_de(fit, contrast = colnames(fit$Beta)[2]) 
	return(de_res)
}

run_alt_DA = function(mat1, mat2, meta1=NULL, meta2=NULL, test) {
	if (test == 'ArchR_wilcoxon') {
		out = sparseMatWilcoxon(mat1=mat1, mat2=mat2, meta1=meta1, meta2=meta2)
	} else if (test == 'ArchR_t') {
		out = sparseMatTTest(mat1=mat1, mat2=mat2, meta1=meta1, meta2=meta2))
	} else if (test == 'glmGamPoi') {
		out = glmGamPoi_de(mat1=mat1, mat2=mat2, meta1=meta1, meta2=meta2)
	}
	out %<>% mutate(
		da_family = 'singlecell',
		da_method = test,
		da_type = test
	)
}