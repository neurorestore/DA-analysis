# Plot a summary figure with performance across all experiments.
setwd("C:/Users/teo/Documents/EPFL/projects/DA-analysis/")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(pals)
library(upstartr)
library(effsize)
source('R/theme.R')

dat = readRDS("data/summaries/biases/matched_bulk_top_k_against_stats.rds") %>% 
    type_convert() %>% 
    mutate(paper = gsub('_scATAC.*', '', basename(filename))) %>%
    dplyr::select(-filename, -seq_type) %>%
    ## no pseudoreplicates
    filter(!pseudo_repl) %>% 
    ## ignore normalization here
    filter(!binarization) %>% 
    ## top-k filter must be none
    filter(top_k_filter == 'None') %>% 
    ## remove ground_truth_* cols
    dplyr::select(-ground_truth_da_method, -ground_truth_da_type) %>% 
    ## distinct
    distinct()
# remove non-unique columns
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

# remove limma and some Gonzalez-Blas
dat %<>%
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    filter(!grepl('limma', da_method)) %>%
    filter(da_method %in% c('LR', 't', 'wilcox'))


# recode normalization
dat %<>% 
    # re-name normalization methods
    mutate(normalization = fct_recode(normalization,
                                      'logTP(10K)' = 'f',
                                      'logTP(median)' = 'log_tp_median',
                                      'smooth GC-full quantile' = 'qsmooth',
                                      'TF-IDF' = 'TFIDF',
                                      'TP(median)' = 'tp_median',
                                      'TP(10K)' = 'tp10k'),
           color = normalization)

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-singlecell') %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           method = fct_recode(method,
                               "'SnapATAC::findDAR'" = 'SnapATAC::findDAR',
                               'Permutation~test' = 'permutation',
                               't~test' = 't',
                               'LR[peaks]' = 'LRpeaks',
                               'LR[clusters]' = 'LR',
                               'Binomial' = 'binomial',
                               'Wilcoxon~rank-sum~test' = 'wilcox',
                               'Fisher~exact~test' = 'FET',
                               'Negative~binomial' = 'negbinom') %>% 
               as.character()) %>% 
    # drop limma
    filter(!grepl('limma', method))
# print 
dplyr::count(dat, color, method) %>% arrange(color, method)

percentile_types = c('mean_expr', 'percentage_open', 'peak_width')
for (percentile_type in percentile_types) {
    if (percentile_type == 'mean_expr') out_plots = list()
    dat$val = dat[,percentile_type]
    
    percentile_type_name = case_when(
        percentile_type == 'mean_expr' ~ 'Read depth',
        percentile_type == 'percentage_open' ~ '% open',
        percentile_type == 'peak_width' ~ 'Peak width (bp)'
    )
    
    # average over bulk methods
    avg1 = dat %>% 
        filter(min_cell_filter == 1, k == 1000, summary_type == 'mean') %>% 
        group_by_at(vars(-starts_with('mean_'), -starts_with('sd_'),
                         -starts_with('percentage_'), -starts_with('-peak_'))) %>% 
        summarise(mean_val = mean(val),
                  n = n()) %>% 
        ungroup()
    
    # plot
    labs = avg1 %>% 
        group_by(k, method, color) %>%
        summarise(mean = mean(val), n = n(),
                  median = median(val), val = median) %>%
        ungroup() %>%
        mutate(label = formatC(median, format = 'f', digits = 2))
    
    grey_pal = colorRampPalette(c('grey10', 'grey85'))(6) %>% 
        setNames(c('none', 'TP(median)', 'TP(10K)', 'logTP(median)', 'smooth GC-full quantile', 'TF-IDF'))
    methods = c('LR[clusters]', 't~test', 'Wilcoxon~rank-sum~test')
    
    curr_out_plots = list()
    for (method in methods) {
        avg2 = avg1 %>% 
            filter(method == !!method)
        labs = avg2 %>% 
            group_by(normalization, method, color) %>%
            summarise(mean = mean(val), n = n(), median = median(val), val = median) %>%
            ungroup() %>%
            mutate(label = formatC(median, format = 'f', digits = 2))
        pal = c('', grey_pal)
        pal[1] = da_analysis_colors[method]
        pal_names = c('logTP(10K)', names(grey_pal))
        names(pal) = pal_names
        
        p1_1 = avg2 %>%
            ggplot(aes(x = reorder_within(normalization, val, method, stats::median), 
                       color = normalization, fill = normalization)) +
            facet_wrap(method ~ ., scales = 'free', ncol = 1, labeller = label_parsed) +
            geom_boxplot(aes(y = val), outlier.shape = NA, alpha = 0.4, width = 0.6,
                         size = 0.4) + 
            coord_flip() +
            scale_y_continuous(percentile_type_name) +
            facet_wrap(~method, scales = 'free_y', labeller = label_parsed) + 
            scale_color_manual('', values = pal) +
            scale_fill_manual('', values = pal) +
            scale_x_reordered() +
            boxed_theme(size_sm = 5, size_lg = 6) +
            theme(axis.title.y = element_blank(),
                  legend.position = 'none',
                  aspect.ratio = 1.2)
        # print(p1_1)
        # if (method != 't~test') {
        #     p1_1 = p1_1 + theme(axis.title.x = element_blank())
        # }
        if (percentile_type == 'mean_expr'){
            p1_1 = p1_1 + 
                geom_text(data = labs, aes(y = -3, label = label), size = 1.75, hjust = 0,
                          color = 'black') +
                scale_y_continuous(percentile_type_name, limits=c(-3, 10))
                
        } else {
            p1_1 = p1_1 + geom_text(data = labs, aes(y = -0.05, label = label), size = 1.75, hjust = 0,
                      color = 'black')
        }
        out_plots[[length(out_plots)+1]] = p1_1
        curr_out_plots[[length(curr_out_plots)+1]] = p1_1
    }
    out_plot_1 = wrap_plots(curr_out_plots, nrow=1)
    ggsave(paste0('fig/Supp_Fig18/bulk-', percentile_type, '-normalization.pdf'), out_plot_1,
           width = 18, height = 6.5, units = "cm", useDingbats = FALSE)
}

full = wrap_plots(out_plots, nrow=3)
# full = plot_grid(out_plots_1[[1]], out_plots_1[[2]], out_plots_1[[3]], nrow=3)
ggsave("fig/Supp_Fig18/full.pdf", full,
       width = 18, height = 15, units = "cm", useDingbats = FALSE)

