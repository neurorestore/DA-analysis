# Plot a summary figure with performance across all experiments.
setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(pals)
library(upstartr)
library(effsize)
source('R/theme.R')

# load data
dat = readRDS("data/summaries/concordance_res/matched_bulk/all_aucc_normalization_only.rds") %>% 
    type_convert()
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

# print some summaries
## bulk DA methods
table(dat$bulk_da_type, dat$bulk_da_method)
## filtering
table(dat$sc_min_cell_filter, dat$sc_percent_filter)
## single-cell methods
dplyr::count(dat, sc_da_method)
## k
table(dat$k)

# average over bulk methods
avg = dat %>% 
    group_by_at(vars(-bulk_da_method, -bulk_da_type,
                     -bulk_features, -overlap, -aucc)) %>% 
    summarise(aucc = mean(aucc),
              n = n()) %>% 
    ungroup()

# is the bulk grid complete?
table(avg$n)
# is the filtering grid complete?
table(avg$sc_min_cell_filter, avg$sc_percent_filter)
# is the comparison grid complete?
dplyr::count(avg, sc_da_method) %>% 
    as.data.frame()

# re-code x-values and colors
avg %<>%
    mutate(method = paste0(sc_da_method, '-singlecell') %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .),
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', 'singlecell') %>% 
               fct_recode('single-cell' = 'singlecell',
                          'other' = 'non_libra') %>% 
               fct_relevel('single-cell', 'pseudobulk', 'other'))
# print 
dplyr::count(avg, color, method) %>% arrange(color, method)

# re-code filtering
avg %<>%
    mutate(filter = ifelse(sc_percent_filter, 
                           paste0(sc_min_cell_filter, '%'),
                           paste0(sc_min_cell_filter, ' cells')) %>% 
               fct_recode('1 cell' = '1 cells') %>% 
               fct_relevel('1 cell', '10 cells', '0.2%', '0.5%', '1%', '5%', '10%', '20%')) %>% 
    filter(!filter %in% c('10%', '20%')) %>% 
    mutate(method = fct_recode(method,
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

###############################################################################-
## AUCC: boxplot, main text ####
###############################################################################-

avg1 = avg %>% 
    filter(!grepl('pseudo-replicates', method)) %>% 
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    filter(k == 5000, filter == '1 cell') %>% 
    # re-name normalization methods
    mutate(normalization = fct_recode(sc_normalization,
                                      'logTP(10K)' = 'f',
                                      'logTP(median)' = 'log_tp_median',
                                      'smooth GC-full quantile' = 'qsmooth',
                                      'TF-IDF' = 'TFIDF',
                                      'TP(median)' = 'tp_median',
                                      'TP(10K)' = 'tp10k'))
dplyr::count(avg1, method, normalization)

methods = c('LR[clusters]', 't~test', 'Wilcoxon~rank-sum~test')
out_plots = list()

grey_pal = colorRampPalette(c('grey10', 'grey85'))(6) %>% 
    setNames(c('none', 'TP(median)', 'TP(10K)', 'logTP(median)', 'TF-IDF', 'smooth GC-full quantile'))

curr_out_plots = list()
for (method in methods) {
    avg2 = avg1 %>% 
        filter(method == !!method)
    labs = avg2 %>% 
        group_by(k, normalization, filter, method, color) %>%
        summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
        ungroup() %>%
        mutate(label = formatC(median, format = 'f', digits = 2))
    pal = c('', grey_pal)
    pal[1] = da_analysis_colors[method]
    pal_names = c('logTP(10K)', names(grey_pal))
    names(pal) = pal_names
    
    p1_1 = avg2 %>%
        ggplot(aes(x = reorder_within(normalization, aucc, method, stats::median), 
                   color = normalization, fill = normalization)) +
        facet_wrap(method ~ ., scales = 'free', ncol = 1, labeller = label_parsed) +
        geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                     size = 0.4) + 
        geom_text(data = labs, aes(y = -0.12, label = label), size = 1.75, hjust = 0,
                  color = 'black') +
        scale_color_manual('', values = pal) +
        scale_fill_manual('', values = pal) +
        scale_x_reordered() +
        scale_y_continuous('AUCC', breaks = seq(0, 1, 0.1)) +
        coord_flip() +
        boxed_theme(size_sm = 5, size_lg = 6) +
        theme(axis.title.y = element_blank(),
              legend.position = 'none',
              legend.key.width = unit(0.4, 'lines'),
              aspect.ratio = 1.2)
    # p1
    out_plots[[length(out_plots)+1]] = p1_1
    curr_out_plots[[length(curr_out_plots)+1]] = p1_1
}

p1 = wrap_plots(curr_out_plots, nrow=1)
ggsave("fig/Supp_Fig17/bulk-aucc-normalization.pdf", p1,
       width = 18, height = 6.5, units = "cm", useDingbats = FALSE)

# load data
dat = readRDS("data/summaries/false_discoveries/false_discoveries_luecken.rds") %>% 
    type_convert() %>% 
    # flag study
    mutate(study = gsub("-.*$", "", filename)) %>% 
    dplyr::select(-filename) %>% 
    # filter out methods we are not using
    filter(pseudo_repl == FALSE) %>% 
    # ignore percentile here
    filter(is.na(percentile), !binarization, is.na(latent_vars)) %>%
    filter(da_method %in% c('LR', 't', 'wilcox'))
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

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
               as.character())

# recode normalization
dat %<>% 
    # re-name normalization methods
    mutate(normalization = fct_recode(normalization,
                                      'logTP(10K)' = 'f',
                                      'logTP(median)' = 'log_tp_median',
                                      'smooth GC-full quantile' = 'qsmooth',
                                      'TF-IDF' = 'TFIDF',
                                      'TP(median)' = 'tp_median',
                                      'TP(10K)' = 'tp10k'))
curr_out_plots = list()
for (method in methods) {
    dat1 = dat %>% 
        filter(method == !!method)
    labs = dat1 %>% 
        group_by(method, normalization) %>%
        summarise(mean = mean(number_of_da_regions), 
                  n = n(), 
                  median = median(number_of_da_regions), 
                  number_of_da_regions = median) %>%
        ungroup() %>%
        mutate(label = round(median))
    pal = c('', grey_pal)
    pal[1] = da_analysis_colors[method]
    pal_names = c('logTP(10K)', names(grey_pal))
    names(pal) = pal_names
    
    p2_1 = dat1 %>%
        ggplot(aes(x = reorder_within(normalization, number_of_da_regions, method, stats::median), 
                   color = normalization, fill = normalization)) +
        facet_wrap(method ~ ., scales = 'free', ncol = 1, labeller = label_parsed) +
        geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, 
                     alpha = 0.4, width = 0.6, size = 0.4) + 
        geom_text(data = labs, aes(y = -2e3, label = label), size = 1.75, hjust = 1,
                  color = 'black') +
        scale_color_manual('', values = pal) +
        scale_fill_manual('', values = pal) +
        scale_x_reordered() +
        scale_y_continuous(expression('# of DA peaks'~(10^3)), 
                           labels = function(x) x / 1e3,
                           breaks = seq(0, 80, 10) * 1e3) +
        coord_flip(ylim = c(-12e3, 55e3)) +
        guides(color = guide_legend(ncol = 1)) +
        boxed_theme(size_sm = 5, size_lg = 6) +
        theme(axis.title.y = element_blank(),
              aspect.ratio = 1.2,
              legend.position = 'none',
              legend.key.width = unit(0.35, 'lines'),
              legend.key.height = unit(0.5, 'lines'))
    p2_1
    out_plots[[length(out_plots)+1]] = p2_1
    curr_out_plots[[length(curr_out_plots)+1]] = p2_1
}

p2 = wrap_plots(curr_out_plots, nrow=1)
ggsave("fig/Supp_Fig17/luecken-false-discoveries-normalization.pdf", p2,
       width = 18, height = 6.5, units = "cm", useDingbats = FALSE)

# load data
dat = readRDS("data/summaries/false_discoveries/false_discoveries_snapatac.rds") %>% 
    type_convert() %>% 
    dplyr::select(-ori_filename) %>%
    # filter out methods we are not using
    filter(da_family != 'mixedmodel',
           pseudo_repl == FALSE) %>% 
    # entire range of expression values only - no binarization
    filter(!binarization) %>% 
    # normalization: only LR, t, wilcox
    filter(da_method %in% c('LR', 't', 'wilcox'), is.na(percentile))
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-singlecell') %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', 'singlecell') %>% 
               fct_recode('single-cell' = 'singlecell',
                          'other' = 'non_libra') %>% 
               fct_relevel('single-cell', 'pseudobulk', 'other'),
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

# recode normalization
dat %<>% 
    # re-name normalization methods
    mutate(normalization = fct_recode(normalization,
                                      'logTP(10K)' = 'f',
                                      'logTP(median)' = 'log_tp_median',
                                      'smooth GC-full quantile' = 'qsmooth',
                                      'TF-IDF' = 'TFIDF',
                                      'TP(median)' = 'tp_median',
                                      'TP(10K)' = 'tp10k'))

curr_out_plots = list()
for (method in methods) {
    dat1 = dat %>% 
        filter(method == !!method)
    labs = dat1 %>% 
        group_by(method, normalization) %>%
        summarise(mean = mean(number_of_da_regions), 
                  n = n(), 
                  median = median(number_of_da_regions), 
                  number_of_da_regions = median) %>%
        ungroup() %>%
        mutate(label = round(median))
    pal = c('', grey_pal)
    pal[1] = da_analysis_colors[method]
    pal_names = c('logTP(10K)', names(grey_pal))
    names(pal) = pal_names
    
    p3_1 = dat1 %>%
        ggplot(aes(x = reorder_within(normalization, number_of_da_regions, method, stats::median), 
                   color = normalization, fill = normalization)) +
        facet_wrap(method ~ ., scales = 'free', ncol = 1, labeller = label_parsed) +
        geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, 
                     alpha = 0.4, width = 0.6, size = 0.4) + 
        geom_text(data = labs, aes(y = -1e3, label = label), size = 1.75, hjust = 1,
                  color = 'black') +
        scale_color_manual('', values = pal) +
        scale_fill_manual('', values = pal) +
        scale_x_reordered() +
        scale_y_continuous(expression('# of DA peaks'~(10^3)), 
                           labels = function(x) x / 1e3,
                           breaks = seq(0, 80, 10) * 1e3) +
        coord_flip(ylim = c(-13e3, 55e3)) +
        guides(color = guide_legend(ncol = 1)) +
        boxed_theme(size_sm = 5, size_lg = 6) +
        theme(axis.title.y = element_blank(),
              aspect.ratio = 1.2,
              legend.position = 'none',
              legend.key.width = unit(0.35, 'lines'),
              legend.key.height = unit(0.5, 'lines'))
    p3_1
    out_plots[[length(out_plots)+1]] = p3_1
    curr_out_plots[[length(curr_out_plots)+1]] = p3_1
}

p3 = wrap_plots(curr_out_plots, nrow=1)
ggsave("fig/Supp_Fig17/snapatac-false-discoveries-normalization.pdf", p3,
       width = 18, height = 6.5, units = "cm", useDingbats = FALSE)


# load data
dat = readRDS("data/summaries/false_discoveries/false_discoveries_simulated.rds") %>% 
    type_convert() %>% 
    dplyr::select(-ori_filename) %>%
    # filter out methods we are not using
    filter(da_family != 'mixedmodel',
           pseudo_repl == FALSE,
           is.na(percentile)) %>% 
    # entire range of expression values only - no binarization
    filter(!binarization) %>% 
    # normalization: only LR, t, wilcox
    filter(da_method %in% c('LR', 't', 'wilcox'))
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)
dat %<>%  filter(de_loc == 1, de_prob == 0.5, n_reps == 6, n_cells == 1000)
# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-singlecell') %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', 'singlecell') %>% 
               fct_recode('single-cell' = 'singlecell',
                          'other' = 'non_libra') %>% 
               fct_relevel('single-cell', 'pseudobulk', 'other'),
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

# recode normalization
dat %<>% 
    # re-name normalization methods
    mutate(normalization = fct_recode(normalization,
                                      'logTP(10K)' = 'f',
                                      'logTP(median)' = 'log_tp_median',
                                      'smooth GC-full quantile' = 'qsmooth',
                                      'TF-IDF' = 'TFIDF',
                                      'TP(median)' = 'tp_median',
                                      'TP(10K)' = 'tp10k'))
curr_out_plots = list()
for (method in methods) {
    dat1 = dat %>% 
        filter(method == !!method)
    labs = dat1 %>% 
        group_by(method, normalization) %>%
        summarise(mean = mean(number_of_da_regions), 
                  n = n(), 
                  median = median(number_of_da_regions), 
                  number_of_da_regions = median) %>%
        ungroup() %>%
        mutate(label = round(median))
    pal = c('', grey_pal)
    pal[1] = da_analysis_colors[method]
    pal_names = c('logTP(10K)', names(grey_pal))
    names(pal) = pal_names
    
    p4_1 = dat1 %>%
        ggplot(aes(x = reorder_within(normalization, number_of_da_regions, method, stats::median), 
                   color = normalization, fill = normalization)) +
        facet_wrap(method ~ ., scales = 'free', ncol = 1, labeller = label_parsed) +
        geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, 
                     alpha = 0.4, width = 0.6, size = 0.4) + 
        geom_text(data = labs, aes(y = -1e2, label = label), size = 1.75, hjust = 1,
                  color = 'black') +
        scale_color_manual('', values = pal) +
        scale_fill_manual('', values = pal) +
        scale_x_reordered() +
        scale_y_continuous(expression('# of DA peaks'~(10^3)), 
                           labels = function(x) x / 1e3,
                           breaks = seq(0, 5, 1) * 1e3) +
        coord_flip(ylim = c(-1e3, 5e3)) +
        guides(color = guide_legend(ncol = 1)) +
        boxed_theme(size_sm = 5, size_lg = 6) +
        theme(axis.title.y = element_blank(),
              aspect.ratio = 1.2,
              legend.position = 'none',
              legend.key.width = unit(0.35, 'lines'),
              legend.key.height = unit(0.5, 'lines'))
    p4_1
    out_plots[[length(out_plots)+1]] = p4_1
    curr_out_plots[[length(curr_out_plots)+1]] = p4_1
}
p4 = wrap_plots(curr_out_plots, nrow=1)
ggsave("fig/Supp_Fig17/splatter-false-discoveries-normalization.pdf", p4,
       width = 18, height = 6.5, units = "cm", useDingbats = FALSE)

# full = wrap_plots(p1, p2, p3, p4, nrow=4)
full = wrap_plots(out_plots, nrow=4)
ggsave("fig/Supp_Fig17/full.pdf", full,
       width = 18, height = 20, units = "cm", useDingbats = FALSE)
