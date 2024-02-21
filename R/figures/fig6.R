setwd("C:/Users/teo/Documents/EPFL/projects/DA-analysis/")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(effsize)
source("R/theme.R")

###############################################################################-
## bulk binarization plot ####
###############################################################################-
# load data
dat = readRDS("data/summaries/concordance_res/matched_bulk/all_aucc_full_best_practice_check.rds") %>% 
    type_convert() %>% 
    # ignore percentiles here
    filter(expr == 'all') 

# remove non-unique columns
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)
# remove limma and some Gonzalez-Blas
dat %<>%
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    filter(!grepl('limma', sc_da_method))

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(sc_da_method, '-', sc_da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(sc_binarization, 'binarized', sc_da_family),
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', color) %>% 
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

# average over bulk methods
avg1 = dat %>% 
    filter(sc_min_cell_filter == 1, sc_percent_filter == FALSE,
           k == 5000) %>%
    group_by_at(vars(-bulk_da_method, -bulk_da_type,
                     -overlap, -aucc, -bulk_features, -sc_features)) %>% 
    summarise(aucc = mean(aucc),
              n = n()) %>% 
    ungroup()
# plot
labs = avg1 %>% 
    group_by(k, method, color) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
pal = c(da_analysis_colors, 'binarized' = "#6E645F")

p1 = avg1 %>%
    mutate(
        color = ifelse(color != 'binarized', method, color)
    ) %>%
    ggplot(aes(x = reorder(method, aucc, stats::median), 
               color = color, fill = color)) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.65,
                 size = 0.35) + 
    geom_text(data = labs, aes(y = -0.08, group = color, label = label),
              position = position_dodge(width = 0.65),
              size = 1.75, hjust = 0, show.legend = FALSE) +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('AUCC', breaks = seq(0, 1, 0.1)) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          # legend.position = 'none',
          legend.key.width = unit(0.18, 'lines'),
          legend.key.height = unit(0.14, 'lines'),
          aspect.ratio = 1.4)
p1
ggsave("fig/Fig6/bulk-aucc-binarization.pdf", p1,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)

p1 = p1 + theme(legend.position = 'none')

###############################################################################-
## luecken binarization plot ####
###############################################################################-

# load data
dat = readRDS("data/summaries/false_discoveries/false_discoveries_luecken.rds") %>% 
    type_convert() %>% 
    # flag study
    mutate(study = gsub("-.*$", "", filename)) %>% 
    dplyr::select(-filename) %>% 
    # filter out methods we are not using
    filter(da_family != 'mixedmodel',
           pseudo_repl == FALSE) %>% 
    # ignore percentile here
    filter(is.na(percentile), normalization == 'f', is.na(latent_vars))
# remove non-unique columns
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)
# remove limma and some Gonzalez-Blas
dat %<>%
    filter(!grepl('limma', da_method))

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-', da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(binarization, 'binarized', da_family),
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', color) %>% 
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

# filter to only methods that _can_ be binarized
dat %<>%
    group_by(method) %>% 
    filter(n_distinct(color) == 2) %>%
    ungroup()
dplyr::count(dat, color, method) %>% arrange(method, color)

labs = dat %>% 
    group_by(method, color) %>%
    summarise(mean = mean(number_of_da_regions), 
              n = n(), 
              median = median(number_of_da_regions), 
              number_of_da_regions = median) %>%
    ungroup() %>%
    mutate(label = round(median))
pal = c(da_analysis_colors, 'binarized' = "#6E645F")
p2 = dat %>%
    mutate(
        color = ifelse(color != 'binarized', method, color)
    ) %>%
    ggplot(aes(x = reorder(method, number_of_da_regions, stats::median), 
               color = color, fill = color)) +
    geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, 
                 alpha = 0.4, width = 0.6, size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -18e3, group = color, label = label),
              position = position_dodge(width = 0.65),
              size = 1.75, hjust = 0, show.legend = FALSE) +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous(expression('# of DA peaks'~(10^3)), 
                       labels = function(x) x / 1e3,
                       breaks = seq(0, 80, 10) * 1e3) +
    coord_flip(ylim = c(-18e3, 55e3)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none'
          )
p2
ggsave("fig/Fig6/luecken-false-discoveries-binarization.pdf", p2,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)


###############################################################################-
## snapatac binarization plot ####
###############################################################################-

# load data
dat = readRDS("data/summaries/false_discoveries/false_discoveries_snapatac.rds") %>% 
    type_convert() %>% 
    dplyr::select(-ori_filename, -data_filename) %>%
    # ignore pseudoreplicates
    filter(!pseudo_repl) %>% 
    # no normalization
    filter(normalization == 'f') %>% 
    # ignore percentile here
    filter(is.na(percentile))
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

# remove limma
dat %<>%
    filter(!grepl('limma', da_method))

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(de_method, '-', de_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(binarization, 'binarized', de_family),
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', color) %>% 
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

# filter to only methods that _can_ be binarized
dat %<>%
    group_by(method) %>% 
    filter(n_distinct(color) == 2) %>%
    ungroup()

labs = dat %>% 
    group_by(method, color) %>%
    summarise(mean = mean(number_of_da_regions), 
              n = n(), 
              median = median(number_of_da_regions), 
              number_of_da_regions = median) %>%
    ungroup() %>%
    mutate(label = round(median))
pal = c(da_analysis_colors, 'binarized' = "#6E645F")
p3 = dat %>%
    mutate(
        color = ifelse(color != 'binarized', method, color)
    ) %>%
    ggplot(aes(x = reorder(method, number_of_da_regions, stats::median), 
               color = color, fill = color)) +
    geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, 
                 alpha = 0.4, width = 0.6, size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -20e3, group = color, label = label),
              position = position_dodge(width = 0.65),
              size = 1.75, hjust = 0, show.legend = FALSE) +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous(expression('# of DA peaks'~(10^3)), 
                       labels = function(x) x / 1e3,
                       breaks = seq(0, 80, 10) * 1e3) +
    coord_flip(ylim = c(-20e3, 62e3)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none'
          )
p3
ggsave("fig/Fig6/snapatac-false-discoveries-binarization.pdf", p3,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)


###############################################################################-
## simulated binarization plot ####
###############################################################################-

# load data
dat = readRDS("data/summaries/false_discoveries/false_discoveries_simulated.rds") %>% 
    type_convert() %>% 
    dplyr::select(-ori_filename) %>%
    # ignore pseudoreplicates
    filter(!pseudo_repl) %>% 
    # ignore percentile here
    filter(is.na(percentile), normalization == 'f')
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep) %>%
    distinct()
dat %<>%  filter(de_loc == 1, de_prob == 0.5, n_reps == 6, n_cells == 1000)
# remove limma
dat %<>%
    filter(!grepl('limma', da_method))

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(de_method, '-', de_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(binarization, 'binarized', de_family),
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', color) %>% 
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

# filter to only methods that _can_ be binarized
dat %<>%
    group_by(method) %>% 
    filter(n_distinct(color) == 2) %>%
    ungroup()

labs = dat %>% 
    group_by(method, color) %>%
    summarise(mean = mean(number_of_da_regions), 
              n = n(), 
              median = median(number_of_da_regions), 
              number_of_da_regions = median) %>%
    ungroup() %>%
    mutate(label = round(median))
pal = c(da_analysis_colors, 'binarized' = "#6E645F")
p4 = dat %>%
    mutate(
        color = ifelse(color != 'binarized', method, color)
    ) %>%
    ggplot(aes(x = reorder(method, number_of_da_regions, stats::median), 
               color = color, fill = color)) +
    geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, 
                 alpha = 0.4, width = 0.6, size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -10e3, group = color, label = label),
              position = position_dodge(width = 0.65),
              size = 1.75, hjust = 0, show.legend = FALSE) +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous(expression('# of DA peaks'~(10^3)), 
                       labels = function(x) x / 1e3,
                       breaks = seq(0, 30, 10) * 1e3) +
    coord_flip(ylim = c(-10e3, 31e3)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none'
          )
p4
ggsave("fig/Fig6/splatter-false-discoveries-binarization.pdf", p4,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)

###############################################################################-
## bulk binarization biases ####
###############################################################################-

# load data
dat = readRDS("data/summaries/biases/matched_bulk_top_k_against_stats.rds") %>% 
    type_convert() %>% 
    mutate(
        paper = gsub('_scATAC.*', '', basename(filename)) 
    ) %>%
    dplyr::select(-filename, -seq_type) %>% 
    ## no pseudoreplicates
    filter(!pseudo_repl) %>% 
    ## ignore normalization here
    filter(normalization == 'f') %>% 
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
    filter(!grepl('limma', da_method))

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-', da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(binarization, 'binarized', da_family),
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', color) %>% 
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

# filter to only methods that _can_ be binarized
dat %<>%
    group_by(method) %>% 
    filter(n_distinct(color) == 2) %>%
    ungroup()

# average over bulk methods
avg1 = dat %>% 
    filter(min_cell_filter == 1, k == 1000, summary_type == 'mean') %>% 
    group_by_at(vars(-starts_with('mean_'), -starts_with('sd_'),
                     -starts_with('percentage_'), -starts_with('-peak_'))) %>% 
    summarise(mean_expr = mean(mean_expr),
              n = n()) %>% 
    ungroup()

# plot
labs = avg1 %>% 
    group_by(k, method, color) %>%
    summarise(mean = mean(mean_expr), n = n(),
              median = median(mean_expr), mean_expr = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
pal = c(da_analysis_colors, 'binarized' = "#6E645F")
p5 = avg1 %>%
    mutate(
        color = ifelse(color != 'binarized', method, color)
    ) %>%
    ggplot(aes(x = reorder(method, mean_expr, stats::median), 
               color = color, fill = color)) +
    geom_boxplot(aes(y = mean_expr), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    geom_text(data = labs, aes(y = -1.5, group = color, label = label),
              position = position_dodge(width = 0.65),
              size = 1.75, hjust = 0, show.legend = FALSE) +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('Read depth') +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none'
          )
p5
ggsave("fig/Fig6/bulk-mean-expr-binarization.pdf", p5,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)

# average over bulk methods
avg3 = dat %>% 
    filter(min_cell_filter == 1, k == 1000, summary_type == 'mean') %>% 
    group_by_at(vars(-starts_with('mean_'), -starts_with('sd_'),
                     -starts_with('percentage_'), -starts_with('-peak_'))) %>% 
    summarise(percentage_open = mean(percentage_open),
              n = n()) %>% 
    ungroup()

# plot
labs = avg3 %>% 
    group_by(k, method, color) %>%
    summarise(mean = mean(percentage_open), n = n(),
              median = median(percentage_open), percentage_open = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
p6 = avg3 %>%
    mutate(
        color = ifelse(color != 'binarized', method, color)
    ) %>%
    ggplot(aes(x = reorder(method, percentage_open, stats::median), 
               color = color, fill = color)) +
    geom_boxplot(aes(y = percentage_open), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    geom_text(data = labs, aes(y = -0.05, group = color, label = label),
              position = position_dodge(width = 0.65),
              size = 1.75, hjust = 0, show.legend = FALSE) +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('% open', labels = ~ . * 100) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none')
p6
ggsave("fig/Fig6/bulk-pct-open-binarization.pdf", p6,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)

# average over bulk methods
avg5 = dat %>% 
    filter(min_cell_filter == 1, k == 1000, summary_type == 'mean') %>% 
    group_by_at(vars(-starts_with('mean_'), -starts_with('sd_'),
                     -starts_with('percentage_'), -starts_with('-peak_'))) %>% 
    summarise(peak_width = mean(peak_width),
              n = n()) %>% 
    ungroup()

# plot
labs = avg5 %>% 
    group_by(k, method, color) %>%
    summarise(mean = mean(peak_width), n = n(),
              median = median(peak_width), peak_width = median) %>%
    ungroup() %>%
    mutate(label = round(median))
pal = c(da_analysis_colors, 'binarized' = "#6E645F")
p7 = avg5 %>%
    mutate(
        color = ifelse(color != 'binarized', method, color)
    ) %>%
    ggplot(aes(x = reorder(method, peak_width, stats::median), 
               color = color, fill = color)) +
    geom_boxplot(aes(y = peak_width), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    geom_text(data = labs, aes(y = -0.05, group = color, label = label),
              position = position_dodge(width = 0.65),
              size = 1.75, hjust = 0, show.legend = FALSE) +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('Peak width (bp)') +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          aspect.ratio = 1.4)
p7
ggsave("fig/Fig6/bulk-peak-width-binarization.pdf", p7,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)


###############################################################################-
## bulk normalization effect ####
###############################################################################-

# load data
dat = readRDS("data/summaries/concordance_res/matched_bulk/all_aucc_normalization_only.rds") %>% 
    type_convert()
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

# average over bulk methods
avg = dat %>% 
    group_by_at(vars(-bulk_da_method, -bulk_da_type,
                     -bulk_features, -overlap, -aucc)) %>% 
    summarise(aucc = mean(aucc),
              n = n()) %>% 
    ungroup()

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
delta = avg1 %>% 
    group_by_at(vars(-sc_normalization, -normalization, -sc_features, -aucc)) %>% 
    mutate(delta = aucc - aucc[normalization == 'logTP(10K)']) %>% 
    ungroup()

grid = distinct(delta, method, normalization) %>%
    filter(normalization != 'logTP(10K)')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    norm = as.character(current$normalization)
    df = filter(avg1, method == current$method, 
                normalization %in% c(norm, 'logTP(10K)')) %>% 
        droplevels() %>% 
        mutate(normalization = factor(normalization, levels = c(norm, 'logTP(10K)')))
    eff = cohen.d(formula = aucc ~ normalization | Subject(cell_type),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})
color_pal = da_analysis_colors[c("LR[clusters]",  "t~test", "Wilcoxon~rank-sum~test")]

p8 = eff %>%
    ggplot(aes(x = reorder(normalization, -d, mean),
               color = method, fill = method, y = d)) +
    geom_hline(aes(yintercept = 0), linetype = 'solid', size = 0.4,
               color = 'grey82') +
    geom_point(position = position_dodge(width = 0.6), shape = 21, size = 0.8) + 
    geom_errorbar(aes(ymin = 0, ymax = d, group = method), width = 0, size = 0.3,
                  position = position_dodge(width = 0.6)) +
    scale_y_continuous('Cohen\'s d,\nrelative to log(TP10K)') +
    ggtitle('Cohen\'s d, AUCC') +
    scale_color_manual(values = color_pal,labels = parse_format()) +
    scale_fill_manual(values = color_pal,labels = parse_format()) +
    boxed_theme() +
    theme(aspect.ratio = 1, 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size=6),
          legend.position = 'bottom',
          legend.text.align = 0,
          legend.title = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), 'cm')
          ) +
    guides(color=guide_legend(ncol=1))
p8
ggsave("fig/Fig6/bulk-aucc-normalization-d.pdf", p8, width = 3.8, height = 4.6, units = "cm", useDingbats = FALSE)

###############################################################################-
## luecken normalization effect ####
###############################################################################-

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

delta = dat %>% 
    group_by(da_method, label) %>% 
    mutate(delta = number_of_da_regions - 
               number_of_da_regions[normalization == 'logTP(10K)'],
           n = n()) %>% 
    ungroup()

grid = distinct(delta, method, normalization) %>%
    filter(normalization != 'logTP(10K)')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    norm = as.character(current$normalization)
    df = filter(dat, method == current$method, 
                normalization %in% c(norm, 'logTP(10K)')) %>% 
        droplevels() %>% 
        mutate(normalization = factor(normalization, levels = c(norm, 'logTP(10K)')))
    eff = cohen.d(formula = number_of_da_regions ~ normalization | Subject(label),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})

p9 = eff %>%
    ggplot(aes(x = reorder(normalization, -d, mean),
               color = method, fill = method, y = d)) +
    geom_hline(aes(yintercept = 0), linetype = 'solid', size = 0.4,
               color = 'grey82') +
    geom_point(position = position_dodge(width = 0.6), shape = 21, size = 0.8) + 
    geom_errorbar(aes(ymin = 0, ymax = d, group = method), width = 0, size = 0.3,
                  position = position_dodge(width = 0.6)) +
    scale_y_continuous('Cohen\'s d,\nrelative to log(TP10K)') +
    ggtitle('Cohen\'s d, false discoveries') +
    scale_color_manual(values = color_pal,labels = parse_format()) +
    scale_fill_manual(values = color_pal,labels = parse_format()) +
    boxed_theme() +
    theme(aspect.ratio = 1, 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size=6),
          legend.position = 'bottom',
          legend.text.align = 0,
          legend.title = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), 'cm')
    ) +
    guides(color=guide_legend(ncol=1))
p9
ggsave("fig/Fig6/luecken-false-discoveries-normalization-d.pdf", p9, width = 3.8, height = 4.6, units = "cm", useDingbats = FALSE)

###############################################################################-
## snapatac normalization effect ####
###############################################################################-
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

delta = dat %>% 
    group_by(da_method, seed) %>% 
    mutate(delta = number_of_da_regions - 
               number_of_da_regions[normalization == 'logTP(10K)'],
           n = n()) %>% 
    ungroup()
grid = distinct(delta, method, normalization) %>%
    filter(normalization != 'logTP(10K)')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    norm = as.character(current$normalization)
    df = filter(dat, method == current$method, 
                normalization %in% c(norm, 'logTP(10K)')) %>% 
        droplevels() %>% 
        mutate(normalization = factor(normalization, levels = c(norm, 'logTP(10K)')))
    eff = cohen.d(formula = number_of_da_regions ~ normalization | Subject(seed),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
}) %>%
    mutate(d = ifelse(is.nan(d), 0, d))

p10 = eff %>%
    ggplot(aes(x = reorder(normalization, -d, mean, na.rm = TRUE),
               color = method, fill = method, y = d)) +
    geom_hline(aes(yintercept = 0), linetype = 'solid', size = 0.4,
               color = 'grey82') +
    geom_point(position = position_dodge(width = 0.6), shape = 21, size = 0.8) + 
    geom_errorbar(aes(ymin = 0, ymax = d, group = method), width = 0, size = 0.3,
                  position = position_dodge(width = 0.6)) +
    scale_y_continuous('Cohen\'s d,\nrelative to log(TP10K)') +
    ggtitle('Cohen\'s d, false discoveries') +
    scale_color_manual(values = color_pal,labels = parse_format()) +
    scale_fill_manual(values = color_pal,labels = parse_format()) +
    boxed_theme() +
    theme(aspect.ratio = 1, 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size=6),
          legend.position = 'bottom',
          legend.text.align = 0,
          legend.title = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), 'cm')
    ) +
    guides(color=guide_legend(ncol=1))
p10
ggsave("fig/Fig6/snapatac-false-discoveries-normalization-d.pdf", p10, width = 3.8, height = 4.6, units = "cm", useDingbats = FALSE)
    
###############################################################################-
## splatter normalization effect ####
###############################################################################-

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

delta = dat %>% 
    group_by(da_method, simulation_iter) %>% 
    mutate(delta = number_of_da_regions - 
               number_of_da_regions[normalization == 'logTP(10K)'],
           n = n()) %>% 
    ungroup()
grid = distinct(delta, method, normalization) %>%
    filter(normalization != 'logTP(10K)')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    norm = as.character(current$normalization)
    df = filter(dat, method == current$method, 
                normalization %in% c(norm, 'logTP(10K)')) %>% 
        droplevels() %>% 
        mutate(normalization = factor(normalization, levels = c(norm, 'logTP(10K)')))
    eff = cohen.d(formula = number_of_da_regions ~ normalization | Subject(simulation_iter),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
}) %>%
    mutate(d = ifelse(is.nan(d), 0, d))

p11 = eff %>%
    ggplot(aes(x = reorder(normalization, -d, mean, na.rm = TRUE),
               color = method, fill = method, y = d)) +
    geom_hline(aes(yintercept = 0), linetype = 'solid', size = 0.4,
               color = 'grey82') +
    geom_point(position = position_dodge(width = 0.6), shape = 21, size = 0.8) + 
    geom_errorbar(aes(ymin = 0, ymax = d, group = method), width = 0, size = 0.3,
                  position = position_dodge(width = 0.6)) +
    scale_y_continuous('Cohen\'s d,\nrelative to log(TP10K)') +
    ggtitle('Cohen\'s d, false discoveries') +
    scale_color_manual(values = color_pal,labels = parse_format()) +
    scale_fill_manual(values = color_pal,labels = parse_format()) +
    boxed_theme() +
    theme(aspect.ratio = 1, 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size=6),
          legend.position = 'bottom',
          legend.text.align = 0,
          legend.title = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), 'cm')
    ) +
    guides(color=guide_legend(ncol=1))
p11
ggsave("fig/Fig6/splatter-false-discoveries-normalization-d.pdf", p11, width = 3.8, height = 4.6, units = "cm", useDingbats = FALSE)

###############################################################################-
## bulk normalization biases ####
###############################################################################-

# load data
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

###############################################################################-
## Mean expression, top-1000 ####
###############################################################################-

# average over bulk methods
avg1 = dat %>% 
    filter(min_cell_filter == 1, k == 1000, summary_type == 'mean') %>% 
    group_by_at(vars(-starts_with('mean_'), -starts_with('sd_'),
                     -starts_with('percentage_'), -starts_with('-peak_'))) %>% 
    summarise(mean_expr = mean(mean_expr),
              n = n()) %>% 
    ungroup()
delta = avg1 %>% 
    group_by_at(vars(-normalization, 
                     -mean_expr, -peak_width, -peak_width_rank, -color)) %>% 
    mutate(delta = mean_expr - mean_expr[normalization == 'logTP(10K)']) %>% 
    ungroup()
grid = distinct(delta, method, normalization) %>%
    filter(normalization != 'logTP(10K)')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    norm = as.character(current$normalization)
    df = filter(avg1, method == current$method, 
                normalization %in% c(norm, 'logTP(10K)')) %>% 
        droplevels() %>% 
        mutate(normalization = factor(normalization, levels = c(norm, 'logTP(10K)')))
    eff = cohen.d(formula = mean_expr ~ normalization | Subject(cell_type),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})
p12 = eff %>%
    ggplot(aes(x = reorder(normalization, -d, mean),
               color = method, fill = method, y = d)) +
    geom_hline(aes(yintercept = 0), linetype = 'solid', size = 0.4,
               color = 'grey82') +
    geom_point(position = position_dodge(width = 0.6), shape = 21, size = 0.8) + 
    geom_errorbar(aes(ymin = 0, ymax = d, group = method), width = 0, size = 0.3,
                  position = position_dodge(width = 0.6)) +
    scale_y_continuous('Cohen\'s d,\nrelative to log(TP10K)') +
    ggtitle('Cohen\'s d, read depth') +
    scale_color_manual(values = color_pal,labels = parse_format()) +
    scale_fill_manual(values = color_pal,labels = parse_format()) +
    boxed_theme() +
    theme(aspect.ratio = 1, 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size=6),
          legend.position = 'bottom',
          legend.text.align = 0,
          legend.title = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), 'cm')
    ) +
    guides(color=guide_legend(ncol=1))
p12
ggsave("fig/Fig6/bulk-mean-expr-normalization-d.pdf", p12, width = 3.8, height = 4.6, units = "cm", useDingbats = FALSE)

###############################################################################-
## Percentage open, top-1000 ####
###############################################################################-

# average over bulk methods
avg2 = dat %>% 
    filter(min_cell_filter == 1, k == 1000, summary_type == 'mean') %>% 
    group_by_at(vars(-starts_with('mean_'), -starts_with('sd_'),
                     -starts_with('percentage_'), -starts_with('-peak_'))) %>% 
    summarise(percentage_open = mean(percentage_open),
              n = n()) %>% 
    ungroup()
delta = avg2 %>% 
    group_by_at(vars(-normalization, -percentage_open, -peak_width, -peak_width_rank, -color)) %>% 
    mutate(delta = percentage_open - percentage_open[normalization == 'logTP(10K)']) %>% 
    ungroup()

grid = distinct(delta, method, normalization) %>%
    filter(normalization != 'logTP(10K)')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    norm = as.character(current$normalization)
    df = filter(avg2, method == current$method, 
                normalization %in% c(norm, 'logTP(10K)')) %>% 
        droplevels() %>% 
        mutate(normalization = factor(normalization, levels = c(norm, 'logTP(10K)')))
    eff = cohen.d(formula = percentage_open ~ normalization | Subject(cell_type),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})
p13 = eff %>%
    ggplot(aes(x = reorder(normalization, -d, mean),
               color = method, fill = method, y = d)) +
    geom_hline(aes(yintercept = 0), linetype = 'solid', size = 0.4,
               color = 'grey82') +
    geom_point(position = position_dodge(width = 0.6), shape = 21, size = 0.8) + 
    geom_errorbar(aes(ymin = 0, ymax = d, group = method), width = 0, size = 0.3,
                  position = position_dodge(width = 0.6)) +
    scale_y_continuous('Cohen\'s d,\nrelative to log(TP10K)') +
    ggtitle('Cohen\'s d, % open') +
    scale_color_manual(values = color_pal,labels = parse_format()) +
    scale_fill_manual(values = color_pal,labels = parse_format()) +
    boxed_theme() +
    theme(aspect.ratio = 1, 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size=6),
          legend.position = 'bottom',
          legend.text.align = 0,
          legend.title = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), 'cm')
    ) +
    guides(color=guide_legend(ncol=1))
p13
ggsave("fig/Fig6/bulk-pct-open-normalization-d.pdf", p13, width = 3.8, height = 4.6, units = "cm", useDingbats = FALSE)

# average over bulk methods
avg3 = dat %>% 
    filter(min_cell_filter == 1, k == 1000, summary_type == 'mean') %>% 
    group_by_at(vars(-starts_with('mean_'), -starts_with('sd_'),
                     -starts_with('percentage_'), -starts_with('-peak_'))) %>% 
    summarise(peak_width = mean(peak_width),
              n = n()) %>% 
    ungroup()
delta = avg3 %>% 
    group_by_at(vars(-normalization, -peak_width, -peak_width_rank, -color)) %>% 
    mutate(delta = peak_width - peak_width[normalization == 'logTP(10K)']) %>% 
    ungroup()
grid = distinct(delta, method, normalization) %>%
    filter(normalization != 'logTP(10K)')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    norm = as.character(current$normalization)
    df = filter(avg3, method == current$method, 
                normalization %in% c(norm, 'logTP(10K)')) %>% 
        droplevels() %>% 
        mutate(normalization = factor(normalization, levels = c(norm, 'logTP(10K)')))
    eff = cohen.d(formula = peak_width ~ normalization | Subject(cell_type),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})
p14 = eff %>%
    ggplot(aes(x = reorder(normalization, -d, mean),
               color = method, fill = method, y = d)) +
    geom_hline(aes(yintercept = 0), linetype = 'solid', size = 0.4,
               color = 'grey82') +
    geom_point(position = position_dodge(width = 0.6), shape = 21, size = 0.8) + 
    geom_errorbar(aes(ymin = 0, ymax = d, group = method), width = 0, size = 0.3,
                  position = position_dodge(width = 0.6)) +
    scale_y_continuous('Cohen\'s d,\nrelative to log(TP10K)') +
    ggtitle('Cohen\'s d, peak width (bp)') +
    scale_color_manual(values = color_pal,labels = parse_format()) +
    scale_fill_manual(values = color_pal,labels = parse_format()) +
    boxed_theme() +
    theme(aspect.ratio = 1, 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size=6),
          legend.position = 'bottom',
          legend.text.align = 0,
          legend.title = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), 'cm')
    ) +
    guides(color=guide_legend(ncol=1))
p14
ggsave("fig/Fig6/bulk-peak-width-normalization-d.pdf", p14, width = 3.8, height = 4.6, units = "cm", useDingbats = FALSE)


###############################################################################-
## ArchR cell counts plot ####
###############################################################################-

# Arch R cell counts
cell_counts = readRDS('data/summaries/misc/comparisons_with_ArchR_counts.rds')
cell_counts %<>%
    filter(!is.na(archR_cells))

dat0 = rbind(
    cell_counts %>% 
        dplyr::select(-archR_cells) %>%
        mutate(type = 'default'),
    cell_counts %>% 
        dplyr::select(-cell_count) %>%
        dplyr::rename(cell_count = archR_cells) %>%
        mutate(type = 'background matching')
    )

labs = dat0 %>% 
    group_by(type) %>% 
    summarise(text_val = round(mean(cell_count), 0),cell_count=14) %>%
    ungroup() %>%
    mutate(text_val = formatC(text_val, format="d", big.mark=","))

pal2 = c('default' = 'black',
         'background matching' = 'grey75')
p15 = dat0 %>%
    mutate(
        type = factor(type, levels=c('default', 'background matching')),
        cell_count = cell_count/1000
        ) %>%
    ggplot(aes(x=type, y=cell_count, color=type, fill=type)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
    geom_label(dat = labs,
               aes(label = text_val), color = 'black',
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0.5, vjust = 0,
               show.legend = FALSE) +
    scale_fill_manual(values = pal2) +
    scale_color_manual(values = pal2) +
    boxed_theme() +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 1.5,
        legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle=45, hjust=1)
    ) +
    ylim(0, 15) +
    ylab(expression('# of cells per comparison ('~10^3~')'))
p15
ggsave('fig/Fig6/ArchR-downsample-cells.pdf', p15, width=3, height=8, units='cm')


###############################################################################-
## ArchR AUCC ####
###############################################################################-
# load data
dat = readRDS("data/summaries/concordance_res/matched_bulk/all_aucc_ArchR_bg.rds") %>% 
    type_convert() %>% 
    # here, we are not interested in terciles
    filter(expr == 'all')
keep = map_lgl(dat, ~ n_distinct(.) > 1)
dat %<>% extract(, keep)
avg = dat %>% 
    group_by_at(vars(-bulk_da_method, -bulk_da_type,
                     -bulk_features, -overlap, -aucc)) %>% 
    summarise(aucc = mean(aucc),
              n = n()) %>% 
    ungroup()

# re-code x-values and colors
avg %<>%
    mutate(method = paste0(sc_da_method, '-', sc_da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .),
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', sc_da_family) %>% 
               fct_recode('single-cell' = 'singlecell',
                          'other' = 'non_libra') %>% 
               fct_relevel('single-cell', 'pseudobulk', 'other'))
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
                               'LR[peaks]' = 'LR_peaks',
                               'LR[clusters]' = 'LR',
                               'Binomial' = 'binomial',
                               'Wilcoxon~rank-sum~test' = 'wilcox',
                               'Fisher~exact~test' = 'FET',
                               'Negative~binomial' = 'negbinom') %>% 
               as.character()) %>% 
    # drop limma
    filter(!grepl('limma', method))
avg2 = avg %>% 
    filter(!grepl('pseudo-replicates', method)) %>% 
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    filter(k == 5000, filter == '1 cell')
# load data
dat = readRDS("data/summaries/concordance_res/matched_bulk/all_aucc_skinny.rds") %>% 
    type_convert() %>% 
    # here, we are not interested in terciles
    filter(expr == 'all')
# average over bulk methods
avg = dat %>% 
    group_by_at(vars(-bulk_da_family, -bulk_da_method, -bulk_da_type,
                     -bulk_features, -overlap, -aucc)) %>% 
    summarise(aucc = mean(aucc),
              n = n()) %>% 
    ungroup()

# re-code x-values and colors
avg %<>%
    mutate(method = paste0(sc_da_method, '-', sc_da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = ifelse(sc_pseudo_repl,
                           paste(method, '(pseudo-replicates)'), method),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .),
           color = ifelse(sc_pseudo_repl, 'pseudo-replicates', sc_da_family),
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', color) %>% 
               fct_recode('single-cell' = 'singlecell',
                          'other' = 'non_libra') %>% 
               fct_relevel('single-cell', 'pseudobulk', 'other'))
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
                               'LR[peaks]' = 'LR_peaks',
                               'LR[clusters]' = 'LR',
                               'Binomial' = 'binomial',
                               'Wilcoxon~rank-sum~test' = 'wilcox',
                               'Fisher~exact~test' = 'FET',
                               'Negative~binomial' = 'negbinom') %>% 
               as.character()) %>% 
    # drop limma
    filter(!grepl('limma', method))

# extract parameters of interest
avg1 = avg %>% 
    filter(!grepl('pseudo-replicates', method)) %>% 
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    filter(k == 5000, filter == '1 cell')

# combine
comb = bind_rows(avg1 %>% mutate(mode = 'default'),
                 avg2 %>% mutate(mode = 'background matching'))
lvls = avg1 %$% 
    reorder(method, aucc, stats::median) %>% 
    levels()
pal2 = c('default' = 'black',
         'background matching' = 'grey75')
p16 = comb %>%
    ggplot(aes(x = factor(method, levels = lvls),
               color = mode, fill = mode)) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    # geom_text(data = labs, aes(y = -0.12, label = label), size = 1.75, hjust = 0,
    #           color = 'black') +
    scale_color_manual('', values = pal2) +
    scale_fill_manual('', values = pal2) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('AUCC', breaks = seq(0, 1, 0.1)) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'right',
          legend.justification = 'bottom',
          legend.key.width = unit(0.4, 'lines'),
          aspect.ratio = 1.7)
p16
ggsave("fig/Fig6/ArchR-downsample-aucc.pdf", p16, width = 7, height = 6, units = "cm", useDingbats = FALSE)

delta = comb %>% 
    group_by(paper, cell_type, method) %>% 
    mutate(delta = aucc[mode == 'background matching'] - 
               aucc[mode == 'default'],
           n = n()) %>% 
    ungroup()

grid = distinct(delta, method, color)
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    print(current)
    df = filter(comb, method == current$method) %>% 
        mutate(mode = factor(mode, levels = c('background matching', 'default')))
    eff = cohen.d(formula = aucc ~ mode | Subject(paste(paper, cell_type)),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})
pal = da_analysis_colors
p17 = eff %>%
    ggplot(aes(
        x = factor(method, levels = lvls),
        # x = reorder(method, -d, mean, na.rm = TRUE),
        color = method, fill = method, y = d)) +
    geom_hline(aes(yintercept = 0), linetype = 'solid', size = 0.4,
               color = 'grey82') +
    geom_point(position = position_dodge(width = 0.6), shape = 21, size = 0.8) + 
    geom_errorbar(aes(ymin = 0, ymax = d, group = method), width = 0, size = 0.3,
                  position = position_dodge(width = 0.6)) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('Cohen\'s d') +
    coord_flip() +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    boxed_theme() +
    theme(aspect.ratio = 1, 
          axis.title.y = element_blank(),
          legend.position = 'none'
          # axis.text.x = element_text(angle = 45, hjust = 1)
    )
p17
ggsave("fig/Fig6/aucc-ArchR-background-matching-d.pdf", p17,
       width = 5, height = 5, units = "cm", useDingbats = FALSE)

top_row = wrap_plots(
    p1 + theme(legend.position = 'none', aspect.ratio = 1.8, plot.margin = unit(c(0.01, 0.04, 0.01, 0.01), "cm")),
    p2 + theme(legend.position = 'none', aspect.ratio = 1.8, plot.margin = unit(c(0.01, 0.04, 0.01, 0.01), "cm")),
    p3 + theme(legend.position = 'none', aspect.ratio = 1.8, plot.margin = unit(c(0.01, 0.04, 0.01, 0.01), "cm")),
    p4 + theme(legend.position = 'none', aspect.ratio = 1.8, plot.margin = unit(c(0.01, 0.04, 0.01, 0.01), "cm")),
    nrow = 1
)
# top_row
ggsave("fig/Fig6/top-row.pdf", top_row, width = 18, height = 15, units = "cm", useDingbats = FALSE)

middle_row = wrap_plots(
    p5 + theme(legend.position = 'none', aspect.ratio = 1.8, plot.margin = unit(c(0.01, 0.04, 0.01, 0.01), "cm")),
    p6 + theme(legend.position = 'none', aspect.ratio = 1.8, plot.margin = unit(c(0.01, 0.04, 0.01, 0.01), "cm")),
    p7 + theme(legend.position = 'none', aspect.ratio = 1.8, plot.margin = unit(c(0.01, 0.04, 0.01, 0.01), "cm")),
    nrow = 1
)
ggsave("fig/Fig6/middle-row.pdf", middle_row, width = 13.5, height = 15, units = "cm", useDingbats = FALSE)

bottom_row = wrap_plots(
    p9 + theme(legend.position = 'none', plot.margin = unit(c(0.01, 0.4, 0.01, 0.01), "cm")),
    p10 + theme(legend.position = 'none', plot.margin = unit(c(0.01, 0.4, 0.01, 0.01), "cm"), axis.title.y = element_blank()),
    p11 + theme(legend.position = 'none', plot.margin = unit(c(0.01, 0.4, 0.01, 0.01), "cm"), axis.title.y = element_blank()),
    p12 + theme(legend.position = 'none', plot.margin = unit(c(0.01, 0.4, 0.01, 0.01), "cm"), axis.title.y = element_blank()),
    p13 + theme(legend.position = 'none', plot.margin = unit(c(0.01, 0.4, 0.01, 0.01), "cm"), axis.title.y = element_blank()),
    p14 + theme(legend.position = 'none', plot.margin = unit(c(0.01, 0.4, 0.01, 0.01), "cm"), axis.title.y = element_blank()),
    ncol = 6
) 
ggsave("fig/Fig6/bottom-row.pdf", bottom_row, width = 18.3, height = 12, units = "cm", useDingbats = FALSE)
    
