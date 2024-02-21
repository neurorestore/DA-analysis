# Plot a summary figure with performance across all experiments.
setwd("C:/Users/teo/Documents/EPFL/projects/DA-analysis/")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

###############################################################################-
## Luecken false discoveries plot ####
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
    filter(is.na(percentile), is.na(latent_vars), normalization=='f', !binarization)

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-', da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = ifelse(pseudo_repl,
                           paste(method, '(pseudo-replicates)'), method),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(pseudo_repl, 'pseudo-replicates', da_family),
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

labs = dat %>% 
    group_by(method, color) %>%
    summarise(mean = mean(number_of_da_regions), 
              n = n(), 
              median = median(number_of_da_regions), 
              number_of_da_regions = median) %>%
    ungroup() %>%
    mutate(label = round(median))
saveRDS(labs, 'data/summaries/meta_summaries/luecken_false_discoveries_summary.rds')

pal = da_analysis_colors

p1 = dat %>%
    ggplot(aes(x = reorder(method, number_of_da_regions, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, 
                 alpha = 0.4, width = 0.6, size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -12.5e3, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous(expression('# of DA peaks'~(10^3)), 
                       labels = function(x) x / 1e3,
                       breaks = seq(0, 80, 10) * 1e3) +
    coord_flip(ylim = c(-12e3, 55e3)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    ggtitle('False discoveries in null comparisons\nof real scATAC-seq data') +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none',
          plot.title = element_text(size=5,hjust = 0.5)
      )
p1
ggsave("fig/Fig3/luecken-false-discoveries.pdf", p1,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)


###############################################################################-
## Snapatac false discoveries plot ####
###############################################################################-

# load data
dat = readRDS("data/summaries/false_discoveries/false_discoveries_snapatac.rds") %>% 
    type_convert() %>% 
    # dplyr::select(-filename) %>% 
    # filter out methods we are not using
    filter(da_family != 'mixedmodel',
           pseudo_repl == FALSE) %>% 
    # entire range of expression values only - no binarization
    filter(is.na(percentile), !binarization, normalization == 'f')

## DA methods
dplyr::count(dat, da_family, da_method, de_type, pseudobulk, pseudo_repl)

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-', de_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = ifelse(pseudo_repl,
                           paste(method, '(pseudo-replicates)'), method),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(pseudo_repl, 'pseudo-replicates', da_family),
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

labs = dat %>% 
    group_by(method, color) %>%
    summarise(mean = mean(number_of_da_regions), 
              n = n(), 
              median = median(number_of_da_regions), 
              number_of_da_regions = median) %>%
    ungroup() %>%
    mutate(label = round(median))
saveRDS(labs, 'data/summaries/meta_summaries/snapatac_false_discoveries_summary.rds')

pal = c('mixed model' = 'grey80',
        'single-cell' = colours.cafe322[1],
        'pseudobulk' = colours.cafe322[2],
        'pseudo-replicates' = "#6E645F",
        'other' = colours.cafe322[3])
pal = da_analysis_colors

p2 = dat %>%
    ggplot(aes(x = reorder(method, number_of_da_regions, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, 
                 alpha = 0.4, width = 0.6, size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -15.5e3, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous(expression('# of DA peaks'~(10^3)), 
                       labels = function(x) x / 1e3,
                       breaks = seq(0, 80e3, 10e3)) +
    coord_flip(ylim = c(-15e3, 75e3)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    ggtitle('False discoveries in downsampled\nbulk ATAC-seq data') +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none',
          plot.title = element_text(size=5, hjust=0.5)
      )
p2
ggsave("fig/Fig3/snapatac-false-discoveries.pdf", p2,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)

###############################################################################-
## Simulated false discoveries plot ####
###############################################################################-

dat = readRDS("data/summaries/false_discoveries/false_discoveries_simulated.rds") %>% 
    type_convert() %>% 
    dplyr::select(-ori_filename) %>% 
    ## remove libra bug
    dplyr::select(-de_family, -de_method) %>% 
    dplyr::rename(da_type = de_type) %>%
    ## cell type is actually simulation replicate
    dplyr::rename(replicate = cell_type) %>% 
    # remove pseudoreplicates and mixed models
    filter(!pseudo_repl, da_family != 'mixedmodel') %>% 
    # filter to new grid
    filter(de_loc >= 0.5) %>%
    filter(is.na(percentile)) %>%
    filter(!binarization, normalization=='f')
    
# all 5% FDR# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-', da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = ifelse(pseudo_repl,
                           paste(method, '(pseudo-replicates)'), method),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(pseudo_repl, 'pseudo-replicates', da_family),
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

## defaults
dat %<>%  filter(de_loc == 1, de_prob == 0.5, n_reps == 6, n_cells == 1000)
labs = dat %>%
    group_by(method, color, n_cells) %>%
    summarise(mean = mean(number_of_da_regions), 
              n = n(), 
              median = median(number_of_da_regions), 
              number_of_da_regions = median) %>%
    ungroup() %>%
    mutate(label = round(median))
saveRDS(labs, 'data/summaries/meta_summaries/splatter_false_discoveries_summary.rds')

pal = da_analysis_colors
p3 = dat %>%
    ggplot(aes(x = reorder(method, number_of_da_regions, stats::median), 
               color = method, fill = method)) +
    # facet_wrap(~ n_cells, nrow = 1, 
    # labeller = as_labeller(~ paste0('n_cells=', .))) +
    geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, 
                 alpha = 0.4, width = 0.6, size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs0, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -5.3e3, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    ggtitle('False discoveries in\nmodel-based simulations') + 
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous(expression('# of DA peaks'~(10^3)), 
                       labels = function(x) x / 1e3,
                       breaks = seq(0, 30, 10) * 1e3) +
    coord_flip(ylim = c(-5.2e3, 31e3)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none',
          plot.title = element_text(size=5, hjust=0.5)
    )
p3
ggsave("fig/Fig3/simulated-false-discoveries.pdf", p3,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)


###############################################################################-
## Simulated depth false discoveries plot ####
###############################################################################-
# load data
dat = readRDS("data/summaries/false_discoveries/false_discoveries_simulated_with_sequencing_depth.rds") %>% 
    type_convert() %>% 
    dplyr::select(-ori_filename) %>%
    # filter out methods we are not using
    filter(da_family != 'mixedmodel',
           pseudo_repl == FALSE) %>% 
    # drop percentiles
    filter(is.na(percentile),
           normalization == 'f',
           !binarization)

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-', de_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = ifelse(pseudo_repl,
                           paste(method, '(pseudo-replicates)'), method),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(pseudo_repl, 'pseudo-replicates', da_family),
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

labs = dat %>% 
    group_by(method, color) %>%
    summarise(mean = mean(number_of_da_regions), 
              n = n(), 
              median = median(number_of_da_regions), 
              number_of_da_regions = median) %>%
    ungroup() %>%
    mutate(label = round(median))
pal = c('mixed model' = 'grey80',
        'single-cell' = colours.cafe322[1],
        'pseudobulk' = colours.cafe322[2],
        'pseudo-replicates' = "#6E645F",
        'other' = colours.cafe322[3])
pal = da_analysis_colors
p4 = dat %>%
    ggplot(aes(x = reorder(method, number_of_da_regions, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, 
                 alpha = 0.4, width = 0.6, size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -10.5e3, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous(expression('# of DA peaks'~(10^3)), 
                       labels = function(x) x / 1e3,
                       breaks = seq(0, 80, 10) * 1e3) +
    coord_flip(ylim = c(-10e3, 62e3)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    ggtitle('False discoveries in model-based\nsimulations with variable sequencing depth') + 
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none',
          plot.title = element_text(size=5, hjust=0.5)
    )
p4
ggsave("fig/Fig3/simulated-depth-false-discoveries.pdf", p4,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)

p5 = p2 | p3 | p4 
ggsave("fig/Fig3/bottom-row.pdf", p5,
       width = 15.5, height = 6.5, units = "cm", useDingbats = FALSE)

