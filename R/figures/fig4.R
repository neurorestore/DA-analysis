# Plot a summary figure with performance across all experiments.
setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

###############################################################################-
## bulk concordance biases plot ####
###############################################################################-

# load data
dat = readRDS("data/summaries/biases/matched_bulk_top_k_against_stats.rds") %>% 
    type_convert() %>% 
    mutate(paper = gsub('_scATAC_.*', '', basename(filename))) %>%
    dplyr::select(-filename, -seq_type) %>% 
    ## no pseudoreplicates
    filter(!pseudo_repl) %>% 
    ## ignore binarization here
    filter(!binarization) %>% 
    ## ignore normalization here
    filter(normalization == 'f') %>% 
    ## top-k filter must be none
    filter(top_k_filter == 'None') %>% 
    ## remove ground_truth_* cols
    dplyr::select(-ground_truth_da_method, -ground_truth_da_type) %>% 
    ## distinct
    distinct() %>% 
    # filter == 0
    filter(min_cell_filter == 1, !percent_filter)

# remove non-unique columns
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

# remove limma and some Gonzalez-Blas
dat %<>%
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>%
    # filter(!(grepl('Melanoma cell line-24 hours', cell_type) | grepl('Melanoma cell line-48 hours', cell_type))) %>% 
    filter(!grepl('limma', da_method))

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-', da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = da_family,
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
# print 
dplyr::count(dat, color, method) %>% arrange(color, method)

###############################################################################-
## Mean expression, top-1000 ####
###############################################################################-

# average over bulk methods
avg1 = dat %>% 
    filter(k == 1000, summary_type == 'mean') %>% 
    group_by_at(vars(-starts_with('mean_'), -starts_with('sd_'),
                     -starts_with('percentage_'), -starts_with('-peak_'))) %>% 
    summarise(mean_expr = mean(mean_expr),
              n = n()) %>% 
    ungroup()
table(avg1$n)
dplyr::count(avg1, da_family, da_method, da_type)

# plot
labs = avg1 %>% 
    group_by(k, method, color) %>%
    summarise(mean = mean(mean_expr), n = n(),
              median = median(mean_expr), mean_expr = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
saveRDS(labs, 'data/summaries/meta_summaries/bulk_aucc_mean_expr_bias_summary.rds')

pal = da_analysis_colors
p1 = avg1 %>%
    ggplot(aes(x = reorder(method, mean_expr, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = mean_expr), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    geom_text(data = labs, aes(y = -1.5, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('Read depth') +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          aspect.ratio = 1.4)
p1
ggsave("fig/Fig4/bulk-bias-mean.pdf", p1,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)

# average over bulk methods
avg3 = dat %>% 
    filter(k == 1000, summary_type == 'mean') %>% 
    group_by_at(vars(-starts_with('mean_'), -starts_with('sd_'),
                     -starts_with('percentage_'), -starts_with('-peak_'))) %>% 
    summarise(percentage_open = mean(percentage_open),
              n = n()) %>% 
    ungroup()
table(avg3$n)
dplyr::count(avg3, da_family, da_method, da_type)

# plot
labs = avg3 %>% 
    group_by(k, method, color) %>%
    summarise(mean = mean(percentage_open), n = n(),
              median = median(percentage_open), percentage_open = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
saveRDS(labs, 'data/summaries/meta_summaries/bulk_aucc_pct_open_bias_summary.rds')

p2 = avg3 %>%
    ggplot(aes(x = reorder(method, percentage_open, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = percentage_open), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    geom_text(data = labs, aes(y = -0.05, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('% open', labels = ~ . * 100) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          aspect.ratio = 1.4)
p2
ggsave("fig/Fig4/bulk-bias-pct-open.pdf", p2,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)

# average over bulk methods
avg5 = dat %>% 
    filter(k == 1000, summary_type == 'mean') %>% 
    group_by_at(vars(-starts_with('mean_'), -starts_with('sd_'),
                     -starts_with('percentage_'), -starts_with('-peak_'))) %>% 
    summarise(peak_width = mean(peak_width),
              n = n()) %>% 
    ungroup()
table(avg5$n)
dplyr::count(avg5, da_family, da_method, da_type)

# plot
labs = avg5 %>% 
    group_by(k, method, color) %>%
    summarise(mean = mean(peak_width), n = n(),
              median = median(peak_width), peak_width = median) %>%
    ungroup() %>%
    mutate(label = round(median))
saveRDS(labs, 'data/summaries/meta_summaries/bulk_aucc_peak_width_bias_summary.rds')

p3 = avg5 %>%
    ggplot(aes(x = reorder(method, peak_width, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = peak_width), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    geom_text(data = labs, aes(y = -0.05, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('Peak width (bp)') +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          aspect.ratio = 1.4)
p3
ggsave("fig/Fig4/bulk-bias-peak-width.pdf", p3,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)

p4 = p1 | p2 + ggtitle('Characteristics of top-1,000 DA peaks from each DA method') + theme(plot.title = element_text(size=6, margin=margin(0,0,20,0))) | p3
ggsave('fig/Fig4/row.pdf', p4, width = 15.5, height=6.5, units='cm')
