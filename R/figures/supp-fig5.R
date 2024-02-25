# Plot a summary figure with performance across all experiments.
setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(pals)
library(upstartr)
source('R/theme.R')

################-

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

# average over bulk methods
avg2 = dat %>% 
    filter(k %in% c(500, 5000),
           summary_type == 'mean') %>% 
    group_by_at(vars(-starts_with('mean_'), -starts_with('sd_'),
                     -starts_with('percentage_'), -starts_with('-peak_'))) %>% 
    summarise(mean_expr = mean(mean_expr),
              n = n()) %>% 
    ungroup()
table(avg2$n, avg2$k)
dplyr::count(avg2, da_family, da_method, da_type)

# plot
labs = avg2 %>% 
    group_by(k, method, color) %>%
    summarise(mean = mean(mean_expr), n = n(),
              median = median(mean_expr), mean_expr = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
pal = da_analysis_colors
p1 = map(c(500, 5000), ~ {
    avg2 %>%
        filter(k == .x) %>% 
        ggplot(aes(x = reorder(method, mean_expr, stats::median), 
                   color = method, fill = method)) +
        ggtitle(paste('k =', .x)) +
        geom_boxplot(aes(y = mean_expr), outlier.shape = NA, alpha = 0.4, width = 0.6,
                     size = 0.4) + 
        geom_text(data = filter(labs, k == .x),
                  aes(y = ifelse(.x == 500, -3.5, -1), label = label), size = 1.75, hjust = 0,
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
}) %>% wrap_plots(nrow = 1)
p1
ggsave("fig/Supp_Fig5/bulk-mean-expr-k.pdf", p1,
       width = 14, height = 6.5, units = "cm", useDingbats = FALSE)


# average over bulk methods
avg4 = dat %>% 
    filter(k %in% c(500, 5000),
           summary_type == 'mean') %>% 
    group_by_at(vars(-starts_with('mean_'), -starts_with('sd_'),
                     -starts_with('percentage_'), -starts_with('-peak_'))) %>% 
    summarise(percentage_open = mean(percentage_open),
              n = n()) %>% 
    ungroup()
table(avg4$n, avg4$k)
dplyr::count(avg4, da_family, da_method, da_type)

# plot
labs = avg4 %>% 
    group_by(k, method, color) %>%
    summarise(mean = mean(percentage_open), n = n(),
              median = median(percentage_open), percentage_open = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
p2 = map(c(500, 5000), ~ {
    avg4 %>%
        filter(k == .x) %>% 
        ggplot(aes(x = reorder(method, percentage_open, stats::median), 
                   color = method, fill = method)) +
        ggtitle(paste('k =', .x)) +
        geom_boxplot(aes(y = percentage_open), outlier.shape = NA, alpha = 0.4,
                     width = 0.6, size = 0.4) + 
        geom_text(data = filter(labs, k == .x),
                  aes(y = -0.1, label = label), size = 1.75, hjust = 0,
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
}) %>% wrap_plots(nrow = 1)
p2
ggsave("fig/Supp_Fig5/bulk-pct-open-k.pdf", p2, width = 14, height = 6.5, units = "cm", useDingbats = FALSE)

# average over bulk methods
avg6 = dat %>% 
    filter(k %in% c(500, 5000),
           summary_type == 'mean') %>% 
    group_by_at(vars(-starts_with('mean_'), -starts_with('sd_'),
                     -starts_with('percentage_'), -starts_with('-peak_'))) %>% 
    summarise(peak_width = mean(peak_width),
              n = n()) %>% 
    ungroup()
table(avg6$n, avg6$k)
dplyr::count(avg6, da_family, da_method, da_type)

# plot
labs = avg6 %>% 
    group_by(k, method, color) %>%
    summarise(mean = mean(peak_width), n = n(),
              median = median(peak_width), peak_width = median) %>%
    ungroup() %>%
    mutate(label = round(median))
p3 = map(c(500, 5000), ~ {
    avg6 %>%
        filter(k == .x) %>% 
        ggplot(aes(x = reorder(method, peak_width, stats::median), 
                   color = method, fill = method)) +
        ggtitle(paste('k =', .x)) +
        geom_boxplot(aes(y = peak_width), outlier.shape = NA, alpha = 0.4,
                     width = 0.6, size = 0.4) + 
        geom_text(data = filter(labs, k == .x),
                  aes(y = -0.1, label = label), size = 1.75, hjust = 0,
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
}) %>% wrap_plots(nrow = 1)
p3
ggsave("fig/Supp_Fig5//bulk-peak-width-k.pdf", p3,
       width = 14, height = 6.5, units = "cm", useDingbats = FALSE)

full = wrap_plots(p1, p2, p3, ncol=2)
# full
ggsave("fig/Supp_Fig5/full_row.pdf", full,
       width = 18.5, height = 9, units = "cm", useDingbats = FALSE)

