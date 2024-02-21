# Plot a summary figure with performance across all experiments.
setwd("C:/Users/teo/Documents/EPFL/projects/DA-analysis/")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(pals)
library(upstartr)
library(effsize)
source('R/theme.R')

# load data
dat = readRDS("data/summaries/biases/matched_bulk_top_k_against_stats.rds") %>% 
    type_convert() %>% 
    mutate(paper = gsub('_scATAC_.*', '', basename(filename))) %>%
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
    distinct() %>%
    filter(min_cell_filter == 1, !percent_filter)
# remove non-unique columns
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

# remove limma and some Gonzalez-Blas
dat %<>%
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    filter(!grepl('limma', da_method))

# print some summaries
## all results on peaks
all(dat$binsize == 'peak') ## TRUE
## peak filter
table(dat$min_cell_percent_filter)
## DA methods
dplyr::count(dat, da_family, da_method, da_type, pseudobulk)
## datasets
dplyr::count(dat, paper, cell_type, paper_da_type)
## k
table(dat$k)
## summary type
table(dat$summary_type)

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
# print 
dplyr::count(dat, color, method) %>% arrange(color, method)

# filter to only methods that _can_ be binarized
dat %<>%
    group_by(method) %>% 
    filter(n_distinct(color) == 2) %>%
    ungroup()
dplyr::count(dat, color, method) %>% arrange(method, color)

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

delta = avg1 %>% 
    group_by_at(vars(-mean_expr, -peak_width, -peak_width_rank, -color, -n, -binarization)) %>% 
    mutate(delta = mean_expr[binarization] - mean_expr[!binarization]) %>% 
    ungroup()

grid = distinct(delta, method, color) %>%
    filter(color != 'binarized')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    df = filter(avg1, method == current$method) %>% 
        droplevels() %>% 
        mutate(binarization = factor(binarization, levels = c(T, F)))
    eff = cohen.d(formula = mean_expr ~ binarization | Subject(cell_type),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})

pal = da_analysis_colors
p1 = eff %>%
    ggplot(aes(x = reorder(method, -d, mean),
               color = method, fill = method, y = d)) +
    geom_hline(aes(yintercept = 0), linetype = 'solid', size = 0.4,
               color = 'grey82') +
    geom_point(position = position_dodge(width = 0.6), shape = 21, size = 0.8) + 
    geom_errorbar(aes(ymin = 0, ymax = d, group = method), width = 0, size = 0.3,
                  position = position_dodge(width = 0.6)) +
    scale_y_continuous('Cohen\'s d') +
    
    boxed_theme() +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    theme(aspect.ratio = 1, 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'none')
p1
ggsave("fig/Supp_Fig13/bulk-mean-expr-binarization-d.pdf", p3,
       width = 5.5, height = 5.5, units = "cm", useDingbats = FALSE)

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
delta = avg3 %>% 
    group_by_at(vars(-percentage_open, -peak_width, -peak_width_rank, -color, -n, -binarization)) %>% 
    mutate(delta = percentage_open[binarization] - percentage_open[!binarization]) %>% 
    ungroup()

grid = distinct(delta, method, color) %>%
    filter(color != 'binarized')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    df = filter(avg3, method == current$method) %>% 
        droplevels() %>% 
        mutate(binarization = factor(binarization, levels = c(T, F)))
    eff = cohen.d(formula = percentage_open ~ binarization | Subject(cell_type),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})

pal = da_analysis_colors
p2 = eff %>%
    ggplot(aes(x = reorder(method, -d, mean),
               color = method, fill = method, y = d)) +
    geom_hline(aes(yintercept = 0), linetype = 'solid', size = 0.4,
               color = 'grey82') +
    geom_point(position = position_dodge(width = 0.6), shape = 21, size = 0.8) + 
    geom_errorbar(aes(ymin = 0, ymax = d, group = method), width = 0, size = 0.3,
                  position = position_dodge(width = 0.6)) +
    scale_y_continuous('Cohen\'s d') +
    
    boxed_theme() +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    theme(aspect.ratio = 1, 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'none')
p2
ggsave("fig/Supp_Fig13/bulk-pct-open-binarization-d.pdf", p2,
       width = 5.5, height = 5.5, units = "cm", useDingbats = FALSE)

# average over bulk methods
avg5 = dat %>% 
    filter(k == 1000, summary_type == 'mean') %>% 
    group_by_at(vars(-starts_with('mean_'), -starts_with('sd_'),
                     -starts_with('percentage_'), -starts_with('-peak_'))) %>% 
    summarise(peak_width = mean(peak_width),
              n = n()) %>% 
    ungroup()

delta = avg5 %>% 
    group_by_at(vars(-peak_width, -peak_width_rank, -color, -n, -binarization)) %>% 
    mutate(delta = peak_width[binarization] - peak_width[!binarization]) %>% 
    ungroup()

grid = distinct(delta, method, color) %>%
    filter(color != 'binarized')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    df = filter(avg5, method == current$method) %>% 
        droplevels() %>% 
        mutate(binarization = factor(binarization, levels = c(T, F)))
    eff = cohen.d(formula = peak_width ~ binarization | Subject(cell_type),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})

pal = da_analysis_colors
p3 = eff %>%
    ggplot(aes(x = reorder(method, -d, mean),
               color = method, fill = method, y = d)) +
    geom_hline(aes(yintercept = 0), linetype = 'solid', size = 0.4,
               color = 'grey82') +
    geom_point(position = position_dodge(width = 0.6), shape = 21, size = 0.8) + 
    geom_errorbar(aes(ymin = 0, ymax = d, group = method), width = 0, size = 0.3,
                  position = position_dodge(width = 0.6)) +
    scale_y_continuous('Cohen\'s d') +
    
    boxed_theme() +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    theme(aspect.ratio = 1, 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'none')
p3
ggsave("fig/Supp_Fig13/bulk-peak-size-binarization-d.pdf", p3,
       width = 5.5, height = 5.5, units = "cm", useDingbats = FALSE)

full = wrap_plots(p1,p2,p3, nrow=1)
ggsave('fig/Supp_Fig13/full.pdf', full, width = 10, height = 12, units='cm')
