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
                               'LR[peaks]' = 'LR_peaks',
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
    filter(sc_min_cell_filter == 1, sc_percent_filter == FALSE,
           k == 5000) %>%
    group_by_at(vars(-bulk_da_method, -bulk_da_type,
                     -overlap, -aucc, -bulk_features, -sc_features)) %>% 
    summarise(aucc = mean(aucc),
              n = n()) %>% 
    ungroup()

delta = avg1 %>% 
    group_by_at(vars(-aucc, -color, -n, -sc_binarization)) %>% 
    mutate(delta = aucc[sc_binarization] - aucc[!sc_binarization]) %>% 
    ungroup()

delta_labs = delta %>% 
    filter(!sc_binarization) %>%
    group_by(method, color) %>%
    summarise(mean = mean(delta), n = n(), median = median(delta), 
              delta = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'e', digits = 1))
grid = distinct(delta, method, color) %>%
    filter(color != 'binarized')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    df = filter(avg1, method == current$method) %>% 
        droplevels() %>% 
        mutate(sc_binarization = factor(sc_binarization, levels = c(T, F)))
    eff = cohen.d(formula = aucc ~ sc_binarization | Subject(cell_type),
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
    scale_x_discrete(labels = scales::parse_format()) +
    boxed_theme() +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    theme(aspect.ratio = 1, 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'none'
          )
p1
ggsave("fig/Supp_Fig12/bulk-aucc-binarization-d.pdf", p1,
       width = 5.5, height = 5.5, units = "cm", useDingbats = FALSE)

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
table(dat$study) ## only Luecken peak matrix

## comparisons (study/label)
dplyr::count(dat, study, label) ## n=20
## bin sizes
table(dat$mapping_size, dat$mapping_type)
## DA methods
dplyr::count(dat, da_family, da_method, da_type, pseudo_repl, pseudobulk)

# remove non-unique columns
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)
# remove limma and some Gonzalez-Blas
dat %<>%
    filter(!grepl('limma', da_method))

# print some summaries
## single-cell DA
dplyr::count(dat, da_family, da_method, da_type)
## binarization
table(dat$binarization)

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

delta = dat %>% 
    group_by_at(vars(-number_of_da_regions, -number_of_nominal_da_regions, -number_of_regions, -color, -binarization)) %>% 
    mutate(delta = number_of_da_regions[binarization] - number_of_da_regions[!binarization]) %>% 
    ungroup()

grid = distinct(delta, method, color) %>%
    filter(color != 'binarized')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    df = filter(delta, method == current$method) %>% 
        droplevels() %>% 
        mutate(binarization = factor(binarization, levels = c(T, F)))
    eff = cohen.d(formula = number_of_da_regions ~ binarization | Subject(label),
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
          legend.position = 'none'
          )
p2
ggsave("fig/Supp_Fig12/false-discoveries-Luecken-binarization-d.pdf", p2,
       width = 5.5, height = 5.5, units = "cm", useDingbats = FALSE)

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

# print some summaries
## iteration
table(dat$seed)
## single-cell DA
dplyr::count(dat, da_family, da_method, pseudobulk)
dplyr::count(dat, de_family, de_method, de_type, pseudobulk, binarization)
## binarization
table(dat$binarization)

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
# print 
dplyr::count(dat, color, method) %>% arrange(color, method)

# filter to only methods that _can_ be binarized
dat %<>%
    group_by(method) %>% 
    filter(n_distinct(color) == 2) %>%
    ungroup()
dplyr::count(dat, color, method) %>% arrange(method, color)

delta = dat %>% 
    group_by_at(vars(-number_of_da_regions, -number_of_nominal_da_regions, -number_of_regions, -color, -binarization)) %>% 
    mutate(delta = number_of_da_regions[binarization] - number_of_da_regions[!binarization]) %>% 
    ungroup()

grid = distinct(delta, method, color) %>%
    filter(color != 'binarized')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    df = filter(delta, method == current$method) %>% 
        droplevels() %>% 
        mutate(binarization = factor(binarization, levels = c(T, F)))
    eff = cohen.d(formula = number_of_da_regions ~ binarization | Subject(seed),
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
ggsave("fig/Supp_Fig12/false-discoveries-snapatac-binarization-d.pdf", p3,
       width = 5.5, height = 5.5, units = "cm", useDingbats = FALSE)


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

# print some summaries
## iteration
table(dat$simulation_iter)
## single-cell DA
dplyr::count(dat, da_family, da_method, pseudobulk)
dplyr::count(dat, de_family, de_method, de_type, pseudobulk, binarization)
## binarization
table(dat$binarization)

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
# print 
dplyr::count(dat, color, method) %>% arrange(color, method)

# filter to only methods that _can_ be binarized
dat %<>%
    group_by(method) %>% 
    filter(n_distinct(color) == 2) %>%
    ungroup()
dplyr::count(dat, color, method) %>% arrange(method, color)

delta = dat %>% 
    group_by_at(vars(-number_of_da_regions, -number_of_nominal_da_regions, -number_of_regions, -color, -binarization)) %>% 
    mutate(delta = number_of_da_regions[binarization] - number_of_da_regions[!binarization]) %>% 
    ungroup()

grid = distinct(delta, method, color) %>%
    filter(color != 'binarized')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    df = filter(delta, method == current$method) %>% 
        droplevels() %>% 
        mutate(binarization = factor(binarization, levels = c(T, F)))
    eff = cohen.d(formula = number_of_da_regions ~ binarization | Subject(simulation_iter),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})

pal = da_analysis_colors
p4 = eff %>%
    mutate(d = ifelse(is.na(d), 0, d)) %>%
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
p4
ggsave("fig/Supp_Fig12/false-discoveries-splatter-binarization-d.pdf", p4,
       width = 5.5, height = 5.5, units = "cm", useDingbats = FALSE)

full = wrap_plots(p1,p2,p3,p4, nrow=1)
ggsave('fig/Supp_Fig12/full.pdf', full, width = 13, height = 12, units='cm')
