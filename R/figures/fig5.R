# Plot a summary figure with performance across all experiments.
setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(effsize)
library(upstartr)
source("R/theme.R")

###############################################################################-
## bulk logFC plot ####
###############################################################################-

# load data
dat = readRDS("data/summaries/concordance_res/matched_bulk/all_aucc_logFC_filter.rds") %>% 
    type_convert() %>% 
    # remove pseudoreplicates
    filter(!sc_pseudo_repl)
# remove fixed columns
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

# restructure 'expr' column
dat %<>% mutate(expr_set = ifelse(grepl('pseudo', expr), 'Random', 'logFC') %>% 
                    fct_relevel('Random', 'logFC'),
                expr = gsub("_pseudo", "", expr))

# average over bulk methods
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
    mutate(logFC = gsub("^.*_", "", expr) %>% 
               replace(. == 'all', '0.0')) %>%
    filter(!logFC %in% c('95%', '99%'))
enr = avg1 %>% 
    group_by_at(vars(-expr_set, -aucc)) %>% 
    summarise(delta = ifelse(n() == 2, 
                             aucc[expr_set == 'logFC'] - aucc[expr_set == 'Random'],
                             aucc),
              nn = n()) %>% 
    ungroup() 

grid = distinct(enr, method, logFC) %>% filter(logFC != '0.0')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    df = filter(avg1, method == current$method, logFC == current$logFC) %>% 
        mutate(expr_set = factor(expr_set, levels = c('logFC', 'Random')))
    eff = cohen.d(formula = aucc ~ expr_set | Subject(paste(paper, cell_type)),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})
labs = eff %>% 
    group_by(logFC) %>% 
    summarise(mean = mean(d), n = n(), median = median(d), d = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))

p1 = eff %>% 
    ggplot(aes(x = factor(logFC), y = d)) +
    geom_hline(aes(yintercept = 0), linetype = 'dotted', size = 0.25) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.6, size = 0.35,
                 fill = 'grey82', color = 'grey20') + 
    geom_label(data = labs, aes(y = Inf, label = label), size = 1.5, hjust = 0.5,
               vjust = 1, color = 'black', fill = NA, label.size = NA,
               label.padding = unit(0.6, 'lines')) +
    # scale_color_manual('', values = pal) +
    # scale_fill_manual('', values = pal) +
    scale_x_discrete('% of peaks removed') +
    scale_y_continuous('Cohen\'s d (vs. random peak removal)', limits = c(-0.65, 0.65)) +
    ggtitle('Effect of log-fold change filtering on concordance\nbetween scATAC-seq and matched bulk ATAC-seq') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = 1,
          plot.title = element_text(size=5, hjust=0.5)
          )
p1
ggsave("fig/Fig5/bulk-concordance-logFC-d.pdf", p1,
       width = 6, height = 6, units = "cm", useDingbats = FALSE)


###############################################################################-
## multiome logFC plot ####
###############################################################################-

# load data
dat = readRDS("data/summaries/concordance_res/multiome/all_aucc_logFC_filter.rds") %>% 
    type_convert() 
# no binarization or normalization
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

# restructure 'expr' column
dat %<>% mutate(expr_set = ifelse(grepl('pseudo', expr), 'Random', 'logFC'),
                expr = gsub("_pseudo", "", expr))
# filter Zhu, Anandon, atlas by label
dat %<>%
    filter(
        !grepl('Zhu', paper),
        !grepl('Anadon', paper),
        !(paper == 'Squair_2022' & grepl('7d_vs_Uninjured', cell_type)),
        !(paper == 'Squair_2022' & grepl('2m_vs_Uninjured', cell_type)),
    )

# remove mixed models and pseudoreplicates
dat %<>% filter(scatac_da_family != 'mixedmodel')

# average over bulk methods
avg = dat %>% 
    group_by_at(vars(-scrna_da_method, -bulk_da_type, -sc_features,
                     -bulk_features, -overlap, -aucc)) %>% 
    summarise(aucc = mean(aucc),
              n = n()) %>% 
    ungroup()

# re-code x-values and colors
avg %<>%
    mutate(method = paste0(scatac_da_method, '-', scatac_da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', scatac_da_family) %>% 
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

# re-code filtering
avg %<>%
    mutate(filter = ifelse(scatac_percent_filter, 
                           paste0(scatac_min_cell_filter, '%'),
                           paste0(scatac_min_cell_filter, ' cells')) %>% 
               fct_recode('1 cell' = '1 cells') %>% 
               fct_relevel('1 cell', '10 cells', '0.2%', '0.5%', '1%', '5%', 
                           '10%', '20%')) %>% 
    filter(!filter %in% c('10%', '20%')) 

avg1 = avg %>% 
    filter(k == 500, scatac_min_cell_filter == 1, scatac_percent_filter == FALSE) %>% 
    mutate(logFC = gsub("^.*_", "", expr) %>% 
               replace(. == 'all', '0.0')) %>%
    filter(!logFC %in% c('95%', '99%'))
enr = avg1 %>% 
    group_by_at(vars(-expr_set, -aucc)) %>% 
    summarise(delta = ifelse(n() == 2, 
                             aucc[expr_set == 'logFC'] - aucc[expr_set == 'Random'],
                             aucc),
              nn = n()) %>% 
    ungroup() 

grid = distinct(enr, method, logFC) %>% filter(logFC != '0.0')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    df = filter(avg1, method == current$method, logFC == current$logFC) %>% 
        mutate(expr_set = factor(expr_set, levels = c('logFC', 'Random')))
    eff = cohen.d(formula = aucc ~ expr_set | Subject(paste(paper, cell_type)),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})
labs = eff %>% 
    group_by(logFC) %>% 
    summarise(mean = mean(d), n = n(), median = median(d), d = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))

p2 = eff %>% 
    ggplot(aes(x = factor(logFC), y = d)) +
    geom_hline(aes(yintercept = 0), linetype = 'dotted', size = 0.25) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.6, size = 0.35,
                 fill = 'grey82', color = 'grey20') + 
    geom_label(data = labs, aes(y = Inf, label = label), size = 1.5, hjust = 0.5,
               vjust = 1, color = 'black', fill = NA, label.size = NA,
               label.padding = unit(0.6, 'lines')) +
    scale_x_discrete('% of genes removed') +
    scale_y_continuous('Cohen\'s d (vs. random gene removal)', limits=c(-0.8, 0.05)) +
    ggtitle('Effect of log-fold change filtering on concordance\nbetween ATAC and RNA modalities of multiome data') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = 1,
          plot.title = element_text(size=5, hjust=0.5)
          )
p2
ggsave("fig/Fig5/multiome-concordance-logFC-d.pdf", p2,
       width = 6, height = 6, units = "cm", useDingbats = FALSE)


###############################################################################-
## luecken logFC plot ####
###############################################################################-
dat = readRDS("data/summaries/false_discoveries/false_discoveries_logFC_filter.rds") %>% 
    type_convert() %>% 
    # flag study
    mutate(study = gsub("-.*$", "", filename)) %>% 
    filter(!grepl("latent_vars=t", filename)) %>% 
    dplyr::select(-filename) %>% 
    # filter out methods we are not using
    filter(da_family != 'mixedmodel',
           pseudo_repl == FALSE)
# restructure 'expr' column
dat %<>% 
    mutate(expr = replace(percentile, is.na(percentile), '0.0')) %>% 
    mutate(expr_set = ifelse(grepl('pseudo', expr), 'Random', 'logFC'),
           expr = gsub("_pseudo", "", expr))

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-', da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', da_family) %>% 
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

# re-code filtering
dat %<>%
    mutate(filter = ifelse(percent_filter, 
                           paste0(min_cell_filter, '%'),
                           paste0(min_cell_filter, ' cells')) %>% 
               fct_recode('1 cell' = '1 cells') %>% 
               fct_relevel('1 cell', '10 cells', '0.2%', '0.5%', '1%', '5%', '10%', '20%')) %>% 
    filter(!filter %in% c('10%', '20%'))

dat %<>% 
    mutate(logFC = gsub("^.*_", "", expr) %>% 
               replace(. == 'all', '0.0')) %>%
    filter(!logFC %in% c('95%', '99%'))
enr = dat %>% 
    group_by_at(vars(-study, -expr_set, -percentile,
                     -number_of_regions, -number_of_da_regions,
                     -number_of_nominal_da_regions)) %>% 
    summarise(delta = ifelse(n() == 2, 
                             number_of_da_regions[expr_set == 'logFC'] - number_of_da_regions[expr_set == 'Random'],
                             number_of_da_regions),
              nn = n()) %>% 
    ungroup() 

grid = distinct(enr, method, logFC) %>% filter(logFC != '0.0')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    df = filter(dat, method == current$method, logFC == current$logFC) %>% 
        mutate(expr_set = factor(expr_set, levels = c('logFC', 'Random')))
    eff = cohen.d(formula = number_of_da_regions ~ expr_set | Subject(label),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})
labs = eff %>% 
    group_by(logFC) %>% 
    summarise(mean = mean(d), n = n(), median = median(d), d = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))

p3 = eff %>% 
    ggplot(aes(x = factor(logFC), y = d)) +
    geom_hline(aes(yintercept = 0), linetype = 'dotted', size = 0.25) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.6, size = 0.35,
                 fill = 'grey82', color = 'grey20') + 
    geom_label(data = labs, aes(y = Inf, label = label), size = 1.5, hjust = 0.5,
               vjust = 1, color = 'black', fill = NA, label.size = NA,
               label.padding = unit(0.6, 'lines')) +
    # scale_color_manual('', values = pal) +
    # scale_fill_manual('', values = pal) +
    scale_x_discrete('% of peaks removed') +
    scale_y_continuous('Cohen\'s d (vs. random peak removal)') +
    ggtitle('Effect of log-fold change filtering on false discoveries\nin null comparisons of real scATAC-seq data') +
    coord_cartesian(ylim = c(0, 1.3)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = 1,
          plot.title = element_text(size=5)
          )
p3
ggsave("fig/Fig5/luecken-concordance-logFC-d.pdf", p3,
       width = 6, height = 6, units = "cm", useDingbats = FALSE)

p4 = p1 + theme(aspect.ratio = 0.9) | p2 + theme(aspect.ratio = 0.9) | p3 + theme(aspect.ratio = 0.9)
ggsave("fig/Fig5/row.pdf", p4,
       width = 16.5, height = 5.5, units = "cm", useDingbats = FALSE)
