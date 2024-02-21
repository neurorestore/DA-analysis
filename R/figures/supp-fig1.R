# Plot a summary figure with performance across all experiments.
setwd("C:/Users/teo/Documents/EPFL/projects/DA-analysis/")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(pals)
library(upstartr)
library(effsize)
source('R/theme.R')

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
                               'LR[peaks]' = 'LR_peaks',
                               'LR[clusters]' = 'LR',
                               'Binomial' = 'binomial',
                               'Wilcoxon~rank-sum~test' = 'wilcox',
                               'Fisher~exact~test' = 'FET',
                               'Negative~binomial' = 'negbinom') %>% 
               as.character()) %>% 
    # drop limma
    filter(!grepl('limma', method))

# check we are testing the same nubmer of features for each method
dplyr::count(avg, paper, cell_type, sc_min_cell_filter, sc_percent_filter, sc_features) %>% 
    arrange(sc_min_cell_filter, sc_percent_filter) %>% 
    as.data.frame()

# drop some columns to keep things clear
# avg %<>% 
#   dplyr::select(-binsize, -sc_min_cell_filter, -sc_percent_filter, 
#                 -sc_pseudobulk, -sc_pseudo_repl, -sc_da_family, -sc_da_method,
#                 -sc_da_type, -expr)

###############################################################################-
## AUCC: effect of k ####
###############################################################################-

avg3 = avg %>% 
    filter(!grepl('pseudo-replicates', method)) %>% 
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    filter(k == 1000, filter == '1 cell')
dplyr::count(avg3, method)
labs = avg3 %>% 
    group_by(k, filter, method, color) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
pal = da_analysis_colors
p1 = avg3 %>%
    ggplot(aes(x = reorder(method, aucc, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -0.12, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('AUCC', breaks = seq(0, 1, 0.1)) +
    coord_flip() +
    ggtitle('k=1000') + 
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          plot.title = element_text(size=5, hjust=0.5),
          aspect.ratio = 1.4)
# p1
ggsave("fig/Supp_Fig1/bulk-aucc-k=1000.pdf", p1,
       width = 4.5, height = 6.5, units = "cm", useDingbats = FALSE)

###############################################################################-
## AUCC: effect of filtering ####
###############################################################################-

avg4 = avg %>% 
    filter(!grepl('pseudo-replicates', method)) %>% 
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    filter(k == 5000, filter == '5%')
labs = avg4 %>% 
    group_by(k, filter, method, color) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
p2 = avg4 %>%
    ggplot(aes(x = reorder(method, aucc, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -0.12, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('AUCC', breaks = seq(0, 1, 0.2)) +
    coord_flip() +
    ggtitle('5% cell filter') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          aspect.ratio = 1.4,
          plot.title = element_text(size=5, hjust=0.5)
          )
p2
ggsave("fig/Supp_Fig1/bulk-aucc-filter=5pct.pdf", p2,
       width = 4.5, height = 6.5, units = "cm", useDingbats = FALSE)

###############################################################################-
## AUCC: leave-bulk-out ####
###############################################################################-

bulk_grid = distinct(dat, bulk_da_method, bulk_da_type)
plots = map(seq_len(nrow(bulk_grid)), ~ {
    bulk_idx = .x
    bulk_row = bulk_grid[bulk_idx, ]
    print(bulk_row)
    avg0 = dat %>% 
        anti_join(bulk_row) %>% 
        group_by_at(vars(-bulk_da_family, -bulk_da_method, -bulk_da_type,
                         -bulk_features, -overlap, -aucc)) %>% 
        summarise(aucc = mean(aucc),
                  n = n()) %>% 
        ungroup()
    
    # re-code methods
    avg0 %<>%
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
    avg0 %<>%
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
                                   'LR[peaks]' = 'LR_peaks',
                                   'LR[clusters]' = 'LR',
                                   'Binomial' = 'binomial',
                                   'Wilcoxon~rank-sum~test' = 'wilcox',
                                   'Fisher~exact~test' = 'FET',
                                   'Negative~binomial' = 'negbinom') %>% 
                   as.character()) %>% 
        # drop limma
        filter(!grepl('limma', method))
    
    # main-text params: k=5000, no filter
    avg0 %<>%
        filter(!grepl('pseudo-replicates', method)) %>% 
        filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type)))  %>% 
        filter(k == 5000, filter == '1 cell')
    
    ## label median
    labs = avg0 %>% 
        group_by(k, filter, method, color) %>%
        summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
        ungroup() %>%
        mutate(label = formatC(median, format = 'f', digits = 2))
    pal = da_analysis_colors
    
    # plot
    title = with(bulk_row, paste0(bulk_da_method, '-', bulk_da_type))
    p = avg0 %>%
        ggplot(aes(x = reorder(method, aucc, stats::median), 
                   color = method, fill = method)) +
        ggtitle(title) +
        geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                     size = 0.4) + 
        # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
        # stroke = 0.3) +
        # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
        geom_text(data = labs, aes(y = -0.12, label = label), size = 1.75, hjust = 0,
                  color = 'black') +
        scale_color_manual('', values = pal) +
        scale_fill_manual('', values = pal) +
        scale_x_discrete(labels = scales::parse_format()) +
        scale_y_continuous('AUCC', breaks = seq(0, 1, 0.1)) +
        coord_flip() +
        boxed_theme(size_sm = 5, size_lg = 6) +
        theme(axis.title.y = element_blank(),
              legend.position = 'none',
              aspect.ratio = 1.4,
              plot.title = element_text(size = 5, hjust = 0.5))
    p
})
p3 = wrap_plots(plots, nrow = 2, ncol=3)
p3
# ggsave("fig/Supp_Fig1/bulk-aucc-leave-bulk-out.pdf", p3,
#        width = 15.5, height = 13, units = "cm", useDingbats = FALSE)
ggsave("fig/Supp_Fig1/bulk-aucc-leave-bulk-out.pdf", p3,
       width = 13.85, height = 15, units = "cm", useDingbats = FALSE)

# correlation matrix
leave_out = pmap_dfr(bulk_grid, function(...) { 
    bulk_row = tibble(...)
    print(bulk_row)
    avg0 = dat %>% 
        anti_join(bulk_row) %>% 
        group_by_at(vars(-bulk_da_family, -bulk_da_method, -bulk_da_type,
                         -bulk_features, -overlap, -aucc)) %>% 
        summarise(aucc = mean(aucc),
                  n = n()) %>% 
        ungroup()
    
    # re-code methods
    avg0 %<>%
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
    avg0 %<>%
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
                                   'LR[peaks]' = 'LR_peaks',
                                   'LR[clusters]' = 'LR',
                                   'Binomial' = 'binomial',
                                   'Wilcoxon~rank-sum~test' = 'wilcox',
                                   'Fisher~exact~test' = 'FET',
                                   'Negative~binomial' = 'negbinom') %>% 
                   as.character()) %>% 
        # drop limma
        filter(!grepl('limma', method))
    
    # main-text params: k=5000, no filter
    avg0 %<>%
        filter(!grepl('pseudo-replicates', method)) %>% 
        filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type)))  %>% 
        filter(k == 5000, filter == '1 cell')
    
    # flag left-out method
    title = with(bulk_row, paste0(bulk_da_method, '-', bulk_da_type))
    avg0 %<>% 
        mutate(left_out = title)
    return(avg0)
})

# correlate
library(broom)
methods = unique(leave_out$left_out)
cor = tidyr::crossing(method1 = methods, 
                      method2 = methods) %>% 
    pmap_dfr(function(...) {
        current = tibble(...)
        x = filter(leave_out, left_out == current$method1) %>% 
            arrange(paper, cell_type, method) %>% 
            pull(aucc)
        y = filter(leave_out, left_out == current$method2) %>% 
            arrange(paper, cell_type, method) %>% 
            pull(aucc)
        tidy(cor.test(x, y)) %>% 
            cbind(current, .)
    })
range = range(cor$estimate)
brks = c(range[1] + 0.1 * diff(range),
         range[2] - 0.1 * diff(range))
labs = format(range, digits = 3)
p4 = ggplot(cor, aes(x = method1, y = method2, fill = estimate)) +
    geom_tile(color = 'white') + 
    coord_fixed() + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37", 
                           name = 'Correlation', 
                           limits = range, breaks = brks, labels = labs) +
    guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
    boxed_theme() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.key.width = unit(0.18, 'lines'),
          legend.key.height = unit(0.14, 'lines'),
          legend.position = 'right',
          legend.justification = 'bottom')
p4
ggsave("fig/Supp_Fig1/bulk-aucc-leave-bulk-out-correlation.pdf",
       p4, width = 5, height = 5, units = "cm", useDingbats = FALSE)

# load data
dat = readRDS("data/summaries/concordance_res/matched_bulk/all_aucc_pseudobulk_peaks.rds") %>% 
    type_convert() %>% 
    # here, we are not interested in terciles
    filter(expr == 'all')
keep = map_lgl(dat, ~ n_distinct(.) > 1)
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
    mutate(sc_min_cell_filter = 1) %>% 
    mutate(filter = ifelse(sc_percent_filter,
                           paste0(sc_min_cell_filter, '%'),
                           paste0(sc_min_cell_filter, ' cells')) %>%
               fct_recode('1 cell' = '1 cells') %>%
               fct_relevel('1 cell', '10 cells', '0.2%', '0.5%', '1%', '5%', '10%', '20%')) %>%
    # filter(!filter %in% c('10%', '20%')) %>% 
    mutate(method = fct_recode(method,
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

###############################################################################-
## AUCC: boxplot, main text ####
###############################################################################-

avg2 = avg %>% 
    filter(!grepl('pseudo-replicates', method)) %>% 
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    filter(# k == 5000, 
        filter == '1 cell')
dplyr::count(avg2, method)
labs = avg2 %>% 
    group_by(# k, 
        pseudobulk_peak_type,
        filter, method, color) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
pal = da_analysis_colors

p5 = avg2 %>%
    filter(pseudobulk_peak_type == 'combined') %>% 
    ggplot(aes(x = reorder(method, aucc, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs %>% filter(pseudobulk_peak_type == 'combined'), 
              aes(y = -0.12, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('AUCC', breaks = seq(0, 1, 0.2)) +
    coord_flip() +
    ggtitle('Pesudobulk peak calling') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          plot.title = element_text(size=5),
          aspect.ratio = 1.4)
p5
ggsave("fig/Supp_Fig1/bulk-aucc-ArchR-peak-calling.pdf", p5, width = 4.5, height = 6.5, units = "cm", useDingbats = FALSE)

###############################################################################-
## AUCC: test for interaction term ####
###############################################################################-

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
comb = bind_rows(avg1 %>% mutate(mode = 'bulk'),
                 avg2 %>% mutate(mode = 'pseudobulk') %>% 
                     filter(pseudobulk_peak_type == 'combined'))

library(lmerTest)
ctrl = lmerControl(check.conv.singular = .makeCC(action = "ignore", 
                                                 tol = 1e-4))
comb %<>% mutate(sample = paste(paper, cell_type))
fit = lmerTest::lmer(aucc ~ mode * method + (1 | sample), comb, 
                     REML = TRUE, control = ctrl)
coefs = coef(summary(fit)) %>% 
    as.data.frame() %>% 
    dplyr::rename(pval = `Pr(>|t|)`) %>% 
    mutate(padj = p.adjust(pval, 'BH'))
filter(coefs, padj < 0.05)
filter(coefs, pval < 0.05)

fit = lm(aucc ~ mode * method, comb)
coefs = coef(summary(fit)) %>% 
    as.data.frame() %>% 
    dplyr::rename(pval = `Pr(>|t|)`) %>% 
    mutate(padj = p.adjust(pval, 'BH'))
filter(coefs, padj < 0.05)
filter(coefs, pval < 0.05)


# load data
dat = readRDS("data/summaries/concordance_res/matched_bulk/all_aucc_noisy_peaks.rds") %>% 
    type_convert() %>% 
    # here, we are not interested in terciles
    filter(expr == 'all')
keep = map_lgl(dat, ~ n_distinct(.) > 1)
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
    mutate(sc_min_cell_filter = 1) %>% 
    mutate(filter = ifelse(sc_percent_filter,
                           paste0(sc_min_cell_filter, '%'),
                           paste0(sc_min_cell_filter, ' cells')) %>%
               fct_recode('1 cell' = '1 cells') %>%
               fct_relevel('1 cell', '10 cells', '0.2%', '0.5%', '1%', '5%', '10%', '20%')) %>%
    # filter(!filter %in% c('10%', '20%')) %>% 
    mutate(method = fct_recode(method,
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

avg2 = avg %>% 
    filter(!grepl('pseudo-replicates', method)) %>% 
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    filter(# k == 5000, 
        filter == '1 cell')

labs = avg2 %>% 
    group_by(filter, method, color) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
pal = da_analysis_colors
p6 = avg2 %>%
    ggplot(aes(x = reorder(method, aucc, stats::median), 
               color = method, fill = method)) +
    # facet_wrap(~ pseudobulk_peak_type, scales = 'free') +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -0.12, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    # scale_x_discrete(labels = scales::parse_format()) +
    scale_x_reordered() + 
    scale_y_continuous('AUCC', breaks = seq(0, 1, 0.1)) +
    ggtitle('Spurious peaks added') +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          plot.title = element_text(size=5),
          aspect.ratio = 1.4)
p6
ggsave("fig/Supp_Fig1/bulk-aucc-noisy-peaks.pdf", p1,
       width = 4.5, height = 6.5, units = "cm", useDingbats = FALSE)

# load original data
dat1 = readRDS("data/summaries/concordance_res/matched_bulk//all_aucc_skinny.rds") %>% 
    type_convert() %>% 
    filter(expr == 'all')

# load noisy data
dat2 = readRDS("data/summaries/concordance_res/matched_bulk/all_aucc_noisy_peaks.rds") %>% 
    type_convert() %>% 
    # here, we are not interested in terciles
    filter(expr == 'all')

# combine
dat = bind_rows(dat1 %>% mutate(mode = 'Original'),
                dat2 %>% mutate(mode = 'Noisy')) %>% 
    extract(, map_lgl(., ~ n_distinct(.) > 1))
# filter
dat %<>%
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    filter(k == 5000, 
           sc_pseudo_repl == FALSE,
           sc_percent_filter == FALSE,
           sc_min_cell_filter == 1) %>% 
    extract(, map_lgl(., ~ n_distinct(.) > 1))
# re-code x-values and colors
dat %<>%
    mutate(method = paste0(sc_da_method, '-', sc_da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .),
           color = sc_da_family,
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', color) %>% 
               fct_recode('single-cell' = 'singlecell',
                          'other' = 'non_libra') %>% 
               fct_relevel('single-cell', 'pseudobulk', 'other'))

###############################################################################-
## AUCC: test for interaction term ####
###############################################################################-

library(lmerTest)
ctrl = lmerControl(check.conv.singular = .makeCC(action = "ignore", 
                                                 tol = 1e-4))
dat %<>% mutate(sample = paste(paper, cell_type))
fit = lmerTest::lmer(aucc ~ mode * method + (1 | sample), dat, 
                     REML = TRUE, control = ctrl)
coefs = coef(summary(fit)) %>% 
    as.data.frame() %>% 
    dplyr::rename(pval = `Pr(>|t|)`) %>% 
    mutate(padj = p.adjust(pval, 'BH'))
filter(coefs, padj < 0.05)
filter(coefs, pval < 0.05)

fit = lm(aucc ~ mode * method , dat)
coefs = coef(summary(fit)) %>% 
    as.data.frame() %>% 
    dplyr::rename(pval = `Pr(>|t|)`) %>% 
    mutate(padj = p.adjust(pval, 'BH'))
filter(coefs, padj < 0.05)
filter(coefs, pval < 0.05)
# load data
dat = readRDS("data/summaries/concordance_res/matched_bulk/all_aucc_enhancer_promoter.rds") %>% 
    type_convert() %>% 
    # here, we are not interested in terciles
    filter(expr == 'all')
keep = map_lgl(dat, ~ n_distinct(.) > 1)
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
    mutate(sc_min_cell_filter = 1) %>% 
    mutate(filter = ifelse(sc_percent_filter,
                           paste0(sc_min_cell_filter, '%'),
                           paste0(sc_min_cell_filter, ' cells')) %>%
               fct_recode('1 cell' = '1 cells') %>%
               fct_relevel('1 cell', '10 cells', '0.2%', '0.5%', '1%', '5%', '10%', '20%')) %>%
    # filter(!filter %in% c('10%', '20%')) %>% 
    mutate(method = fct_recode(method,
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

###############################################################################-
## AUCC: boxplot, main text ####
###############################################################################-

avg2 = avg %>% 
    filter(!grepl('pseudo-replicates', method)) %>% 
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    filter(# k == 5000, 
        filter == '1 cell')
# inspect number of rows
filter(avg2, cell_type == 'CMP_vs_GMP', sc_da_method == 'negbinom')
dplyr::count(avg2, method)
labs = avg2 %>% 
    group_by(filter, method, color, ann) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
pal = da_analysis_colors
p7_1 = avg2 %>%
    filter(ann == 'enhancer') %>%
    ggplot(aes(x = reorder(method, aucc, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs %>% filter(ann == 'enhancer'), aes(y = -0.12, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('AUCC', breaks = seq(0, 1, 0.1)) +
    coord_flip() +
    ggtitle('Enhancer regions only') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          plot.title = element_text(size=5),
          aspect.ratio = 1.4)
p7_1

p7_2 = avg2 %>%
    filter(ann == 'promoter') %>%
    ggplot(aes(x = reorder(method, aucc, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs %>% filter(ann == 'promoter'), aes(y = 0.7, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('AUCC', limits = c(0.7, 0.81), breaks = seq(0.75, 1, 0.05)) +
    coord_flip() +
    ggtitle('Promoter regions only') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          plot.title = element_text(size=5),
          aspect.ratio = 1.4)
p7_2
p7 = p7_2 | p7_1
p7
ggsave("fig/Supp_Fig1/bulk-aucc-enhancer-promoter.pdf", p7,
       width = 9.3, height = 6.5, units = "cm", useDingbats = FALSE)

middle_row = p1 | p2 | p5 | p6
middle_row
ggsave("fig/Supp_Fig1/middle-row.pdf", middle_row,
       width = 18.35, height = 6.5, units = "cm", useDingbats = FALSE)
