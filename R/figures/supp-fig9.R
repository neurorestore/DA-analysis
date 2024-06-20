# Plot a summary figure with performance across all experiments.
setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(pals)
library(upstartr)
source('R/theme.R')

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
## AUCC: boxplot, overview ####
###############################################################################-

avg1 = avg %>% 
    filter(!grepl('pseudo-replicates', method)) %>% 
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    filter(k == 5000, filter == '1 cell') %>% 
    filter(!expr %in% c('logFC_filter_95%', 'logFC_filter_99%')) %>%
    mutate(logFC = gsub("^.*_", "", expr) %>% 
               replace(. == 'all', '0%'))
hlines = avg1 %>% filter(logFC == '0%') %>% 
    group_by(method, color) %>% 
    summarise(aucc = median(aucc))
dplyr::count(avg1, method)
labs = avg1 %>% 
    group_by(logFC, expr_set, method, color) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
pal = c('logFC' = 'black',
        'Random' = 'grey75') 
p1 = avg1 %>%
    ggplot(aes(x = logFC, y = aucc,
               color = expr_set, fill = expr_set)) +
    facet_wrap(~ method, scales = 'free', ncol = 5, labeller = label_parsed) +
    geom_hline(data = hlines, aes(yintercept = aucc), color = 'grey88',
               size = 0.4) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.7,
                 size = 0.3) + 
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    # geom_text(data = labs, aes(y = -0.12, label = label), size = 1.75, hjust = 0,
    #           color = 'black') +
    scale_color_manual('Filter by:', values = pal) +
    scale_fill_manual('Filter by:', values = pal) +
    scale_x_reordered('% of peaks removed') +
    scale_y_continuous('AUCC') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(legend.position = 'top',
          legend.key.width = unit(0.4, 'lines'),
          aspect.ratio = 0.7,
          axis.text.x = element_text(angle=45, hjust=1))
p1
ggsave("fig/Supp_Fig9/bulk-logFC-filter-overview.pdf", p1,
       width = 18, height = 15, units = "cm", useDingbats = FALSE)

###############################################################################-
## AUCC: boxplot, enrichment over random ####
###############################################################################-

enr = avg1 %>% 
    group_by_at(vars(-expr_set, -aucc)) %>% 
    summarise(delta = ifelse(n() == 2, 
                             aucc[expr_set == 'logFC'] - aucc[expr_set == 'Random'],
                             aucc),
              nn = n()) %>% 
    ungroup() 
pal = da_analysis_colors
p2 = enr %>%
    filter(logFC != '0%') %>% 
    ggplot(aes(x = logFC, y = delta,
               color = method, fill = method)) +
    facet_wrap(~ method, scales = 'free', nrow = 3, labeller = label_parsed) +
    geom_boxplot(aes(y = delta), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    # geom_text(data = labs, aes(y = -0.12, label = label), size = 1.75, hjust = 0,
    #           color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_reordered('% of peaks removed') +
    scale_y_continuous(expression(Delta~AUCC)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(legend.position = 'none',
          aspect.ratio = 0.8,
          axis.text.x = element_text(angle=45, hjust=1),
          axis.title.y = element_text(size=6)
          )
p2
ggsave("fig/Supp_Fig9/bulk-logFC-filter-enrichment.pdf", p2,
       width = 18, height = 15, units = "cm", useDingbats = FALSE)

