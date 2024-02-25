# Plot a summary figure with performance across all experiments.
setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(pals)
library(upstartr)
source('R/theme.R')

# load data
dat = readRDS("data/summaries/concordance_res/multiome/all_aucc_logFC_filter.rds") %>% 
    type_convert() 
# no binarization or normalization
table(dat$scatac_binarization, dat$scatac_normalization)
# simplify
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

# restructure 'expr' column
dat %<>% mutate(expr_set = ifelse(grepl('pseudo', expr), 'Random', 'logFC'),
                expr = gsub("_pseudo", "", expr))
dplyr::count(dat, expr, expr_set)

# print comparisons
distinct(dat, paper, cell_type)

# filter Zhu, Anandon, atlas by label
dat %<>%
    filter(
        !grepl('Zhu', paper),
        !grepl('Anadon', paper),
        !(paper == 'Squair_2022' & grepl('7d_vs_Uninjured', cell_type)),
        !(paper == 'Squair_2022' & grepl('2m_vs_Uninjured', cell_type)),
    )
distinct(dat, paper, cell_type)
distinct(dat, paper, cell_type) %>% dplyr::count(paper)

# remove mixed models and pseudoreplicates
dat %<>% filter(scatac_da_family != 'mixedmodel')

## scRNA DE methods
dplyr::count(dat, scrna_da_method, bulk_da_type)
## filtering grid
with(dat, table(scatac_min_cell_filter, scatac_percent_filter))
## single-cell DA methods
dplyr::count(dat, scatac_da_family, scatac_da_method, scatac_da_type, scatac_pseudobulk)
## k
table(dat$k)

# average over bulk methods
avg = dat %>% 
    group_by_at(vars(-scrna_da_method, -bulk_da_type, -sc_features,
                     -bulk_features, -overlap, -aucc)) %>% 
    summarise(aucc = mean(aucc),
              n = n()) %>% 
    ungroup()

# is the bulk grid complete?
table(avg$n)
# is the filtering grid complete?
table(avg$scatac_min_cell_filter, avg$scatac_percent_filter)
# is the comparison grid complete?
dplyr::count(avg, scatac_da_family, scatac_da_method, scatac_da_type) %>% 
    as.data.frame()

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
# print 
dplyr::count(avg, color, method) %>% arrange(color, method)

# re-code filtering
avg %<>%
    mutate(filter = ifelse(scatac_percent_filter, 
                           paste0(scatac_min_cell_filter, '%'),
                           paste0(scatac_min_cell_filter, ' cells')) %>% 
               fct_recode('1 cell' = '1 cells') %>% 
               fct_relevel('1 cell', '10 cells', '0.2%', '0.5%', '1%', '5%', 
                           '10%', '20%')) %>% 
    filter(!filter %in% c('10%', '20%'))

# check we are testing the same nubmer of features for each method
# dplyr::count(avg, paper, cell_type, scatac_min_cell_filter, 
#              scatac_percent_filter, sc_features) %>% 
#   arrange(scatac_min_cell_filter, scatac_percent_filter) %>% 
#   as.data.frame() %>% 
#   arrange(sc_features)
# hist(avg$sc_features)

###############################################################################-
## AUCC: boxplot, main-text ####
###############################################################################-

avg1 = avg %>% 
    filter(k == 500, scatac_min_cell_filter == 1, scatac_percent_filter == FALSE) %>% 
    filter(!expr %in% c('logFC_filter_95%', 'logFC_filter_99%')) %>%
    mutate(logFC = gsub("^.*_", "", expr) %>% 
               replace(. == 'all', '0%'))
hlines = avg1 %>% filter(logFC == '0.0') %>% 
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
    scale_x_reordered() +
    scale_y_continuous('AUCC') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(legend.position = 'top',
          legend.key.width = unit(0.4, 'lines'),
          aspect.ratio = 0.7,
          axis.text.x = element_text(angle=45, hjust=1))
p1
ggsave("fig/Supp_Fig10/multiome-logFC-filter-overview.pdf", p1,
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
    facet_wrap(~ method, scales = 'free', ncol = 5, labeller = label_parsed) +
    geom_boxplot(aes(y = delta), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    # geom_text(data = labs, aes(y = -0.12, label = label), size = 1.75, hjust = 0,
    #           color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_reordered() +
    scale_y_continuous(expression(Delta~AUCC)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(legend.position = 'none',
          aspect.ratio = 0.8,
          axis.text.x = element_text(angle=45, hjust=1),
          axis.title.y = element_text(size=6)
          )
p2
ggsave("fig/Supp_Fig10/multiome-logFC-filter-enrichment.pdf", p2,
       width = 18, height = 15, units = "cm", useDingbats = FALSE)

