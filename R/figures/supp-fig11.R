# Plot a summary figure with performance across all experiments.
setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(pals)
library(upstartr)
source('R/theme.R')

# load data
dat = readRDS("data/summaries/false_discoveries/false_discoveries_logFC_filter.rds") %>% 
    type_convert() %>% 
    # flag study
    mutate(study = gsub("-.*$", "", basename(filename))) %>% 
    filter(!grepl("latent_vars=t", filename)) %>% 
    # dplyr::select(-filename) %>% 
    # filter out methods we are not using
    filter(da_family != 'mixedmodel',
           pseudo_repl == FALSE)

# restructure 'expr' column
dat %<>% 
    mutate(expr = replace(percentile, is.na(percentile), '0%')) %>% 
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

###############################################################################-
## AUCC: boxplot, overview ####
###############################################################################-
dat %<>% 
    mutate(logFC = gsub("^.*_", "", expr)) %>%
    filter(!logFC %in% c('95%', '99%'))
hlines = dat %>% filter(logFC == '0%') %>% 
    group_by(method, color) %>% 
    summarise(number_of_da_regions = median(number_of_da_regions))
dplyr::count(dat, method)
labs = dat %>% 
    group_by(logFC, expr_set, method, color) %>%
    summarise(mean = mean(number_of_da_regions),
              n = n(), 
              median = median(number_of_da_regions), 
              number_of_da_regions = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
pal = c('logFC' = 'black',
        'Random' = 'grey75') 
p1 = dat %>%
    ggplot(aes(x = logFC, y = number_of_da_regions,
               color = expr_set, fill = expr_set)) +
    facet_wrap(~ method, scales = 'free', ncol = 5, labeller = label_parsed) +
    geom_hline(data = hlines, aes(yintercept = number_of_da_regions),
               color = 'grey88', size = 0.4) +
    geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, alpha = 0.4,
                 width = 0.7, size = 0.3) + 
    scale_color_manual('Filter by:', values = pal) +
    scale_fill_manual('Filter by:', values = pal) +
    scale_x_reordered() +
    scale_y_continuous('# of DA peaks') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(legend.position = 'top',
          legend.key.width = unit(0.4, 'lines'),
          aspect.ratio = 0.7,
          axis.text.x = element_text(angle=45, hjust=1))
p1
ggsave("fig/Supp_Fig11/luecken-logFC-filter-overview.pdf", p1,
       width = 18, height = 15, units = "cm", useDingbats = FALSE)

###############################################################################-
## AUCC: boxplot, enrichment over random ####
###############################################################################-

enr = dat %>% 
    group_by_at(vars(-study, -expr_set, -percentile,
                     -number_of_regions, -number_of_da_regions,
                     -number_of_nominal_da_regions)) %>% 
    summarise(delta = ifelse(n() == 2, 
                             number_of_da_regions[expr_set == 'logFC'] - number_of_da_regions[expr_set == 'Random'],
                             number_of_da_regions),
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
    scale_y_continuous(expression(Delta~'# of DA peaks')) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(legend.position = 'none',
          aspect.ratio = 0.8,,
          axis.text.x = element_text(angle=45, hjust=1))
p2
ggsave("fig/Supp_Fig11/luecken-logFC-filter-enrichment.pdf", p2,
       width = 18, height = 15, units = "cm", useDingbats = FALSE)

