# Plot a summary figure with performance across all experiments.
setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(pals)
library(upstartr)
library(effsize)
source('R/theme.R')

# load data
dat = readRDS("data/summaries/concordance_res/matched_bulk/all_aucc_downsampled.rds") %>% 
    type_convert() %>% 
    filter(!sc_pseudo_repl)
# no binarization or normalization
table(dat$sc_binarization, dat$sc_normalization)
# simplify
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
distinct(dat[, !keep])
dat %<>% extract(, keep)

dat %<>%
    filter(downsample_proportion > 200)

# print all analyses
distinct(dat, paper, cell_type)

# print some summaries
## bulk DA methods
table(dat$bulk_da_type, dat$bulk_da_method)
## filtering
table(dat$sc_min_cell_filter, dat$sc_percent_filter)
## single-cell methods
dplyr::count(dat, sc_da_family, sc_da_method, sc_da_type, sc_pseudobulk)
## k
table(dat$k)
## downsample_proportion
table(dat$downsample_proportion)

# average over bulk methods
avg = dat %>% 
    group_by_at(vars(-bulk_da_method, -bulk_da_type,
                     -bulk_features, -overlap, -aucc)) %>% 
    summarise(aucc = mean(aucc),
              n = n()) %>% 
    ungroup()

# is the bulk grid complete?
table(avg$n)
# is the filtering grid complete?
table(avg$sc_min_cell_filter, avg$sc_percent_filter)
# is the comparison grid complete?
dplyr::count(avg, sc_da_family, sc_da_method, sc_da_type) %>% 
    as.data.frame()

# re-code x-values and colors
avg %<>%
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

# check we are testing the same number of features for each method
dplyr::count(avg, paper, cell_type, sc_min_cell_filter, sc_percent_filter, sc_features) %>% 
    arrange(sc_min_cell_filter, sc_percent_filter) %>% 
    as.data.frame()

###############################################################################-
## AUCC: boxplot ####
###############################################################################-

avg0 = filter(avg, k == 5000, sc_min_cell_filter == 1, !sc_percent_filter)
avg0$method = factor(avg0$method, levels=names(da_analysis_colors))
labs = avg0 %>% 
    group_by(method, color, downsample_proportion) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
lvls = labs %>% 
    filter(downsample_proportion == 1000) %>% 
    arrange(median) %>% 
    pull(method)
pal = da_analysis_colors
p1 = avg0 %>%
    ggplot(aes(x = factor(downsample_proportion), y = aucc,
               color = method, fill = method)) +
    facet_wrap(~ method, scales = 'free_y', nrow = 3, labeller = label_parsed) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.35) + 
    geom_text(data = labs, aes(y = -0.05, label = label), size = 1.75, hjust = 0.5,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_reordered('Read depth') +
    geom_hline(yintercept = 0, linetype = 'dashed', size=0.2, color='grey80') +
    scale_y_continuous(name='AUCC', limits=c(-0.1, 0.5), breaks=seq(0, 1, 0.2)) +
    # coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(
        # axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 0.8,
        legend.position = 'none'
        )
p1
ggsave("fig/Supp_Fig27/downsample-reads-aucc.pdf", p1,
       width = 15, height = 15, units = "cm", useDingbats = FALSE)



# load data
dat = readRDS("data/summaries/concordance_res/multiome/all_aucc_downsampled.rds") %>% 
    type_convert() 
# no binarization or normalization
table(dat$scatac_binarization, dat$scatac_normalization)
# simplify
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
distinct(dat[, !keep])
dat %<>% extract(, keep)

# print comparisons
distinct(dat, cell_type)
# remove mixed models and pseudoreplicates
dat %<>% filter(!scatac_pseudo_repl, scatac_da_family != 'mixedmodel')

# print some summaries
## scRNA DE methods
dplyr::count(dat, scrna_da_method, bulk_da_type)
## filtering grid
with(dat, table(scatac_min_cell_filter, scatac_percent_filter))
## single-cell DA methods
dplyr::count(dat, scatac_da_family, scatac_da_method, scatac_da_type, scatac_pseudobulk, scatac_pseudo_repl)
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
# is the downsample grid complete?
table(avg$downsample_proportion)
# is the comparison grid complete?
dplyr::count(avg, scatac_da_family, scatac_pseudo_repl,
             scatac_da_method, scatac_da_type) %>% 
    as.data.frame()

# re-code x-values and colors
avg %<>%
    mutate(method = paste0(scatac_da_method, '-', scatac_da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = ifelse(scatac_pseudo_repl,
                           paste(method, '(pseudo-replicates)'), method),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(scatac_pseudo_repl, 'pseudo-replicates', scatac_da_family),
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

###############################################################################-
## AUCC: boxplot, main-text ####
###############################################################################-

avg0 = filter(avg, k == 500, scatac_min_cell_filter == 1, !scatac_percent_filter)
avg0$method = factor(avg0$method, levels=names(da_analysis_colors))
labs = avg0 %>% 
    group_by(method, color, downsample_proportion) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
lvls = labs %>% 
    filter(downsample_proportion == 1000) %>% 
    arrange(median) %>% 
    pull(method)
pal = da_analysis_colors
p2 = avg0 %>%
    ggplot(aes(x = factor(downsample_proportion), y = aucc,
               color = method, fill = method)) +
    facet_wrap(~ method, scales = 'free_y', nrow = 4, labeller = label_parsed) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.35) + 
    geom_text(data = labs, aes(y = -0.05, label = label), size = 1.75, hjust = 0.5,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_reordered('# of cells') +
    geom_hline(yintercept = 0, linetype = 'dashed', size=0.2, color='grey80') +
    scale_y_continuous(name='AUCC', limits=c(-0.1, 0.5), breaks=seq(0, 1, 0.2)) +
    # coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(
        # axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 0.62,
        legend.position = 'none'
    )
p2
ggsave("fig/Supp_Fig27/downsample-cells-aucc.pdf", p2,
       width = 15, height = 15, units = "cm", useDingbats = FALSE)

