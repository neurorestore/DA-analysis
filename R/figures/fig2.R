# Plot a summary figure with performance across all experiments.
setwd("C:/Users/teo/Documents/EPFL/projects/DA-analysis/")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(pals)
source("R/theme.R")

###############################################################################-
## Bulk concordance plot ####
###############################################################################-

dat = readRDS("data/summaries/concordance_res/matched_bulk/all_aucc_skinny.rds") %>% 
    type_convert() %>% 
    # here, we are not interested in terciles
    filter(expr == 'all')

# summarise across all bulk DA methods
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

# remove all unncessary stuff
avg %<>%
    filter(!grepl('pseudo-replicates', method)) %>%  
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    filter(k == 5000, filter == '1 cell')

# check we are testing the same nubmer of features for each method
dplyr::count(avg, paper, cell_type, sc_min_cell_filter, sc_percent_filter, sc_features) %>% 
    arrange(sc_min_cell_filter, sc_percent_filter) %>% 
    as.data.frame()

labs = avg %>% 
    group_by(k, filter, method, color) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
saveRDS(labs, 'data/summaries/meta_summaries/bulk_aucc_summary.rds')

pal = c('mixed model' = 'grey80',
        'single-cell' = colours.cafe322[1],
        'pseudobulk' = colours.cafe322[2],
        'pseudo-replicates' = "#6E645F",
        'other' = colours.cafe322[3])
pal = da_analysis_colors

p1 = avg %>%
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
    ggtitle('Concordance between scATAC-seq and\nmatched bulk ATAC-seq') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          plot.title = element_text(size=5,hjust = 0.5),
          aspect.ratio = 1.4)
p1
ggsave("fig/Fig2/bulk-concordance.pdf", p1, width = 5, height = 6.5, units = "cm", useDingbats = FALSE)

###############################################################################-
## Multiome concordance plot ####
###############################################################################-
# load data
dat = readRDS("data/summaries/concordance_res/multiome/all_aucc_skinny.rds") %>% 
    type_convert() %>% 
    filter(paper != 'Flochlay_2022') %>%
    dplyr::rename(scrna_da_type = bulk_da_type) %>%
    filter(
        !grepl('Zhu', paper),
        !grepl('Anadon', paper),
        !(paper == 'Squair_2022' & grepl('7d_vs_Uninjured', cell_type)),
        !(paper == 'Squair_2022' & grepl('2m_vs_Uninjured', cell_type)),
    )

# average over bulk methods
avg = dat %>% 
    group_by_at(vars(-scrna_da_family, -scrna_da_method, -scrna_da_type,
                     -bulk_features, -overlap, -aucc)) %>% 
    summarise(aucc = mean(aucc),
              n = n()) %>% 
    ungroup()

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

pal = da_analysis_colors
avg %<>% 
    filter(k == 500, scatac_min_cell_filter == 1, !scatac_percent_filter, !scatac_pseudo_repl)
labs = avg %>% 
    group_by(method, color) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
saveRDS(labs, 'data/summaries/meta_summaries/multiome_aucc_summary.rds')

p2 = avg %>%
    ggplot(aes(x = reorder(method, aucc, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.35) + 
    geom_text(data = labs, aes(y = -0.12, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('AUCC') +
    coord_flip() +
    ggtitle('Concordance between ATAC and RNA\nmodalities of multiome data, gene level') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none',
          plot.title = element_text(size=5, hjust=0.5))
p2
ggsave("fig/Fig2/multiome-concordance.pdf", p2,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)


###############################################################################-
## Multiome GO concordance plot ####
###############################################################################-

dat = readRDS("data/summaries/concordance_res/multiome/all_aucc_gsea.rds") %>% 
    type_convert() %>% 
    filter(paper != 'Flochlay_2022')

# filter Zhu, Anandon, atlas by label
dat %<>%
    filter(
        !grepl('Zhu', paper),
        !grepl('Anadon', paper),
        !(paper == 'Squair_2022' & grepl('7d_vs_Uninjured', cell_type)),
        !(paper == 'Squair_2022' & grepl('2m_vs_Uninjured', cell_type)),
    )
# remove mixed models and pseudoreplicates
dat %<>% filter(!scatac_pseudo_repl, scatac_da_family != 'mixedmodel')

# average over bulk methods
avg = dat %>% 
    group_by_at(vars(-scrna_da_family, -scrna_da_method, -scrna_da_type,
                     -bulk_features, -overlap, -aucc)) %>% 
    summarise(aucc = mean(aucc),
              n = n()) %>% 
    ungroup()

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


avg %<>% filter(k == 100, scatac_min_cell_filter == 1, !scatac_percent_filter)
labs = avg %>% 
    group_by(method, color) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
saveRDS(labs, 'data/summaries/meta_summaries/multiome_gsea_aucc_summary.rds')

pal = da_analysis_colors
p3 = avg %>%
    ggplot(aes(x = reorder(method, aucc, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -0.12, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('AUCC') +
    coord_flip() +
    ggtitle('Concordance between ATAC and RNA\nmodalities of multiome data, GO term level') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none',
          plot.title = element_text(size=5, hjust=0.5))
p3
ggsave("fig/Fig2/multiome-gsea-concordance.pdf", p3,
       width = 5, height = 6.5, units = "cm", useDingbats = FALSE)

p4 = p2 | p3
ggsave("fig/Fig2/combined-concordance.pdf", p4,
       width = 10.5, height = 7, units = "cm", useDingbats = FALSE)
