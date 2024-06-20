# Plot a summary figure with performance across all experiments.
setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(pals)
library(upstartr)
source('R/theme.R')

# load data
dat = readRDS("data/summaries/concordance_res/multiome/all_aucc_non_overlapping_gsea.rds") %>% 
    type_convert() %>% 
    filter(paper != 'Flochlay_2022')
# no percentiles
# no binarization or normalization
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
# unique columns
dat %<>% extract(, map_lgl(., ~ n_distinct(.x) > 1))

# average over bulk methods
avg = dat %>% 
    group_by_at(vars(-scrna_da_method, -scrna_da_type,
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

###############################################################################-
## AUCC: boxplot, main-text ####
###############################################################################-

avg0 = filter(avg, k == 100, scatac_min_cell_filter == 1, !scatac_percent_filter)
labs = avg0 %>% 
    group_by(method, color) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
saveRDS(labs, 'data/summaries/meta_summaries/multiome_gsea_aucc-bidirectional_promoter-summary.rds')

pal = da_analysis_colors
p1 = avg0 %>%
    mutate(title = 'Removal of overlapping\npromoter regions') %>%
    ggplot(aes(x = reorder(method, aucc, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -0.2, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('AUCC', limits = c(-0.2, 0.6), breaks = seq(0, 0.6, 0.2)) +
    coord_flip() +
    # ggtitle('Removal of overlapping\npromoters regions') +
    facet_wrap(~title) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none',
          plot.title = element_text(size=5)
          )
p1
ggsave("fig/Supp_Fig3/multiome-gsea-aucc-bidirectional-promoters.pdf", p1,
       width = 4.5, height = 6.5, units = "cm", useDingbats = FALSE)

# load data
dat = readRDS("data/summaries/concordance_res/multiome/all_aucc_gsea.rds") %>% 
    type_convert() %>% 
    filter(paper != 'Flochlay_2022')
# no percentiles
# no binarization or normalization
table(dat$scatac_binarization, dat$scatac_normalization)

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

# re-code filtering
avg %<>%
    mutate(filter = ifelse(scatac_percent_filter, 
                           paste0(scatac_min_cell_filter, '%'),
                           paste0(scatac_min_cell_filter, ' cells')) %>% 
               fct_recode('1 cell' = '1 cells') %>% 
               fct_relevel('1 cell', '10 cells', '0.2%', '0.5%', '1%', '5%', 
                           '10%', '20%')) %>% 
    filter(!filter %in% c('10%', '20%'))

grid = distinct(dat, scrna_da_family, scrna_da_method, scrna_da_type)
full_labs = data.frame()
plots = list()

for (i in 1:nrow(grid)) {
    current = grid[i,]
    # average over bulk methods
    avg = dat %>% 
        anti_join(current) %>% 
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
    
    # highlight held-out method
    held_out = with(current, paste0(scrna_da_method, '-', scrna_da_type))
    avg %<>% mutate(held_out = held_out)
    
    # plot
    avg0 = filter(avg, k == 100, scatac_min_cell_filter == 1, !scatac_percent_filter)
    labs = avg0 %>% 
        group_by(method, color) %>%
        summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
        ungroup() %>%
        mutate(label = formatC(median, format = 'f', digits = 2))
    full_labs %<>%
        rbind(labs %>% tidyr::crossing(current[2:3]))
    
    
    pal = da_analysis_colors
    p = avg0 %>%
        mutate(title = held_out) %>% 
        ggplot(aes(x = reorder(method, aucc, stats::median), 
                   color = method, fill = method)) +
        facet_grid(~ title) +
        geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                     size = 0.35) + 
        # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
        # stroke = 0.3) +
        # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
        geom_text(data = labs, aes(y = -0.2, label = label), size = 1.75, hjust = 0,
                  color = 'black') +
        scale_color_manual('', values = pal) +
        scale_fill_manual('', values = pal) +
        scale_x_discrete(labels = scales::parse_format()) +
        scale_y_continuous('AUCC', limits = c(-0.2, 0.6), breaks = seq(0, 0.6, 0.2)) +
        coord_flip() +
        boxed_theme(size_sm = 5, size_lg = 6) +
        theme(axis.title.y = element_blank(),
              plot.title = element_text(size = 6),
              aspect.ratio = 1.4,
              legend.position = 'none')
    plots[[length(plots)+1]] = p
}
# plot
p2 = wrap_plots(plots, nrow = 2, ncol=3)
saveRDS(full_labs, 'data/summaries/meta_summaries/multiome_gsea_aucc-leave_bulk_out-summary.rds')
p2
ggsave("fig/Supp_Fig3/multiome-gsea-aucc-leave-bulk-out.pdf", p2, width = 13.85, height = 15, units = "cm", useDingbats = FALSE)

# correlation matrix
bulk_grid = distinct(dat, scrna_da_method, scrna_da_type)
leave_out = pmap_dfr(bulk_grid, function(...) { 
    bulk_row = tibble(...)
    print(bulk_row)
    avg0 = dat %>% 
        anti_join(bulk_row) %>% 
        group_by_at(vars(-scrna_da_family, -scrna_da_method, -scrna_da_type,
                         -bulk_features, -overlap, -aucc)) %>% 
        summarise(aucc = mean(aucc),
                  n = n()) %>% 
        ungroup()
    
    # re-code methods
    avg0 %<>%
        mutate(method = paste0(scatac_da_method, '-', scatac_da_type) %>% 
                   gsub("-singlecell|-binomial|-exact", "", .),
               method = ifelse(scatac_pseudo_repl,
                               paste(method, '(pseudo-replicates)'), method),
               method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
                   gsub("fisher", "FET", .) %>% 
                   gsub("-LR_.*|-perm.*$", "", .),
               color = ifelse(scatac_pseudo_repl, 'pseudo-replicates', scatac_da_family),
               color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                              'singlecell', color) %>% 
                   fct_recode('single-cell' = 'singlecell',
                              'other' = 'non_libra') %>% 
                   fct_relevel('single-cell', 'pseudobulk', 'other'))
    
    # re-code filtering
    avg0 %<>%
        mutate(filter = ifelse(scatac_percent_filter, 
                               paste0(scatac_min_cell_filter, '%'),
                               paste0(scatac_min_cell_filter, ' cells')) %>% 
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
    
    # main-text params: k=5000, no filter
    avg0 %<>%
        filter(!grepl('pseudo-replicates', method)) %>%  
        filter(k == 100, filter == '1%')
    
    # flag left-out method
    title = with(bulk_row, paste0(scrna_da_method, '-', scrna_da_type))
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
p3 = ggplot(cor, aes(x = method1, y = method2, fill = estimate)) +
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
          legend.key.width = unit(0.14, 'lines'),
          legend.key.height = unit(0.14, 'lines'),
          legend.position = 'right',
          legend.justification = 'bottom')
p3
ggsave("fig/Supp_Fig3/multiome-gsea-aucc-leave-bulk-out-correlation.pdf",
       p3, width = 5, height = 5, units = "cm", useDingbats = FALSE)

# load data
dat = readRDS("data/summaries/concordance_res/multiome/all_aucc_gsea.rds") %>% 
    type_convert() %>% 
    filter(paper != 'Flochlay_2022')
# no percentiles
# no binarization or normalization
table(dat$scatac_binarization, dat$scatac_normalization)

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

# re-code filtering
avg %<>%
    mutate(filter = ifelse(scatac_percent_filter, 
                           paste0(scatac_min_cell_filter, '%'),
                           paste0(scatac_min_cell_filter, ' cells')) %>% 
               fct_recode('1 cell' = '1 cells') %>% 
               fct_relevel('1 cell', '10 cells', '0.2%', '0.5%', '1%', '5%', 
                           '10%', '20%')) %>% 
    filter(!filter %in% c('10%', '20%'))

avg0 = filter(avg, k == 50, scatac_min_cell_filter == 1, !scatac_percent_filter)
labs = avg0 %>% 
    group_by(method, color) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
saveRDS(labs, 'data/summaries/meta_summaries/multiome_gsea_aucc-k=50-summary.rds')

pal = da_analysis_colors
p4 = avg0 %>%
    mutate(title = paste('k = 50')) %>% 
    ggplot(aes(x = reorder(method, aucc, stats::median), 
               color = method, fill = method)) +
    facet_grid(~ title) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -0.2, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('AUCC', limits = c(-0.2, 0.6), breaks = seq(0, 0.6, 0.2)) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none'
          )
p4
ggsave("fig/Supp_Fig3/multiome-gsea-aucc-k=50.pdf", p4,
       width = 4.5, height = 6.5, units = "cm", useDingbats = FALSE)

avg0 = filter(avg, k == 500, scatac_min_cell_filter == 1, !scatac_percent_filter)
labs = avg0 %>% 
    group_by(method, color) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
saveRDS(labs, 'data/summaries/meta_summaries/multiome_gsea_aucc-k=500-summary.rds')

pal = da_analysis_colors
p5 = avg0 %>%
    mutate(title = paste('k = 500')) %>% 
    ggplot(aes(x = reorder(method, aucc, stats::median), 
               color = method, fill = method)) +
    facet_grid(~ title) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -0.2, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('AUCC', limits = c(-0.2, 0.6), breaks = seq(0, 0.6, 0.2)) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none'
          )
p5
ggsave("fig/Supp_Fig3/multiome-gsea-aucc-k=500.pdf", p4,
       width = 4.5, height = 6.5, units = "cm", useDingbats = FALSE)

avg0 = filter(avg, k == 100, scatac_min_cell_filter == 1, scatac_percent_filter)
labs = avg0 %>% 
    group_by(method, color) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
saveRDS(labs, 'data/summaries/meta_summaries/multiome_gsea_aucc-filter=1-summary.rds')

p6 = avg0 %>%
    mutate(title = paste('1% cell filter')) %>% 
    ggplot(aes(x = reorder(method, aucc, stats::median), 
               color = method, fill = method)) +
    facet_grid(~ title) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -0.2, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('AUCC', limits = c(-0.2, 0.6), breaks = seq(0, 0.6, 0.2)) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none'
          )
p6
ggsave("fig/Supp_Fig3/multiome-gsea-aucc-1pct.pdf", p5,
       width = 4.5, height = 6.5, units = "cm", useDingbats = FALSE)

# load data
dat = readRDS("data/summaries/concordance_res/multiome/all_aucc_scrna_terciles_gsea.rds") %>% 
    type_convert() 
# no binarization or normalization
table(dat$scatac_binarization, dat$scatac_normalization)
# simplify
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

# filter to terciles 2/3 only
dat %<>% filter(expr == 'mean_percentile_2_3')

# filter Zhu, Anandon, atlas by label
dat %<>%
    filter(
        !grepl('Zhu', paper),
        !grepl('Anadon', paper),
        !grepl('Flochlay', paper),
        !(paper == 'Squair_2022' & grepl('7d_vs_Uninjured', cell_type)),
        !(paper == 'Squair_2022' & grepl('2m_vs_Uninjured', cell_type)),
    )
# remove mixed models and pseudoreplicates
dat %<>% filter(!scatac_pseudo_repl, scatac_da_family != 'mixedmodel')

# average over bulk methods
avg = dat %>% 
    group_by_at(vars(-scrna_da_method, -scrna_da_type, -sc_features,
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

avg$scatac_min_cell_filter = 1

# re-code filtering
avg %<>%
    mutate(filter = ifelse(scatac_percent_filter, 
                           paste0(scatac_min_cell_filter, '%'),
                           paste0(scatac_min_cell_filter, ' cells')) %>% 
               fct_recode('1 cell' = '1 cells') %>% 
               fct_relevel('1 cell', '10 cells', '0.2%', '0.5%', '1%', '5%', 
                           '10%', '20%')) %>% 
    filter(!filter %in% c('10%', '20%'))

avg0 = filter(avg, k == 100, scatac_min_cell_filter == 1, !scatac_percent_filter)
labs = avg0 %>% 
    group_by(method, color, expr) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
saveRDS(labs, 'data/summaries/meta_summaries/multiome_gsea_aucc-scrna_terciles-summary.rds')

lvls = labs %>% 
    filter(expr == 'all') %>% 
    arrange(median) %>% 
    pull(method)
pal = da_analysis_colors
p7 = avg0 %>%
    mutate(title = 'Removal of\nlowly-expressed genes') %>%
    ggplot(aes(x = reorder_within(method, aucc, expr, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.35) + 
    facet_wrap(~title) +
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -0.2, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    # scale_x_discrete(labels = scales::parse_format()) +
    scale_x_reordered() +
    scale_y_continuous('AUCC', limits = c(-0.2, 0.6), breaks = seq(0, 0.6, 0.2)) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none',
          plot.title = element_text(size=5)
          )
p7
ggsave("fig/Supp_Fig3/multiome-gsea-aucc-scRNA-terciles.pdf", p7,
       width = 4.5, height = 6.5, units = "cm", useDingbats = FALSE)

bottom_row = p4 | p5 | p6 | p7
bottom_row
ggsave("fig/Supp_Fig3/bottom-row.pdf", bottom_row,
       width = 18.35, height = 6.5, units = "cm", useDingbats = FALSE)

