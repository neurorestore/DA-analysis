# Plot a summary figure with performance across all experiments.
setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(pals)
library(upstartr)
source('R/theme.R')

dat = readRDS('data/summaries/false_discoveries/false_discoveries_simulated.rds') %>%
    type_convert() %>% 
    dplyr::select(-ori_filename) %>% 
    ## remove libra bug
    dplyr::select(-de_family, -de_method) %>% 
    dplyr::rename(da_type = de_type) %>%
    ## cell type is actually simulation replicate
    dplyr::rename(replicate = cell_type) %>% 
    # remove pseudoreplicates and mixed models
    filter(!pseudo_repl, da_family != 'mixedmodel') %>% 
    # filter to new grid
    filter(de_loc >= 0.5) %>%
    filter(!binarization, normalization=='f') %>%
    filter(is.na(percentile))

# simulation parameters
table(dat$n_cells, dat$n_reps)
table(dat$de_prob, dat$de_loc)
## DA methods
dplyr::count(dat, da_family, da_method, da_type, pseudobulk, pseudo_repl)
# replicate
table(dat$replicate)

## n should be 150 for everything
dplyr::count(dat, n_cells, n_reps, de_prob, de_loc)
filter(dat, n_cells == 5000) %>% 
    dplyr::count(da_family, da_method, da_type, pseudobulk, pseudo_repl)

# all 5% FDR# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-', da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = ifelse(pseudo_repl,
                           paste(method, '(pseudo-replicates)'), method),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(pseudo_repl, 'pseudo-replicates', da_family),
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

###############################################################################-
## Number of DA regions, boxplot, n_cells ####
###############################################################################-

## defaults
dat0 = filter(dat, de_loc == 1, de_prob == 0.5, n_reps == 6) %>% 
    filter(n_cells != 100) %>% 
    mutate(n_cells = n_cells / 2)
labs0 = dat0 %>%
    group_by(method, color, n_cells) %>%
    summarise(mean = mean(number_of_da_regions), 
              n = n(), 
              median = median(number_of_da_regions), 
              number_of_da_regions = median) %>%
    ungroup() %>%
    mutate(label = round(median))
pal = da_analysis_colors
p1 = dat0 %>%
    ggplot(aes(x = reorder(method, number_of_da_regions, stats::median), 
               color = method, fill = method)) +
    facet_wrap(~ n_cells, nrow = 1, # scales = 'free_x',
               labeller = as_labeller(~ paste(format(., big.mark = ',') %>% trimws(), 'cells'))) +
    geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, 
                 alpha = 0.4, width = 0.6, size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs0, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs0, aes(y = -9e3, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous(expression('# of DA peaks'~(10^3)), 
                       labels = function(x) x / 1e3,
                       breaks = seq(0, 80, 10) * 1e3
    ) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none'
          )
p1
ggsave("fig/Supp_Fig4/false-discoveries-n_cells.pdf", p1,
       width = 18, height = 6.5, units = "cm", useDingbats = FALSE)

## defaults
dat0 = filter(dat, de_loc == 1, de_prob == 0.5, n_cells == 1000) %>% 
    mutate(n_reps = n_reps / 2)
labs0 = dat0 %>%
    group_by(method, color, n_reps) %>%
    summarise(mean = mean(number_of_da_regions), 
              n = n(), 
              median = median(number_of_da_regions), 
              number_of_da_regions = median) %>%
    ungroup() %>%
    mutate(label = round(median))
pal = da_analysis_colors
p2 = dat0 %>%
    ggplot(aes(x = reorder(method, number_of_da_regions, stats::median), 
               color = method, fill = method)) +
    facet_wrap(~ n_reps, nrow = 1, # scales = 'free_x',
               labeller = as_labeller(~ paste(., 'vs.', .))) +
    geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, 
                 alpha = 0.4, width = 0.6, size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs0, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs0, aes(y = -6.5e3, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous(expression('# of DA peaks'~(10^3)), 
                       labels = function(x) x / 1e3,
                       # breaks = seq(0, 30, 10) * 1e3
    ) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none'
          )
p2
ggsave("fig/Supp_Fig4//false-discoveries-n_reps.pdf", p2,
       width = 18, height = 6.5, units = "cm", useDingbats = FALSE)

###############################################################################-
## Number of DA regions, boxplot, de_loc ####
###############################################################################-

## defaults
dat0 = filter(dat, de_prob == 0.5, n_cells == 1000, n_reps == 6) %>% 
    filter(de_loc %in% c(0.5, 1, 1.5, 2))
labs0 = dat0 %>%
    group_by(method, color, de_loc) %>%
    summarise(mean = mean(number_of_da_regions), 
              n = n(), 
              median = median(number_of_da_regions), 
              number_of_da_regions = median) %>%
    ungroup() %>%
    mutate(label = round(median))
pal = da_analysis_colors
p3 = dat0 %>%
    ggplot(aes(x = reorder(method, number_of_da_regions, stats::median), 
               color = method, fill = method)) +
    facet_wrap(~ de_loc, nrow = 1, # scales = 'free_x',
               labeller = as_labeller(~ paste0('de.loc=', .))) +
    geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, 
                 alpha = 0.4, width = 0.6, size = 0.35) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs0, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs0, aes(y = -11e3, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous(expression('# of DA peaks'~(10^3)), 
                       labels = function(x) x / 1e3,
                       breaks = seq(0, 80, 10) * 1e3
    ) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.4,
          legend.position = 'none'
          )
p3
ggsave("fig/Supp_Fig4/false-discoveries-de_loc.pdf", p3,
       width = 12, height = 6.5, units = "cm", useDingbats = FALSE)

data_file = 'data/final_datasets/simulated/peaks_with_sequencing_depth/n_cells=1000-n_reps=6-de_loc=0.5-de_prob=0.5-simulation_iter=1.rds'
sc = readRDS(data_file)

expr_df = data.frame(
    'Cell'=colnames(sc),
    'val' = colSums(GetAssayData(sc, slot='counts'))
) %>%
    left_join(sc@meta.data, by='Cell') %>%
    mutate(
        label = ifelse(label == 'unst', 'control', 'perturbed') 
    )

labs = expr_df %>%
    group_by(replicate) %>%
    summarize(
        # stats_val = median(val),
        stats_val = mean(val)/10000,
        val = 3.5
    ) %>%
    ungroup() %>%
    mutate(
        text_val = ifelse(stats_val < 0, 'NA', format(stats_val, digits = 2)),
        text_y = ifelse(val < 0, 0, val)
    )

med = function(x) stats::median(x) %>% replace(is.na(.), 0)
color_pal = c("#088980","#FE9E53") %>% setNames(c('control', 'perturbed'))
p4 = expr_df %>%
    mutate(val = val/10000) %>%
    ggplot(aes(x = reorder(replicate, val, med), y = val, color=label, fill=label)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
    geom_label(dat = labs,
               aes(label = text_val, y = text_y), color = 'black',
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0, vjust = 0.5,
               show.legend = FALSE) +
    coord_flip() +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    boxed_theme() +
    xlab('Library') +
    ylab(expression('# of fragments observed ('~10^4~')')) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 1.5,
        legend.title = element_blank(),
        legend.position = 'top',
        legend.key.width = unit(0.18, 'lines'),
        legend.key.height = unit(0.18, 'lines')
    ) + 
    ylim(0, 4.5)
p4
ggsave('fig/Supp_Fig4/splat_depth.pdf', p4, height=4.5, width=6.5, units='cm')

dat = readRDS("data/summaries/false_discoveries/false_discoveries_luecken.rds") %>% 
    type_convert() %>% 
    # flag study
    mutate(study = gsub("-.*$", "", filename)) %>% 
    dplyr::select(-filename) %>% 
    # filter out methods we are not using
    filter(pseudo_repl == FALSE) %>% 
    # ignore percentile here
    filter(is.na(percentile)) %>%
    filter(normalization == 'f', !binarization)
table(dat$study) ## only Luecken peak matrix

# count latent vars
dplyr::count(dat, latent_vars, da_type, da_method)

# filter here to NB, LRclusters, NB-GLMM
dat0 = filter(dat, da_method %in% c('LR', 'negbinom', 'nebula')) %>% 
    # re-color
    mutate(color = ifelse(da_type == 'NBGMM', 'random effect (GLMM)',
                          ifelse(latent_vars, 'fixed effect (latent.vars)',
                                 'default')))
# remove fixed columns
keep = map_lgl(dat0, ~ n_distinct(.x) > 1)
dat0 %<>% extract(, keep)

# re-code x-values and colors
dat0 %<>%
    mutate(method = paste0(da_method, '-', da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
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
               as.character())
# print 
dplyr::count(dat0, color, method) %>% arrange(color, method)

# check number of regions is the same 
dat0 %>% 
    group_by(label) %>% 
    summarise(n_distinct_n_regions = n_distinct(number_of_regions)) %>% 
    filter(n_distinct_n_regions > 1)
arrange(dat0, label, method, color) %>% 
    dplyr::select(-number_of_da_regions, -number_of_nominal_da_regions) ## OK

###############################################################################-
## Number of DA regions, boxplot ####
###############################################################################-

labs = dat0 %>% 
    group_by(method, color) %>%
    summarise(mean = mean(number_of_da_regions), 
              n = n(), 
              median = median(number_of_da_regions), 
              number_of_da_regions = median) %>%
    ungroup() %>%
    mutate(label = round(median))
pal = c('fixed effect (latent.vars)' = colours.cafe425[1],
        'random effect (GLMM)' = 'grey80',
        da_analysis_colors
        )
p5 = dat0 %>%
    mutate(color = ifelse(is.na(latent_vars), method, color)) %>%
    ggplot(aes(x = reorder(method, number_of_da_regions, stats::median), 
               color = color, fill = color)) +
    geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, 
                 alpha = 0.4, width = 0.6, size = 0.35,
                 position = position_dodge(width = 0.75)) + 
    # geom_jitter(aes(y = aucc), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -12.5e3, group = color, label = label), size = 1.75, hjust = 0,
              color = 'black', position = position_dodge(width = 0.75)) +
    scale_color_manual('', values = pal, breaks = c('fixed effect (latent.vars)', 'random effect (GLMM)'), labels = Hmisc::capitalize) +
    scale_fill_manual('', values = pal, breaks = c('fixed effect (latent.vars)', 'random effect (GLMM)'), labels = Hmisc::capitalize) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous(expression('# of DA peaks'~(10^3)), 
                       labels = function(x) x / 1e3,
                       breaks = seq(0, 80, 10) * 1e3) +
    coord_flip(ylim = c(-12e3, 55e3)) +
    guides(color = guide_legend(ncol = 1)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 0.6,
          legend.position = 'top',
          legend.key.width = unit(0.35, 'lines'),
          legend.key.height = unit(0.5, 'lines')) +
    guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))
p5
ggsave("fig/Supp_Fig4/false-discoveries-Luecken-peaks-latent.vars.pdf", p5,
       width = 4.5, height = 6.5, units = "cm", useDingbats = FALSE)

top_row = plot_grid(p1, p2, nrow=2)
ggsave('fig/Supp_Fig4/top_row.pdf', top_row, width = 18, height = 9.5, units = 'cm')

