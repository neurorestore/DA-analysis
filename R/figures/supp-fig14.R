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
dat = readRDS("data/summaries/false_discoveries/false_discoveries_luecken.rds") %>% 
    type_convert() %>% 
    # flag study
    mutate(study = gsub("-.*$", "", filename)) %>% 
    dplyr::select(-filename) %>% 
    # no pseudoreplicates
    filter(!pseudo_repl) %>% 
    # only percentiles
    filter(!is.na(percentile), normalization == 'f', is.na(latent_vars))
# remove non-unique columns
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

# remove limma
dat %<>% filter(!grepl('limma', da_method))

# print some summaries
## DA methods
dplyr::count(dat, da_family, da_method, da_type, pseudobulk, binarization)
## datasets
dplyr::count(dat, percentile, percentile_type)

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

dat$percentile_type = factor(dat$percentile_type, levels = c('mean', 'percentage_open', 'peak_size'))
out_plots = list()
# iterate through peak types
for (percentile_type in levels(dat$percentile_type)) {
    dat0 = dat %>% 
        filter(percentile_type == !!percentile_type) %>%
        mutate(method_color = ifelse(color == 'binarized', 'binarized', method))
    
    ###############################################################################-
    ## Boxplots, un-normalized number of DA regions ####
    ###############################################################################-
    pal = c(da_analysis_colors, 'binarized' = "#6E645F")
    p1 = dat0 %>% 
        mutate(decile = gsub("^.*_", "", percentile) %>% as.integer() %>% factor()) %>% 
        ggplot(aes(x = decile, y = number_of_da_regions, fill = method_color, color = method_color)) +
        facet_wrap(~ method, nrow = 2, labeller = label_parsed) +
        geom_boxplot(size = 0.35, width = 0.6, alpha = 0.4, outlier.shape = NA) +
        scale_x_discrete('Decile', breaks = c(1,10)) +
        scale_y_continuous('# of DA peaks') +
        scale_color_manual('', values = pal) +
        scale_fill_manual('', values = pal) +
        boxed_theme() +
        theme(aspect.ratio = 1,
              legend.position = 'none')
    p1
    ggsave(paste0("fig/Supp_Fig14/luecken-boxplot-number-", percentile_type, ".pdf"), p1,
           width = 14, height = 7, units = "cm", useDingbats = FALSE)
    out_plots[[length(out_plots)+1]] = p1
    
    ###############################################################################-
    ## Boxplots, proportion number of DA regions in each bin ####
    ###############################################################################-
    
    prop = dat0 %>% 
        mutate(decile = gsub("^.*_", "", percentile) %>% as.integer() %>% factor()) %>% 
        group_by(method, label, binarization) %>% 
        mutate(prop = number_of_da_regions / sum(number_of_da_regions)) %>% 
        ungroup() %>% 
        replace_na(list(prop = 0))
    p2 = prop %>% 
        mutate(decile = gsub("^.*_", "", percentile) %>% as.integer() %>% factor()) %>% 
        ggplot(aes(x = decile, y = prop, fill = method_color, color = method_color)) +
        facet_wrap(~ method, nrow = 2, labeller = label_parsed) +
        geom_boxplot(size = 0.35, width = 0.6, alpha = 0.4, outlier.shape = NA) +
        scale_x_discrete('Decile', breaks = c(1,10)) +
        scale_y_continuous('% of DA peaks', labels = x100) +
        scale_color_manual('', values = pal) +
        scale_fill_manual('', values = pal) +
        boxed_theme() +
        theme(aspect.ratio = 1,
              legend.position = 'none')
    # p2
    ggsave(paste0("fig/Supp_Fig14/luecken-boxplot-proportion-", percentile_type, ".pdf"), p2,
           width = 14, height = 7, units = "cm", useDingbats = FALSE)
    out_plots[[length(out_plots)+1]] = p2
}

full = wrap_plots(out_plots, ncol=2)
ggsave(paste0("fig/Supp_Fig14/full.pdf"), full,
       width = 19, height = 15, units = "cm", useDingbats = FALSE)
