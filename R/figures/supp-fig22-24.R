# Plot a summary figure with performance across all experiments.
setwd("C:/Users/teo/Documents/EPFL/projects/DA-analysis/")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(pals)
library(upstartr)
library(effsize)
source('R/theme.R')

# load data
dat = readRDS("data/summaries/false_discoveries/false_discoveries_snapatac.rds") %>% 
    type_convert() %>% 
    # flag study
    dplyr::select(-ori_filename) %>% 
    # filter out methods we are not using
    filter(da_family != 'mixedmodel',
           pseudo_repl == FALSE,
           !is.na(percentile)) %>% 
    # entire range of expression values only - no binarization
    filter(!binarization) %>% 
    # normalization: only LR, t, wilcox
    filter(da_method %in% c('LR', 't', 'wilcox'))

keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-singlecell') %>% 
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
dplyr::count(dat, method) 

# recode normalization
dat %<>% 
    # re-name normalization methods
    mutate(normalization = fct_recode(normalization,
                                      'logTP(10K)' = 'f',
                                      'logTP(median)' = 'log_tp_median',
                                      'smooth GC-full quantile' = 'qsmooth',
                                      'TF-IDF' = 'TFIDF',
                                      'TP(median)' = 'tp_median',
                                      'TP(10K)' = 'tp10k'))
# remove limma
dat %<>% filter(!grepl('limma', da_method))

## datasets
dplyr::count(dat, percentile, percentile_type)

percentile_types = c('mean', 'percentage_open', 'peak_size')
grey_pal = colorRampPalette(c('grey10', 'grey85'))(6) %>% 
    setNames(c('none', 'TP(median)', 'TP(10K)', 'logTP(median)', 'smooth GC-full quantile', 'TF-IDF'))
pal = c(da_analysis_colors[c('LR[clusters]', 't~test', 'Wilcoxon~rank-sum~test')], grey_pal)
# iterate through peak types
for (percentile_type in percentile_types) {
    dat0 = filter(dat, percentile_type == !!percentile_type)
    dat0$method = factor(dat0$method, levels = c('LR[clusters]', 't~test', 'Wilcoxon~rank-sum~test'))
    
    ###############################################################################-
    ## Boxplots, un-normalized number of DA regions ####
    ###############################################################################-
    
    percentile_type_name = case_when(
        percentile_type == 'mean' ~ 'Read depth',
        percentile_type == 'percentage_open' ~ '% open',
        percentile_type == 'peak_size' ~ 'Peak width (bp)'
    )
    
    normalization_methods = c('logTP(10K)', names(grey_pal))
    
    curr_out_plots = list()
    for (normalization in normalization_methods) {
        if (normalization == 'logTP(10K)') {
            p1_1 = dat0 %>% 
                filter(normalization == !!normalization) %>%
                mutate(decile = gsub("^.*_", "", percentile) %>% as.integer() %>% factor()) %>% 
                mutate(color = ifelse(normalization == 'logTP(10K)', as.character(method), as.character(normalization))) %>%
                ggplot(aes(x = decile, y = number_of_da_regions, fill = color, color = color)) +
                # facet_grid(normalization ~ method, labeller = label_parsed) +
                facet_grid(~ method, labeller = label_parsed) +
                # ggtitle(normalization) +
                geom_boxplot(size = 0.35, width = 0.6, alpha = 0.4, outlier.shape = NA) +
                scale_x_discrete('Decile', breaks=c(1,10)) +
                scale_y_continuous('# of DA peaks') +
                scale_color_manual('', values = pal) +
                scale_fill_manual('', values = pal) +
                boxed_theme() +
                theme(aspect.ratio = 1,
                      legend.position = 'none',
                      plot.title = element_text(size=5)
                )
        } else {
            p1_1 = dat0 %>% 
                filter(normalization == !!normalization) %>%
                mutate(decile = gsub("^.*_", "", percentile) %>% as.integer() %>% factor()) %>% 
                mutate(color = ifelse(normalization == 'logTP(10K)', as.character(method), as.character(normalization))) %>%
                ggplot(aes(x = decile, y = number_of_da_regions, fill = color, color = color)) +
                facet_grid(normalization ~ method) +
                # facet_grid(~ method, labeller = label_parsed) +
                # ggtitle(normalization) +
                geom_boxplot(size = 0.35, width = 0.6, alpha = 0.4, outlier.shape = NA) +
                scale_x_discrete('Decile', breaks=c(1,10)) +
                scale_y_continuous('# of DA peaks') +
                scale_color_manual('', values = pal) +
                scale_fill_manual('', values = pal) +
                boxed_theme() +
                theme(aspect.ratio = 1,
                      legend.position = 'none',
                      plot.title = element_text(size=5)
                )   
        }
        p1_1
        if (normalization != first(normalization_methods)) {
            p1_1 = p1_1 + theme(strip.text.x = element_blank())
        }
        if (normalization != last(normalization_methods)) {
            p1_1 = p1_1 + theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank()
                )
        }
        if (normalization != 'TP(10K)') {
            p1_1 = p1_1 + theme(
                axis.title.y = element_blank()
            )
        }
        # print(p1_1)
        curr_out_plots[[length(curr_out_plots)+1]] = p1_1
    }
    p1 = wrap_plots(curr_out_plots, ncol=1)
    ggsave(paste0("fig/Supp_Fig22-24/snapatac-normalization-boxplot-number-", percentile_type, ".pdf"), p1, width = 9, height = 15, units = "cm", useDingbats = FALSE)
    
    delta = dat0 %>% 
        group_by(da_method, seed, percentile) %>% 
        mutate(delta = number_of_da_regions - 
                   number_of_da_regions[normalization == 'logTP(10K)'],
               n = n()) %>% 
        ungroup()
    
    grid = distinct(delta, method, normalization, percentile) %>%
        filter(normalization != 'logTP(10K)')
    eff = pmap_dfr(grid, function(...) {
        current = tibble(...)
        norm = as.character(current$normalization)
        df = filter(dat0, method == current$method, percentile == current$percentile,
                    normalization %in% c(norm, 'logTP(10K)')) %>% 
            droplevels() %>% 
            mutate(normalization = factor(normalization, levels = c(norm, 'logTP(10K)')))
        eff = cohen.d(formula = number_of_da_regions ~ normalization | Subject(seed),
                      data = df,
                      paired = TRUE)$estimate
        mutate(current, d = eff)
    })
    
    eff %<>% mutate(d = ifelse(is.na(d), 0, d))
    
    curr_out_plots = list()
    normalization_methods = c(names(grey_pal))
    for (normalization in normalization_methods) {
        if (normalization == first(normalization_methods)) {
            p2_1 = eff %>% 
                filter(normalization == !!normalization) %>%
                mutate(decile = gsub("^.*_", "", percentile) %>% as.integer() %>% factor()) %>% 
                mutate(color = ifelse(normalization == 'logTP(10K)', as.character(method), as.character(normalization))) %>%
                ggplot(aes(x = decile, y = d, fill = color, color = color)) +
                facet_grid(normalization ~ method, labeller = label_parsed) +
                # facet_grid(~ method, labeller = label_parsed) +
                # ggtitle(normalization) +
                geom_point(position = position_dodge(width = 0.6), shape = 21, size = 0.8) + 
                geom_errorbar(aes(ymin = 0, ymax = d, group = method), width = 0, size = 0.3,
                              position = position_dodge(width = 0.6)) +
                scale_y_continuous('Cohen\'s d, relative to log(TP10K)') +
                scale_x_discrete('Decile', breaks=c(1,10)) +
                geom_hline(yintercept=0, linetype='dotted', col = 'black') +
                scale_color_manual('', values = pal) +
                scale_fill_manual('', values = pal) +
                boxed_theme() +
                theme(aspect.ratio = 1,
                      legend.position = 'none',
                      plot.title = element_text(size=5)
                )
        } else {
            p2_1 = eff %>% 
                filter(normalization == !!normalization) %>%
                mutate(decile = gsub("^.*_", "", percentile) %>% as.integer() %>% factor()) %>% 
                mutate(color = ifelse(normalization == 'logTP(10K)', as.character(method), as.character(normalization))) %>%
                ggplot(aes(x = decile, y = d, fill = color, color = color)) +
                facet_grid(normalization ~ method) +
                geom_point(position = position_dodge(width = 0.6), shape = 21, size = 0.8) + 
                geom_errorbar(aes(ymin = 0, ymax = d, group = method), width = 0, size = 0.3,
                              position = position_dodge(width = 0.6)) +
                scale_y_continuous('Cohen\'s d, relative to log(TP10K)') +
                scale_x_discrete('Decile', breaks=c(1,10)) +
                geom_hline(yintercept=0, linetype='dotted', col = 'black') +
                scale_color_manual('', values = pal) +
                scale_fill_manual('', values = pal) +
                boxed_theme() +
                theme(aspect.ratio = 1,
                      legend.position = 'none',
                      plot.title = element_text(size=5)
                )   
        }
        p2_1
        if (normalization != first(normalization_methods)) {
            p2_1 = p2_1 + theme(strip.text.x = element_blank())
        }
        if (normalization != last(normalization_methods)) {
            p2_1 = p2_1 + theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank()
            )
        }
        
        if (normalization != 'TP(10K)') {
            p2_1 = p2_1 + theme(
                axis.title.y = element_blank()
            )
        }
        # print(p2_1)
        curr_out_plots[[length(curr_out_plots)+1]] = p2_1
    }
    p2 = wrap_plots(curr_out_plots, ncol=1)
    p2
    ggsave(paste0("fig/Supp_Fig22-24/snapatac-normalization-boxplot-number-", percentile_type, "-d.pdf"), p2, width = 9, height = 13, units = "cm", useDingbats = FALSE)
}
