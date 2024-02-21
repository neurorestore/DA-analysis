# Plot a summary figure with performance across all experiments.
setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(pals)
library(upstartr)
source('R/theme.R')

# load data
dat = readRDS("data/summaries/false_discoveries/false_discoveries_snapatac.rds") %>% 
    type_convert() %>% 
    # dplyr::select(-filename) %>% 
    # filter out methods we are not using
    filter(da_family != 'mixedmodel',
           pseudo_repl == FALSE) %>%
    # ignore results on full dataset
    drop_na(percentile, percentile_type) %>% 
    # ignore normalization, binarization
    filter(normalization == 'f', binarization == FALSE)

# simulation parameters
table(dat$ncells)
table(dat$fragments_prop)
table(dat$seed)
## DA methods
dplyr::count(dat, da_family, da_method, de_type, pseudobulk, pseudo_repl)

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-', de_type) %>% 
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
    filter(!grepl('limma', de_method))
# print 
dplyr::count(dat, color, method) %>% arrange(color, method)
dat$percentile_type = factor(dat$percentile_type, levels = c('mean', 'percentage_open', 'peak_size'))

pal = da_analysis_colors
# iterate through peak types
for (percentile_type in levels(dat$percentile_type)) {
    if (percentile_type == 'mean') out_plots = list()
    dat0 = dat %>%
        filter(
            percentile_type == !!percentile_type
        ) %>%
        mutate(
            decile = gsub("^.*_", "", percentile) %>% as.integer() %>% factor(),
            method = factor(method, levels = names(da_analysis_colors))
        ) 
    
    percentile_type_name = case_when(
        percentile_type == 'mean' ~ 'Read depth',
        percentile_type == 'percentage_open' ~ '% open',
        percentile_type == 'peak_size' ~ 'Peak width (bp)'
    )
    
    #############################################################################-
    ## Heatmap, un-normalized number of DA peaks ####
    #############################################################################-
    
    avg = dat0 %>% 
        group_by_at(vars(-seed, -number_of_regions, -number_of_da_regions,
                         -number_of_nominal_da_regions)) %>% 
        summarise(
            mean_number_of_da_regions = mean(number_of_da_regions),
            median_number_of_da_regions = median(number_of_da_regions),
            number_of_da_regions = median_number_of_da_regions,
            n = n()) %>% 
        ungroup()
    
    summary = avg %>% dplyr::select("da_family", "method", "decile", "median_number_of_da_regions", "mean_number_of_da_regions")
    saveRDS(summary, paste0('data/summaries/meta_summaries/snapatac_false_discoveries_', percentile_type, '_bias_summary.rds'))
    
    dat0 %<>%
        mutate(
            method = fct_recode(
                method, 
                "'Fisher\nexact test'" = 'Fisher~exact~test',
                "'Permutation\ntest'" = 'Permutation~test',
                "'Negative\nbinomial'" = 'Negative~binomial',
                "'Wilcoxon\nrank-sum\ntest'" = 'Wilcoxon~rank-sum~test',
                "'DESeq2\nLRT'" = 'DESeq2-LRT',
                "'DESeq2\nWald'" = 'DESeq2-Wald',
                "'edgeR\nLRT'" = 'edgeR-LRT',
                "'edgeR\nQLF'" = 'edgeR-QLF',
                "'SnapATAC::\nfindDAR'" = "'SnapATAC::findDAR'"
            )
        )
    
    new_color_names = vapply(names(da_analysis_colors), function(x){
        case_when(
            x == 'Fisher~exact~test' ~ "'Fisher\nexact test'",
            x == 'Permutation~test' ~ "'Permutation\ntest'",
            x == 'Negative~binomial' ~ "'Negative\nbinomial'",
            x == 'Wilcoxon~rank-sum~test' ~ "'Wilcoxon\nrank-sum\ntest'",
            x == 'DESeq2-LRT' ~ "'DESeq2\nLRT'",
            x == 'DESeq2-Wald' ~ "'DESeq2\nWald'",
            x == 'edgeR-LRT' ~ "'edgeR\nLRT'",
            x == 'edgeR-QLF' ~ "'edgeR\nQLF'",
            x == "'SnapATAC::findDAR'" ~"'SnapATAC::\nfindDAR'",
            T ~ as.character(x)
        )
    }, as.character(1)) %>% unname()
    
    names(pal) = new_color_names
    
    range = range(avg$number_of_da_regions)
    brks = range
    p1 = dat0 %>% 
        ggplot(aes(x = decile, y = number_of_da_regions, fill = method, color = method)) +
        facet_wrap(~ method, nrow = 1, labeller = label_parsed) +
        geom_boxplot(size = 0.35, width = 0.6, alpha = 0.4, outlier.shape = NA) +
        scale_x_discrete('Decile', expand = c(0,0)) +
        scale_y_continuous('# of DA regions') +
        scale_color_manual('', values = pal) +
        scale_fill_manual('', values = pal) +
        boxed_theme() +
        theme(
            aspect.ratio = 1,
            legend.position = 'none'
        )
    # p1
    ggsave(paste0("fig/Supp_Fig7/snapatac-boxplot-number-", percentile_type, ".pdf"), p1,
           width = 17, height = 3, units = "cm", useDingbats = FALSE)
    if (percentile_type != first(levels(dat$percentile_type))) {
        p1 = p1 + theme(
            strip.text.x = element_blank()
        )
    }
    p1 = p1 + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
    )
    
    out_plots[[length(out_plots)+1]] = p1
    
    ###############################################################################-
    ## Boxplots, proportion number of DA regions in each bin ####
    ###############################################################################-
    
    prop = dat0 %>% 
        group_by(method, seed) %>% 
        mutate(prop = number_of_da_regions / sum(number_of_da_regions)) %>% 
        ungroup() %>% 
        replace_na(list(prop = 0))
    p2 = prop %>% 
        ggplot(aes(x = decile, y = prop, fill = method, color = method)) +
        facet_wrap(~ method, nrow = 1, labeller = label_parsed) +
        geom_boxplot(size = 0.35, width = 0.6, alpha = 0.4, outlier.shape = NA) +
        scale_x_discrete('Decile', expand = c(0,0)) +
        scale_y_continuous('% of DA regions', labels = x100) +
        scale_color_manual('', values = pal) +
        scale_fill_manual('', values = pal) +
        boxed_theme() +
        theme(
            aspect.ratio = 1,
            legend.position = 'none'
        )
    # p2
    ggsave(paste0("fig/Supp_Fig7/snapatac-boxplot-proportion-", percentile_type, ".pdf"), p2,
           width = 17, height = 3, units = "cm", useDingbats = FALSE)
    if (percentile_type != last(levels(dat$percentile_type))) {
        p2 = p2 + theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank()
        )
    } else {
        p2 = p2 + scale_x_discrete('Decile', expand = c(0,0), breaks = c(1,10))
    }
    p2 = p2 + theme(
        strip.text.x = element_blank()
    )
    out_plots[[length(out_plots)+1]] = p2
}

full = wrap_plots(out_plots, ncol=1)
ggsave(paste0("fig/Supp_Fig7/full.pdf"), full,
       width = 18.5, height = 13, units = "cm", useDingbats = FALSE)
