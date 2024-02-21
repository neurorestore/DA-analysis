# Plot a summary figure with performance across all experiments.
setwd("C:/Users/teo/Documents/EPFL/projects/DA-analysis/")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(ggsignif)
source("R/theme.R")

cell_count_df = readRDS('data/summaries/misc/comparisons.rds') %>%
    mutate(
        cell_type = paste0(compar1, '_vs_', compar2),
        cell_type = ifelse(paper == 'GonzalesBlas_2018', paste0('Melanoma cell line-', cell_type), cell_type),
        cell_type = ifelse(paper == 'Pliner_2018', paste0('HSMM-', cell_type), cell_type)
    )

dat0 = readRDS('data/summaries/misc/mem_summary.rds') %>%
    mutate(method = paste0(da_method, '-', da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = da_family,
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', color) %>% 
               fct_recode('single-cell' = 'singlecell',
                          'other' = 'non_libra') %>% 
               fct_relevel('single-cell', 'pseudobulk', 'other'),
           method = fct_recode(method,
                               "'SnapATAC::findDAR'" = 'SnapATAC::findDAR',
                               'Permutation~test' = 'permutation',
                               't~test' = 't-t',
                               'LR[peaks]' = 'LRpeaks',
                               'LR[clusters]' = 'LR-LR',
                               'Binomial' = 'binomial',
                               'Wilcoxon~rank-sum~test' = 'wilcox-wilcox',
                               'Fisher~exact~test' = 'FET',
                               'Negative~binomial' = 'negbinom-negbinom',
                               "'ArchR::t test'" = 'ArchR_t-ArchR_t',
                               "'ArchR::Wilcoxon rank-sum test'" = 'ArchR_wilcoxon-ArchR_wilcoxon',
                               "'glmGamPoi::Negative binomial'"='glmGamPoi-glmGamPoi'
                               ) %>%
               as.character()
           ) %>% 
    # drop limma
    filter(!grepl('limma', method), exp == 'scATAC') %>%
    mutate(
        time = Elapsed_Time_sec/60,
        mem = Peak_RAM_Used_MiB/1000
    )
pal = da_analysis_colors

dat0 %>% dplyr::select(method) %>% table()
dat0 %<>%
    left_join(
        cell_count_df %>%
            dplyr::select(paper, cell_type, cell_count),
        by = c('paper', 'cell_type')
    )

# Experiment 1 first
dat1 = dat0 %>% 
    filter(standard == 'golden', exp == 'scATAC') %>%
    filter(!grepl('ArchR', method), !grepl('glmGamPoi', method))
dat1 %>% dplyr::select(method) %>% table()

# check time

p1 = dat1 %>%
    ggplot(aes(x = log10(cell_count), y = log10(time), group=method, color=method, fill=method, label = method)) +
    # geom_line() +
    geom_point(size=0.01) +
    # geom_smooth(se=F, method='lm', size=0.5) +
    # geom_text(data=labs, aes(label = method, x = cell_count, y = log10(time)), hjust=0) +
    boxed_theme() +
    scale_fill_manual(values=pal, labels = parse_format()) +
    scale_color_manual(values=pal, labels = parse_format()) +
    theme(
        aspect.ratio = 1,
        legend.title = element_blank(),
        legend.position = 'none'
    ) +
    ylab(expression(Log[10](runtime~'in'~mins))) +
    xlab(expression(Log[10]('#'~of~cells)))
# p1

p2 = dat1 %>%
    mutate(method = factor(method, levels = names(da_analysis_colors))) %>%
    ggplot(aes(x = log10(cell_count), y = log10(mem), group=method, color=method, fill=method, label = method)) +
    # geom_line() +
    geom_point(size=0.01) +
    # geom_smooth(se=F, method='lm', size=0.5) +
    # geom_text(data=labs, aes(label = method, x = cell_count, y = log10(time)), hjust=0) +
    boxed_theme() +
    scale_fill_manual(values=pal, labels = parse_format()) +
    scale_color_manual(values=pal, labels = parse_format()) +
    theme(
        aspect.ratio = 1,
        legend.title = element_blank(),
        legend.position = 'bottom',
        # legend.justification = 'bottom',
        legend.spacing.y = unit(0.0, 'cm'),
        legend.spacing.x = unit(0.0, 'cm'),
        legend.text.align = 0,
    ) +
    # ggtitle('Scalability of memory') +
    guides(
        color = guide_legend(nrow = 2, byrow=T),
        fill = guide_legend(nrow = 2, byrow=T),
        plot.title = element_text(size=5)
    ) +
    ylab(expression(Log[10](memory~'in'~GB))) +
    xlab(expression(Log[10]('#'~of~cells)))
p2

dat2 = dat0 %>% filter(standard == 'silver', exp == 'scATAC') %>%
    filter(!grepl('ArchR', method), !grepl('glmGamPoi', method))
dat2 %>% dplyr::select(method) %>% table()

# check time
p3 = dat2 %>%
    ggplot(aes(x = log10(cell_count), y = log10(time), group=method, color=method, fill=method, label = method)) +
    # geom_line() +
    geom_point(size=0.01) +
    # geom_smooth(se=F, method='lm', size=0.5) +
    # geom_text(data=labs, aes(label = method, x = cell_count, y = log10(time)), hjust=0) +
    boxed_theme() +
    scale_fill_manual(values=pal, labels = parse_format()) +
    scale_color_manual(values=pal, labels = parse_format()) +
    theme(
        aspect.ratio = 1,
        legend.title = element_blank(),
        legend.position = 'none'
    ) +
    # ggtitle('Scalability of time') +
    guides(
        color = guide_legend(ncol = 1, byrow=T),
        fill = guide_legend(ncol = 1, byrow=T),
        plot.title = element_text(size=5)
    ) +
    ylab(expression(Log[10](runtime~'in'~mins))) +
    xlab(expression(Log[10]('#'~of~cells)))
p3

p4 = dat2 %>%
    ggplot(aes(x = log10(cell_count), y = log10(mem), group=method, color=method, fill=method, label = method)) +
    # geom_line() +
    geom_point(size=0.01) +
    # geom_smooth(se=F, method='lm', size=0.5) +
    # geom_text(data=labs, aes(label = method, x = cell_count, y = log10(time)), hjust=0) +
    boxed_theme() +
    scale_fill_manual(values=pal, labels = parse_format()) +
    scale_color_manual(values=pal, labels = parse_format()) +
    theme(
        aspect.ratio = 1,
        legend.title = element_blank(),
        legend.position = 'none'
    ) +
    ylab(expression(Log[10](memory~'in'~GB))) +
    xlab(expression(Log[10]('#'~of~cells)))
p4

p5 = p1 | p2 | p3 | p4
ggsave("fig/Supp_Fig28/full.pdf", p5, width = 18, height = 8, units = "cm", useDingbats = FALSE)
 