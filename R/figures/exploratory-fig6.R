setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

bin_bulk_aucc_rank = readRDS('data/summaries/meta_summaries/bulk_aucc_binarized_summary.rds')
bin_bulk_aucc_rank %<>%
    mutate(bin = ifelse(color == 'binarized', 'binarized', 'non-binarized')) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, bin) %>%
    mutate(
        norm = '',
        method_full = paste0(method, '_', bin, '_', norm)
    )

bulk_aucc_rank = readRDS('data/summaries/meta_summaries/bulk_aucc_summary.rds')
bulk_aucc_rank %<>%
    filter(!method %in% bin_bulk_aucc_rank$method) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val) %>%
    mutate(
        bin = 'n/a',
        norm = '',
        method_full = paste0(method, '_', bin, '_', norm)
    )

norm_bulk_aucc_rank = readRDS('data/summaries/meta_summaries/bulk_aucc_normalized_summary.rds')
norm_bulk_aucc_rank %<>%
    mutate(norm = as.character(normalization)) %>% 
    filter(norm != 'logTP(10K)') %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, norm) %>%
    mutate(
        bin = ifelse(method %in% bin_bulk_aucc_rank$method, 'non-binarized', 'n/a'),
        method_full = paste0(method, '_', bin, '_', norm)
    ) %>%
    dplyr::relocate(bin, .before=norm)

dat1 = rbind(bulk_aucc_rank, bin_bulk_aucc_rank, norm_bulk_aucc_rank) %>%
    distinct() %>%
    mutate(
        method = ifelse(method == 'LR_peaks', 'LR[peaks]', method),
        norm = ifelse(norm == '', 'n/a', norm),
        norm = ifelse(method %in% norm_bulk_aucc_rank$method & norm == 'n/a', 'logTP(10K)', norm)) %>%
    distinct() %>%
    mutate(
        rank = rank(-val, ties.method='first'),
        tercile = ntile(-rank, 3),
        type = 'Matched bulk', 
        experiment = 'Experiment 1: Bulk concordance'
    ) %>%
    arrange(rank)

luecken_false_discoveries_rank = readRDS('data/summaries/meta_summaries/luecken_false_discoveries_summary.rds')
luecken_false_discoveries_rank %<>%
    filter(!method %in% bin_bulk_aucc_rank$method) %>% 
    dplyr::rename(val = median) %>%
    dplyr::select(method, val) %>%
    mutate(
        bin = 'n/a',
        norm = '',
        method_full = paste0(method, '_', bin, '_', norm)
    )

bin_luecken_false_discoveries_rank = readRDS('data/summaries/meta_summaries/luecken_false_discoveries_binarized_summary.rds')
bin_luecken_false_discoveries_rank %<>%
    mutate(bin = ifelse(color == 'binarized', 'binarized', 'non-binarized')) %>%
    dplyr::rename(val = median) %>%
    dplyr::select(method, val, bin) %>%
    mutate(
        norm = '',
        method_full = paste0(method, '_', bin, '_', norm)
    )

norm_luecken_false_discoveries_rank = readRDS('data/summaries/meta_summaries/luecken_false_discoveries_normalized_summary.rds')
norm_luecken_false_discoveries_rank %<>%
    mutate(norm = as.character(normalization)) %>%
    filter(norm != 'logTP(10K)') %>%
    dplyr::rename(val = median) %>%
    dplyr::select(method, val, norm) %>%
    mutate(
        bin = ifelse(method %in% bin_bulk_aucc_rank$method, 'non-binarized', 'n/a'),
        method_full = paste0(method, '_', bin, '_', norm)
    ) %>%
    dplyr::relocate(bin, .before=norm)

dat2 = rbind(luecken_false_discoveries_rank, bin_luecken_false_discoveries_rank, norm_luecken_false_discoveries_rank) %>%
    distinct() %>%
    mutate(
        method = ifelse(method == 'LR_peaks', 'LR[peaks]', method),
        norm = ifelse(norm == '', 'n/a', norm),
        norm = ifelse(method %in% norm_bulk_aucc_rank$method & norm == 'n/a', 'logTP(10K)', norm)) %>%
    distinct() %>%
    mutate(
        rank = rank(val, ties.method='first'),
        tercile = ntile(-rank, 3),
        type = 'Published data',
        experiment = 'Experiment 3: False discoveries'
    ) %>%
    arrange(rank)

snapatac_false_discoveries_rank = readRDS('data/summaries/meta_summaries/snapatac_false_discoveries_summary.rds')
snapatac_false_discoveries_rank %<>%
    filter(!method %in% bin_bulk_aucc_rank$method) %>% 
    dplyr::rename(val = median) %>%
    dplyr::select(method, val) %>%
    mutate(
        bin = 'n/a',
        norm = '',
        method_full = paste0(method, '_', bin, '_', norm)
    )

bin_snapatac_false_discoveries_rank = readRDS('data/summaries/meta_summaries/snapatac_false_discoveries_binarized_summary.rds')
bin_snapatac_false_discoveries_rank %<>%
    mutate(bin = ifelse(color == 'binarized', 'binarized', 'non-binarized')) %>%
    dplyr::rename(val = median) %>%
    dplyr::select(method, val, bin) %>%
    mutate(
        norm = '',
        method_full = paste0(method, '_', bin, '_', norm)
    )

norm_snapatac_false_discoveries_rank = readRDS('data/summaries/meta_summaries/snapatac_false_discoveries_normalized_summary.rds')
norm_snapatac_false_discoveries_rank %<>%
    mutate(norm = as.character(normalization)) %>%
    filter(norm != 'logTP(10K)') %>%
    dplyr::rename(val = median) %>%
    dplyr::select(method, val, norm) %>%
    mutate(
        bin = ifelse(method %in% bin_bulk_aucc_rank$method, 'non-binarized', 'n/a'),
        method_full = paste0(method, '_', bin, '_', norm)
    ) %>%
    dplyr::relocate(bin, .before=norm)

dat3 = rbind(snapatac_false_discoveries_rank, bin_snapatac_false_discoveries_rank, norm_snapatac_false_discoveries_rank) %>%
    distinct() %>%
    mutate(
        method = ifelse(method == 'LR_peaks', 'LR[peaks]', method),
        norm = ifelse(norm == '', 'n/a', norm),
        norm = ifelse(method %in% norm_bulk_aucc_rank$method & norm == 'n/a', 'logTP(10K)', norm)) %>%
    distinct() %>%
    mutate(
        rank = rank(val, ties.method='first'),
        tercile = ntile(-rank, 3),
        type = 'Downsampled data',
        experiment = 'Experiment 3: False discoveries'
    ) %>%
    arrange(rank)

splatter_false_discoveries_rank = readRDS('data/summaries/meta_summaries/splatter_false_discoveries_summary.rds')
splatter_false_discoveries_rank %<>%
    filter(!method %in% bin_bulk_aucc_rank$method) %>% 
    dplyr::rename(val = median) %>%
    dplyr::select(method, val) %>%
    mutate(
        bin = 'n/a',
        norm = '',
        method_full = paste0(method, '_', bin, '_', norm)
    )

bin_splatter_false_discoveries_rank = readRDS('data/summaries/meta_summaries/splatter_false_discoveries_binarized_summary.rds')
bin_splatter_false_discoveries_rank %<>%
    mutate(bin = ifelse(color == 'binarized', 'binarized', 'non-binarized')) %>%
    dplyr::rename(val = median) %>%
    dplyr::select(method, val, bin) %>%
    mutate(
        norm = '',
        method_full = paste0(method, '_', bin, '_', norm)
    )

norm_splatter_false_discoveries_rank = readRDS('data/summaries/meta_summaries/splatter_false_discoveries_normalized_summary.rds')
norm_splatter_false_discoveries_rank %<>%
    mutate(norm = as.character(normalization)) %>%
    filter(norm != 'logTP(10K)') %>%
    dplyr::rename(val = median) %>%
    dplyr::select(method, val, norm) %>%
    mutate(
        bin = ifelse(method %in% bin_bulk_aucc_rank$method, 'non-binarized', 'n/a'),
        method_full = paste0(method, '_', bin, '_', norm)
    ) %>%
    dplyr::relocate(bin, .before=norm)

dat4 = rbind(splatter_false_discoveries_rank, bin_splatter_false_discoveries_rank, norm_splatter_false_discoveries_rank) %>%
    distinct() %>%
    mutate(
        method = ifelse(method == 'LR_peaks', 'LR[peaks]', method),
        norm = ifelse(norm == '', 'n/a', norm),
        norm = ifelse(method %in% norm_bulk_aucc_rank$method & norm == 'n/a', 'logTP(10K)', norm)) %>%
    distinct() %>%
    mutate(
        rank = rank(val, ties.method='first'),
        tercile = ntile(-rank, 3),
        type = 'Simulations',
        experiment = 'Experiment 3: False discoveries'
    ) %>%
    arrange(rank)

dat = rbind(
    dat1, 
    dat2,
    dat3,
    dat4
)

method_order = dat1 %>% arrange(rank) %>% pull(method_full) %>% rev()
dat$method_full = factor(dat$method_full, levels=method_order)
dat1$method_full = factor(dat1$method_full, levels=method_order)
# final_rank$method_bin = factor(final_rank$method_bin, levels=method_order)

factor_levels = rev(unique(dat$type))
dat$type = factor(dat$type, levels=factor_levels)

dat %<>%
    mutate(
        tercile = as.character(tercile),
        tercile = fct_recode(tercile,
            'Top' = '3',
            'Middle' = '2',
            'Bottom' = '1'
        )
    )
dat$tercile = factor(dat$tercile, levels = c('Bottom', 'Middle', 'Top'))
dat$experiment = factor(dat$experiment, levels = c(
    'Experiment 1: Bulk concordance',
    'Experiment 3: False discoveries'
    )
)

pal = colorRampPalette(c(brewer.pal(3, 'Reds')[1], "#B30F17"))(3)

experiments = levels(dat$experiment)
table(dat$tercile)
bin_colors = c('#F38400', '#2B3D26', 'white') %>% setNames(c('binarized', 'non-binarized', 'n/a'))
p0_1 = dat1 %>%
    mutate(bin = factor(bin, levels = names(bin_colors))) %>%
    ggplot(aes(x = method_full, y = 1, fill = bin)) +
    geom_tile(color = 'white') +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual('Binarization', values = bin_colors) +
    # coord_fixed() +
    coord_fixed(ratio=1) +
    boxed_theme() +
    theme(
        plot.margin = unit(c(0.01, 0.01, 0.2, 0.01), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face='bold'),
        legend.position = 'bottom',
        legend.justification = 'right',
        legend.title = element_blank(),
        # legend.justification = 'right',
        legend.key.size = unit(0.35, 'lines'),
        plot.title = element_text(size=7)
    ) +
    guides(fill = guide_legend(nrow = 1))
p0_1
norm_colors = c(
    pals::kelly(14)[8:14] %>% setNames(c('none', 'logTP(10K)', 'TP(median)', 'TP(10K)', 'logTP(median)', 'TF-IDF', 'smooth GC-full quantile')),
    c('white') %>% setNames('n/a')
)
p0_2 = dat1 %>%
    mutate(norm = factor(norm, levels = names(norm_colors))) %>%
    ggplot(aes(x = method_full, y = 1, fill = norm)) +
    geom_tile(color = 'white') +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual('Normalization', values = norm_colors) +
    # coord_fixed() +
    coord_fixed(ratio=1) +
    boxed_theme() +
    theme(
        plot.margin = unit(c(0.01, 0.01, 0.2, 0.01), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face='bold'),
        legend.position = 'bottom',
        legend.justification = 'right',
        legend.title = element_blank(),
        # legend.justification = 'right',
        legend.key.size = unit(0.35, 'lines'),
        plot.title = element_text(size=7)
    ) +
    guides(fill = guide_legend(nrow = 1))
p0_2

for (experiment in experiments) {
    if (experiment == first(experiments)) out_plots = list(p0_1, p0_2)
    dat_tmp = dat %>%
        filter(experiment == !!experiment)
    p1_1 = dat_tmp %>%
        ggplot(aes(x = method_full, y = type, fill=tercile)) +
        # facet_grid(experiment ~ ., scales = 'free', space = 'free') +
        geom_tile(color = 'white') +
        scale_x_discrete(expand = c(0, 0), labels = scales::parse_format()(rev(dat1$method))) +
        # scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_fill_manual('Tercile', values = pal,
                          labels = c('1' = 'Bottom', '2' = 'Middle', '3' = 'Top')) +
        # coord_fixed() +
        ggtitle(experiment) +
        coord_fixed(ratio=1) +
        boxed_theme() +
        theme(
            plot.margin = unit(c(0.01, 0.01, 0.2, 0.01), "cm"),
            axis.text.x = element_text(angle = 45, hjust = 1),  
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.text = element_text(face='bold'),
            legend.position = 'bottom',
            legend.justification = 'right',
            legend.key.size = unit(0.35, 'lines'),
            plot.title = element_text(size=7)
            )
    if (experiment != last(experiments)) {
        p1_1 = p1_1 + theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = 'none'
            )
    }
    out_plots[[length(out_plots)+1]] = p1_1
}

p1 = wrap_plots(out_plots, ncol=1)
p1
ggsave('fig/Exp_Fig6/binarization-normalization-terciles.pdf', p1, width=15, height=15, units='cm')

summary_rank = dat1 %>%
    group_by(method) %>%
    arrange(method, -val) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    mutate(
        rank = rank(-val, ties.method='first'),
        tercile = ntile(-rank, 3),
        type = 'Matched bulk', 
        experiment = 'Experiment 1: Bulk concordance'
    ) %>%
    arrange(rank)

method_order = summary_rank %>% arrange(rank) %>% pull(method_full) %>% rev()
summary_dat = rbind(
    dat1 %>% 
        filter(method_full %in% method_order) %>% 
        mutate(
            rank = rank(-val, ties.method='first'),
            tercile = ntile(-rank, 3)
        ),
    dat2 %>% 
        filter(method_full %in% method_order) %>% 
        mutate(
            rank = rank(val, ties.method='first'),
            tercile = ntile(-rank, 3)
        ),
    dat3 %>% 
        filter(method_full %in% method_order) %>% 
        mutate(
            rank = rank(val, ties.method='first'),
            tercile = ntile(-rank, 3)
        ),
    dat4 %>% 
        filter(method_full %in% method_order) %>% 
        mutate(
            rank = rank(val, ties.method='first'),
            tercile = ntile(-rank, 3)
        )
)
summary_dat$method_full = factor(summary_dat$method_full, levels=method_order)
summary_dat %<>%
    mutate(
        tercile = as.character(tercile),
        tercile = fct_recode(tercile,
                             'Top' = '3',
                             'Middle' = '2',
                             'Bottom' = '1'
        )
    )
summary_dat$tercile = factor(summary_dat$tercile, levels = c('Bottom', 'Middle', 'Top'))
summary_dat$experiment = factor(summary_dat$experiment, levels = c(
    'Experiment 1: Bulk concordance',
    'Experiment 3: False discoveries'
    )
)

pal = colorRampPalette(c(brewer.pal(3, 'Reds')[1], "#B30F17"))(3)

p0_1 = summary_rank %>%
    mutate(bin = factor(bin, levels = names(bin_colors))) %>%
    ggplot(aes(x = method_full, y = 1, fill = bin)) +
    geom_tile(color = 'white') +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual('Binarization', values = bin_colors) +
    # coord_fixed() +
    coord_fixed(ratio=1) +
    boxed_theme() +
    theme(
        plot.margin = unit(c(0.01, 0.01, 0.2, 0.01), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face='bold'),
        legend.position = 'right',
        legend.title = element_blank(),
        # legend.justification = 'right',
        legend.key.size = unit(0.35, 'lines'),
        plot.title = element_text(size=7)
    )
p0_1

p0_2 = summary_rank %>%
    mutate(norm = factor(norm, levels = names(norm_colors))) %>%
    ggplot(aes(x = method_full, y = 1, fill = norm)) +
    geom_tile(color = 'white') +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual('Normalization', values = norm_colors) +
    # coord_fixed() +
    coord_fixed(ratio=1) +
    boxed_theme() +
    theme(
        plot.margin = unit(c(0.01, 0.01, 0.2, 0.01), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face='bold'),
        legend.position = 'right',
        legend.title = element_blank(),
        # legend.justification = 'right',
        legend.key.size = unit(0.35, 'lines'),
        plot.title = element_text(size=7)
    )
p0_2

for (experiment in experiments) {
    if (experiment == first(experiments)) out_plots = list(p0_1, p0_2)
    dat_tmp = summary_dat %>%
        filter(experiment == !!experiment)
    p1_1 = dat_tmp %>%
        ggplot(aes(x = method_full, y = type, fill=tercile)) +
        # facet_grid(experiment ~ ., scales = 'free', space = 'free') +
        geom_tile(color = 'white') +
        scale_x_discrete(expand = c(0, 0), labels = scales::parse_format()(rev(summary_rank$method))) +
        # scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_fill_manual('Tercile', values = pal,
                          labels = c('1' = 'Bottom', '2' = 'Middle', '3' = 'Top')) +
        # coord_fixed() +
        ggtitle(experiment) +
        coord_fixed(ratio=1) +
        boxed_theme() +
        theme(
            plot.margin = unit(c(0.01, 0.01, 0.2, 0.01), "cm"),
            axis.text.x = element_text(angle = 45, hjust = 1),  
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.text = element_text(face='bold'),
            legend.position = 'bottom',
            legend.justification = 'right',
            legend.key.size = unit(0.35, 'lines'),
            plot.title = element_text(size=7)
        )
    if (experiment != last(experiments)) {
        p1_1 = p1_1 + theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = 'none'
        )
    }
    out_plots[[length(out_plots)+1]] = p1_1
}

p2 = wrap_plots(out_plots, ncol=1)
p2
ggsave('fig/Exp_Fig6/best-combination.pdf', p2, width=8, height=7, units='cm')
