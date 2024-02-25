setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

### PANEL A ###
bulk_aucc_rank = readRDS('data/summaries/meta_summaries/bulk_aucc_summary.rds')
bulk_aucc_rank %<>%
    mutate(
        rank = rank(-aucc, ties.method='first'),
        method = ifelse(method == 'LR_peaks', 'LR[peaks]', method)
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Matched bulk', 
        # experiment = 'Matched bulk concordance',
        experiment = 'Experiment 1: Bulk concordance',
        tercile = ntile(-rank, 3)
    )

multiome_aucc_rank = readRDS('data/summaries/meta_summaries/multiome_aucc_summary.rds')
multiome_aucc_rank %<>%
    mutate(
        rank = rank(-aucc, ties.method='first')
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Multiome', 
        # experiment = 'Multiome concordance',
        experiment = 'Experiment 2: Multiome concordance',
        tercile = ntile(-rank, 3)
    )

multiome_gsea_aucc_rank = readRDS('data/summaries/meta_summaries/multiome_gsea_aucc_summary.rds')
multiome_gsea_aucc_rank %<>%
    mutate(
        rank = rank(-aucc, , ties.method='first')
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Multiome GO', 
        # experiment = 'Multiome concordance',
        experiment = 'Experiment 2: Multiome concordance',
        tercile = ntile(-rank, 3)
    ) 

luecken_false_discoveries_rank = readRDS('data/summaries/meta_summaries/luecken_false_discoveries_summary.rds')
luecken_false_discoveries_rank %<>%
    filter(!grepl('nebula', method)) %>%
    mutate(
        rank = rank(median, ties.method='first')
    ) %>%
    dplyr::rename(val = median) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Published data',
        # experiment = 'False discoveries',
        experiment = 'Experiment 3: False discoveries',
        tercile = ntile(-rank, 3)
    ) 

snapatac_false_discoveries_rank = readRDS('data/summaries/meta_summaries/snapatac_false_discoveries_summary.rds')
snapatac_false_discoveries_rank %<>%
    mutate(
        rank = rank(median, ties.method='first')
    ) %>%
    dplyr::rename(val = median) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Downsampled data', 
        # experiment = 'False discoveries',
        experiment = 'Experiment 3: False discoveries',
        tercile = ntile(-rank, 3)
    ) 

splatter_false_discoveries_rank = readRDS('data/summaries/meta_summaries/splatter_false_discoveries_summary.rds')
splatter_false_discoveries_rank %<>%
    mutate(
        rank = rank(median, ties.method='first')
    ) %>%
    dplyr::rename(val = median) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Simulations', 
        # experiment = 'False discoveries',
        experiment = 'Experiment 3: False discoveries',
        tercile = ntile(-rank, 3)
    ) 

bulk_mean_expr_bias_rank = readRDS('data/summaries/meta_summaries/bulk_aucc_mean_expr_bias_summary.rds')
bulk_mean_expr_bias_rank %<>%
    mutate(
        rank = rank(median, ties.method='first')
    ) %>%
    dplyr::rename(val = median) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(type = 'mean_expr')

bulk_pct_open_bias_rank = readRDS('data/summaries/meta_summaries/bulk_aucc_pct_open_bias_summary.rds')
bulk_pct_open_bias_rank %<>%
    mutate(
        rank = rank(median, ties.method='first')
    ) %>%
    dplyr::rename(val = median) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(type = 'pct_open') 

bulk_peak_width_bias_rank = readRDS('data/summaries/meta_summaries/bulk_aucc_peak_width_bias_summary.rds')
bulk_peak_width_bias_rank %<>%
    mutate(
        rank = rank(median, ties.method='first')
    ) %>%
    dplyr::rename(val = median) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(type = 'peak_width') 

bulk_bias = rbind(bulk_mean_expr_bias_rank, bulk_pct_open_bias_rank, bulk_peak_width_bias_rank) 
bulk_bias %<>%
    group_by(method) %>%
    summarise(val = mean(rank)) %>%
    ungroup() %>%
    mutate(
        rank = rank(val, ties.method='first')
    ) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Matched bulk', 
        # experiment = 'Biases',
        experiment = 'Experiment 4: Biases',
        tercile = ntile(-rank, 3)
    )

luecken_mean_expr_bias_rank = readRDS('data/summaries/meta_summaries/luecken_false_discoveries_mean_bias_summary.rds')
luecken_mean_expr_bias_rank %<>%
    ungroup() %>%
    filter(decile == 10) %>%
    mutate(
        rank = rank(median_number_of_da_regions, ties.method='first')
    ) %>%
    dplyr::rename(val = median_number_of_da_regions) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(type = 'mean_expr')

luecken_pct_open_bias_rank = readRDS('data/summaries/meta_summaries/luecken_false_discoveries_percentage_open_bias_summary.rds')
luecken_pct_open_bias_rank %<>%
    ungroup() %>%
    filter(decile == 10) %>%
    mutate(
        rank = rank(median_number_of_da_regions, ties.method='first')
    ) %>%
    dplyr::rename(val = median_number_of_da_regions) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(type = 'pct_open')

luecken_peak_width_bias_rank = readRDS('data/summaries/meta_summaries/luecken_false_discoveries_peak_size_bias_summary.rds')
luecken_peak_width_bias_rank %<>%
    ungroup() %>%
    filter(decile == 10) %>%
    mutate(
        rank = rank(median_number_of_da_regions, ties.method='first')
    ) %>%
    dplyr::rename(val = median_number_of_da_regions) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(type = 'peak_width')

luecken_bias = rbind(luecken_mean_expr_bias_rank, luecken_pct_open_bias_rank, luecken_peak_width_bias_rank) 
luecken_bias %<>%
    group_by(method) %>%
    summarise(val = mean(rank)) %>%
    ungroup() %>%
    mutate(
        rank = rank(val, ties.method='first')
    ) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Published data', 
        # experiment = 'Biases',
        experiment = 'Experiment 4: Biases',
        tercile = ntile(-rank, 3)
    ) 

snapatac_mean_expr_bias_rank = readRDS('data/summaries/meta_summaries/snapatac_false_discoveries_mean_bias_summary.rds')
snapatac_mean_expr_bias_rank %<>%
    ungroup() %>%
    filter(decile == 10) %>%
    mutate(
        rank = rank(median_number_of_da_regions, ties.method='first')
    ) %>%
    dplyr::rename(val = median_number_of_da_regions) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(type = 'mean_expr')

snapatac_pct_open_bias_rank = readRDS('data/summaries/meta_summaries/snapatac_false_discoveries_percentage_open_bias_summary.rds')
snapatac_pct_open_bias_rank %<>%
    ungroup() %>%
    filter(decile == 10) %>%
    mutate(
        rank = rank(median_number_of_da_regions, ties.method='first')
    ) %>%
    dplyr::rename(val = median_number_of_da_regions) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(type = 'pct_open')

snapatac_peak_width_bias_rank = readRDS('data/summaries/meta_summaries/snapatac_false_discoveries_peak_size_bias_summary.rds')
snapatac_peak_width_bias_rank %<>%
    ungroup() %>%
    filter(decile == 10) %>%
    mutate(
        rank = rank(median_number_of_da_regions, ties.method='first')
    ) %>%
    dplyr::rename(val = median_number_of_da_regions) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(type = 'peak_width')

snapatac_bias = rbind(snapatac_mean_expr_bias_rank, snapatac_pct_open_bias_rank, snapatac_peak_width_bias_rank) 
snapatac_bias %<>%
    group_by(method) %>%
    summarise(val = mean(rank)) %>%
    ungroup() %>%
    mutate(
        rank = rank(val, ties.method='first')
    ) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Downsampled data', 
        # experiment = 'Biases',
        experiment = 'Experiment 4: Biases',
        tercile = ntile(-rank, 3)
    ) 

splatter_mean_expr_bias_rank = readRDS('data/summaries/meta_summaries/splatter_false_discoveries_mean_bias_summary.rds')
splatter_mean_expr_bias_rank %<>%
    ungroup() %>%
    filter(decile == 10) %>%
    mutate(
        rank = rank(median_number_of_da_regions, ties.method='first')
    ) %>%
    dplyr::rename(val = median_number_of_da_regions) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(type = 'mean_expr')

splatter_pct_open_bias_rank = readRDS('data/summaries/meta_summaries/splatter_false_discoveries_percentage_open_bias_summary.rds')
splatter_pct_open_bias_rank %<>%
    ungroup() %>%
    filter(decile == 10) %>%
    mutate(
        rank = rank(median_number_of_da_regions, ties.method='first')
    ) %>%
    dplyr::rename(val = median_number_of_da_regions) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(type = 'pct_open')

splatter_bias = rbind(splatter_mean_expr_bias_rank, splatter_pct_open_bias_rank) 
splatter_bias %<>%
    group_by(method) %>%
    summarise(val = mean(rank)) %>%
    ungroup() %>%
    mutate(
        rank = rank(val, ties.method='first')
    ) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Simulations', 
        # experiment = 'Biases',
        experiment = 'Experiment 4: Biases',
        tercile = ntile(-rank, 3)
    ) 

bulk_aucc_time = readRDS('data/summaries/meta_summaries/bulk_aucc_time.rds') %>%
    mutate(
        rank = rank(val, ties.method='first')
    ) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Runtime, matched bulk', 
        experiment = 'Experiment 8: Scalability',
        tercile = ntile(-rank, 3)
    )

bulk_aucc_mem = readRDS('data/summaries/meta_summaries/bulk_aucc_mem.rds') %>%
    mutate(
        rank = rank(val, ties.method='first')
    ) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Memory, matched bulk', 
        experiment = 'Experiment 8: Scalability',
        tercile = ntile(-rank, 3)
    )

multiome_aucc_time = readRDS('data/summaries/meta_summaries/multiome_aucc_time.rds') %>%
    mutate(
        rank = rank(val, ties.method='first')
    ) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Runtime, multiome', 
        experiment = 'Experiment 8: Scalability',
        tercile = ntile(-rank, 3)
    )

multiome_aucc_mem = readRDS('data/summaries/meta_summaries/multiome_aucc_mem.rds') %>%
    mutate(
        rank = rank(val, ties.method='first')
    ) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Memory, multiome', 
        experiment = 'Experiment 8: Scalability',
        tercile = ntile(-rank, 3)
    )

dat = rbind(
    bulk_aucc_rank,
    multiome_aucc_rank,
    multiome_gsea_aucc_rank,
    luecken_false_discoveries_rank,
    snapatac_false_discoveries_rank,
    splatter_false_discoveries_rank,
    bulk_bias,
    luecken_bias,
    snapatac_bias,
    splatter_bias,
    bulk_aucc_time,
    multiome_aucc_time,
    bulk_aucc_mem,
    multiome_aucc_mem
)

final_rank = dat %>%
    group_by(method) %>%
    summarise(
        val = mean(tercile)
    ) %>%
    ungroup() %>%
    mutate(
        rank = rank(-val),
        tercile = ntile(-rank, 3),
        # type = '',
        # experiment = 'Summary'
        experiment = 'Summary',
        type = 'Summary'
    ) %>%
    arrange(rank)

dat %<>% rbind(final_rank)
method_order = rev(unique(final_rank$method))
dat$method = factor(dat$method, levels=method_order)

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
# dat$experiment = factor(dat$experiment, levels = c('Matched bulk concordance', 'Multiome concordance', 'False discoveries', 'Biases', 'Summary'))
dat$experiment = factor(dat$experiment, levels = c(
    'Experiment 1: Bulk concordance',
    'Experiment 2: Multiome concordance',
    'Experiment 3: False discoveries',
    'Experiment 4: Biases', 
    'Experiment 8: Scalability',
    'Summary'))


[1] "#000000" "#0E0B0C" "#1C1717" "#292223" "#362D2E" "#423839" "#4E4244" "#594D4F" "#645759" "#6E6063" "#786A6C" "#827376" "#8B7C7F" "#938588"
[15] "#9B8E91" "#A39699" "#AA9EA1" "#B1A6A9" "#B8AEB1" "#BEB5B8" "#C4BCBF" "#C9C3C6" "#CECACC" "#D3D0D2" "#D7D6D8" "#DBDCDE" "#DFE1E3" "#E3E7E8"
[29] "#E6ECED" "#E9F0F2" "#EBF5F6" "#EEF9FA" "#F0FCFD" "#F2FFFF" "#F4FFFF" "#F5FFFF" "#F7FFFF" "#F8FFFF" "#F9FFFF" "#FAFFFF" "#FBFFFF" "#FBFFFF"
[43] "#FCFFFF" "#FCFFFF" "#FDFFFF" "#FDFFFF" "#FEFFFF" "#FEFFFF" "#FEFFFF" "#FEFFFF" "#FFFDFD" "#FFF8F9" "#FFF4F4" "#FFEEEF" "#FFE9EA" "#FFE3E4"
[57] "#FFDDDE" "#FFD6D8" "#FFCFD2" "#FFC8CB" "#FFC1C4" "#FFB9BD" "#FFB2B6" "#FFAAAE" "#FFA2A7" "#FF9A9F" "#FF9298" "#FF8A90" "#FF8288" "#FF7A81"
[71] "#FE7279" "#FD6A72" "#FC626A" "#FB5B63" "#FA535C" "#F84C55" "#F7454F" "#F53F48" "#F33842" "#F1323C" "#EF2C37" "#ED2732" "#EA222D" "#E71D28"
[85] "#E51924" "#E11620" "#DE121D" "#DB0F1A" "#D70D18" "#D30B16" "#CF0A14" "#CB0913" "#C70913" "#C20A13" "#BD0B14" "#B90C15" "#B30F17" "#AE1219"
[99] "#A8161D" "#A31B21"
# pal = brewer.pal(3, 'Blues')
pal = colorRampPalette(c(brewer.pal(3, 'Reds')[1], "#B30F17"))(3)
pal = colorRampPalette(c("#FF9298", "#DB0F1A"))(3)

experiments = levels(dat$experiment)
for (experiment in experiments) {
    if (experiment == first(experiments)) out_plots = list()
    dat1 = dat %>%
        filter(experiment == !!experiment)
    p1_1 = dat1 %>%
        ggplot(aes(x = method, y = type, fill=tercile)) +
        # facet_grid(experiment ~ ., scales = 'free', space = 'free') +
        geom_tile(color = 'white') +
        scale_x_discrete(expand = c(0, 0), labels = scales::parse_format()) +
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
    if (experiment == last(experiments)) {
        p1_1 = p1_1 + 
            theme(
                axis.text.y = element_blank()
            )
    }
    out_plots[[length(out_plots)+1]] = p1_1
}

p1 = wrap_plots(out_plots, ncol=1)
p1
ggsave('fig/Fig9/full.pdf', p1, width=9, height=15, units='cm')

