setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

### PANEL A ###
bulk_aucc_rank1 = readRDS('data/summaries/meta_summaries/bulk_aucc-leave_out_bulk-summary.rds')
bulk_aucc_rank1 %<>%
    group_by(bulk_da_method, bulk_da_type) %>%
    mutate(
        rank = rank(-aucc, ties.method='first'),
        method = ifelse(method == 'LR_peaks', 'LR[peaks]', method)
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(bulk_da_method, bulk_da_type, method, val, rank) %>%
    arrange(bulk_da_method, bulk_da_type, rank) %>%
    mutate(
        type = paste0('Matched bulk (leave out: ', bulk_da_method, '-', bulk_da_type, ')'), 
        tercile = ntile(-rank, 3),
    ) %>%
    ungroup() %>%
    dplyr::select(-bulk_da_method, -bulk_da_type)

bulk_aucc_rank2 = readRDS('data/summaries/meta_summaries/bulk_aucc-k=1000-summary.rds')
bulk_aucc_rank2 %<>%
    mutate(
        rank = rank(-aucc, ties.method='first'),
        method = ifelse(method == 'LR_peaks', 'LR[peaks]', method)
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Matched bulk (k = 1000)', 
        tercile = ntile(-rank, 3)
    )

bulk_aucc_rank3 = readRDS('data/summaries/meta_summaries/bulk_aucc-filter=5-summary.rds')
bulk_aucc_rank3 %<>%
    mutate(
        rank = rank(-aucc, ties.method='first'),
        method = ifelse(method == 'LR_peaks', 'LR[peaks]', method)
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Matched bulk  (5% cell filter)', 
        tercile = ntile(-rank, 3)
    )

bulk_aucc_rank4 = readRDS('data/summaries/meta_summaries/bulk_aucc-pseudobulk_peaks-summary.rds')
bulk_aucc_rank4 %<>%
    mutate(
        rank = rank(-aucc, ties.method='first'),
        method = ifelse(method == 'LR_peaks', 'LR[peaks]', method)
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Matched bulk  (Pseudobulk peak calling)', 
        tercile = ntile(-rank, 3)
    )


bulk_aucc_rank5 = readRDS('data/summaries/meta_summaries/bulk_aucc-noisy_peaks-summary.rds')
bulk_aucc_rank5 %<>%
    mutate(
        rank = rank(-aucc, ties.method='first'),
        method = ifelse(method == 'LR_peaks', 'LR[peaks]', method)
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Matched bulk  (Spurious peaks added)', 
        tercile = ntile(-rank, 3)
    )

bulk_aucc_rank6 = readRDS('data/summaries/meta_summaries/bulk_aucc-enhancer_promoter-summary.rds')
bulk_aucc_rank6 %<>%
    group_by(ann) %>%
    mutate(
        rank = rank(-aucc, ties.method='first'),
        method = ifelse(method == 'LR_peaks', 'LR[peaks]', method)
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(ann, method, val, rank) %>%
    arrange(ann, rank) %>%
    mutate(
        type = paste0('Matched bulk: Bulk concordance (', str_to_title(ann), ' regions only)'),
        tercile = ntile(-rank, 3)
    ) %>%
    ungroup() %>%
    dplyr::select(-ann)

bulk_aucc_rank = rbind(
    bulk_aucc_rank1,
    bulk_aucc_rank2,
    bulk_aucc_rank3,
    bulk_aucc_rank4,
    bulk_aucc_rank5,
    bulk_aucc_rank6
) %>%
    mutate(
        experiment = 'Experiment 1: Bulk concordance'
    )

multiome_aucc_rank1 = readRDS('data/summaries/meta_summaries/multiome_aucc-bidrectional_promoter-summary.rds')
multiome_aucc_rank1 %<>%
    mutate(
        rank = rank(-aucc, ties.method='first')
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Multiome (Removal of overlapping promoter regions)',
        tercile = ntile(-rank, 3)
    )


multiome_aucc_rank2 = readRDS('data/summaries/meta_summaries/multiome_aucc-leave_out_bulk-summary.rds')
multiome_aucc_rank2 %<>%
    group_by(scrna_da_method, scrna_da_type) %>%
    mutate(
        rank = rank(-aucc, ties.method='first')
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(scrna_da_method, scrna_da_type, method, val, rank) %>%
    arrange(scrna_da_method, scrna_da_type, rank) %>%
    mutate(
        type = paste0('Multiome (leave out: ', scrna_da_method, '-', scrna_da_type, ')'),
        tercile = ntile(-rank, 3)
    ) %>%
    ungroup() %>% 
    dplyr::select(-scrna_da_method, -scrna_da_type)

multiome_aucc_rank3 = readRDS('data/summaries/meta_summaries/multiome_aucc-k=100-summary.rds')
multiome_aucc_rank3 %<>%
    mutate(
        rank = rank(-aucc, ties.method='first')
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Multiome (k = 100)',
        tercile = ntile(-rank, 3)
    ) 

multiome_aucc_rank4 = readRDS('data/summaries/meta_summaries/multiome_aucc-k=1000-summary.rds')
multiome_aucc_rank4 %<>%
    mutate(
        rank = rank(-aucc, ties.method='first')
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Multiome (k = 1000)',
        tercile = ntile(-rank, 3)
    ) 

multiome_aucc_rank5 = readRDS('data/summaries/meta_summaries/multiome_aucc-filter=1-summary.rds')
multiome_aucc_rank5 %<>%
    mutate(
        rank = rank(-aucc, ties.method='first')
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Multiome (1% cell filter)',
        tercile = ntile(-rank, 3)
    )

multiome_aucc_rank6 = readRDS('data/summaries/meta_summaries/multiome_aucc-scrna_terciles-summary.rds')
multiome_aucc_rank6 %<>%
    mutate(
        rank = rank(-aucc, ties.method='first')
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Multiome (Removal of lowly expressed genes)',
        tercile = ntile(-rank, 3)
    ) 

multiome_gsea_aucc_rank1 = readRDS('data/summaries/meta_summaries/multiome_gsea_aucc-bidirectional_promoter-summary.rds')
multiome_gsea_aucc_rank1 %<>%
    mutate(
        rank = rank(-aucc, ties.method='first')
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Multiome GO  (Removal of overlapping promoter regions)', 
        tercile = ntile(-rank, 3)
    )

multiome_gsea_aucc_rank2 = readRDS('data/summaries/meta_summaries/multiome_gsea_aucc-leave_bulk_out-summary.rds')
multiome_gsea_aucc_rank2 %<>%
    group_by(scrna_da_method, scrna_da_type) %>%
    mutate(
        rank = rank(-aucc, ties.method='first')
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(scrna_da_method, scrna_da_type, method, val, rank) %>%
    arrange(scrna_da_method, scrna_da_type, rank) %>%
    mutate(
        type = paste0('Multiome GO (leave out: ', scrna_da_method, '-', scrna_da_type, ')'), 
        tercile = ntile(-rank, 3)
    ) %>%
    ungroup() %>% 
    dplyr::select(-scrna_da_method, -scrna_da_type)

multiome_gsea_aucc_rank3 = readRDS('data/summaries/meta_summaries/multiome_gsea_aucc-k=50-summary.rds')
multiome_gsea_aucc_rank3 %<>%
    mutate(
        rank = rank(-aucc, ties.method='first')
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Multiome GO (k = 50)', 
        tercile = ntile(-rank, 3)
    ) 

multiome_gsea_aucc_rank4 = readRDS('data/summaries/meta_summaries/multiome_gsea_aucc-k=500-summary.rds')
multiome_gsea_aucc_rank4 %<>%
    mutate(
        rank = rank(-aucc, ties.method='first')
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Multiome GO (k = 500)',
        tercile = ntile(-rank, 3)
    ) 

multiome_gsea_aucc_rank5 = readRDS('data/summaries/meta_summaries/multiome_gsea_aucc-filter=1-summary.rds')
multiome_gsea_aucc_rank5 %<>%
    mutate(
        rank = rank(-aucc, ties.method='first')
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Multiome GO (1% cell filter)', 
        tercile = ntile(-rank, 3)
    )

multiome_gsea_aucc_rank6 = readRDS('data/summaries/meta_summaries/multiome_gsea_aucc-scrna_terciles-summary.rds')
multiome_gsea_aucc_rank6 %<>%
    mutate(
        rank = rank(-aucc, ties.method='first')
    ) %>%
    dplyr::rename(val = aucc) %>%
    dplyr::select(method, val, rank) %>%
    arrange(rank) %>%
    mutate(
        type = 'Multiome GO (Removal of lowly expressed genes)', 
        tercile = ntile(-rank, 3)
    ) 

multiome_aucc_rank = rbind(
    multiome_aucc_rank1,
    multiome_aucc_rank2,
    multiome_aucc_rank3,
    multiome_aucc_rank4,
    multiome_aucc_rank5,
    multiome_aucc_rank6,
    multiome_gsea_aucc_rank1,
    multiome_gsea_aucc_rank2,
    multiome_gsea_aucc_rank3,
    multiome_gsea_aucc_rank4,
    multiome_gsea_aucc_rank5,
    multiome_gsea_aucc_rank6
) %>%
    mutate(
        experiment = 'Experiment 2: Multiome concordance'
    )

dat = rbind(
    bulk_aucc_rank,
    multiome_aucc_rank
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

# dat %<>% rbind(final_rank)
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
    'Experiment 2: Multiome concordance')
)

# [1] "#000000" "#0E0B0C" "#1C1717" "#292223" "#362D2E" "#423839" "#4E4244" "#594D4F" "#645759" "#6E6063" "#786A6C" "#827376" "#8B7C7F" "#938588"
# [15] "#9B8E91" "#A39699" "#AA9EA1" "#B1A6A9" "#B8AEB1" "#BEB5B8" "#C4BCBF" "#C9C3C6" "#CECACC" "#D3D0D2" "#D7D6D8" "#DBDCDE" "#DFE1E3" "#E3E7E8"
# [29] "#E6ECED" "#E9F0F2" "#EBF5F6" "#EEF9FA" "#F0FCFD" "#F2FFFF" "#F4FFFF" "#F5FFFF" "#F7FFFF" "#F8FFFF" "#F9FFFF" "#FAFFFF" "#FBFFFF" "#FBFFFF"
# [43] "#FCFFFF" "#FCFFFF" "#FDFFFF" "#FDFFFF" "#FEFFFF" "#FEFFFF" "#FEFFFF" "#FEFFFF" "#FFFDFD" "#FFF8F9" "#FFF4F4" "#FFEEEF" "#FFE9EA" "#FFE3E4"
# [57] "#FFDDDE" "#FFD6D8" "#FFCFD2" "#FFC8CB" "#FFC1C4" "#FFB9BD" "#FFB2B6" "#FFAAAE" "#FFA2A7" "#FF9A9F" "#FF9298" "#FF8A90" "#FF8288" "#FF7A81"
# [71] "#FE7279" "#FD6A72" "#FC626A" "#FB5B63" "#FA535C" "#F84C55" "#F7454F" "#F53F48" "#F33842" "#F1323C" "#EF2C37" "#ED2732" "#EA222D" "#E71D28"
# [85] "#E51924" "#E11620" "#DE121D" "#DB0F1A" "#D70D18" "#D30B16" "#CF0A14" "#CB0913" "#C70913" "#C20A13" "#BD0B14" "#B90C15" "#B30F17" "#AE1219"
# [99] "#A8161D" "#A31B21"
# pal = brewer.pal(3, 'Blues')
pal = colorRampPalette(c(brewer.pal(3, 'Reds')[1], "#B30F17"))(3)
# pal = colorRampPalette(c("#FF9298", "#DB0F1A"))(3)

experiments = levels(dat$experiment)
table(dat$tercile)

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
    out_plots[[length(out_plots)+1]] = p1_1
}

p1 = wrap_plots(out_plots, ncol=1)
p1
ggsave('fig/Exp_Fig5/sensitivity-terciles.pdf', p1, width=15, height=20, units='cm')

