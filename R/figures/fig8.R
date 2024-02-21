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
dat1$val = dat1$time
labs = dat1 %>% 
    group_by(method, color) %>%
    summarise(mean = mean(val), n = n(), median = median(val), val = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
saveRDS(labs, 'data/summaries/meta_summaries/bulk_aucc_time.rds')

p1 = dat1 %>%
    ggplot(aes(x = reorder(method, val, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = val), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    geom_text(data = labs, aes(y = -70, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('Runtime (min)', limits=c(-80, 170), breaks=seq(0,200,50)) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          aspect.ratio = 1.4)
p1

dat1$val = dat1$mem
labs = dat1 %>% 
    group_by(method, color) %>%
    summarise(mean = mean(val), n = n(), median = median(val), val = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
saveRDS(labs, 'data/summaries/meta_summaries/bulk_aucc_mem.rds')

p2 = dat1 %>%
    ggplot(aes(x = reorder(method, val, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = val), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    # geom_jitter(aes(y = val), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -1.75, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('Memory (GB)', limits=c(-2, 6), breaks=c(0:6)) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          aspect.ratio = 1.4)
p2

# Experiment 2
dat2 = dat0 %>% filter(standard == 'silver', exp == 'scATAC') %>%
    filter(!grepl('ArchR', method), !grepl('glmGamPoi', method))
dat2 %>% dplyr::select(method) %>% table()

# check time
dat2$val = dat2$time
labs = dat2 %>% 
    group_by(method, color) %>%
    summarise(mean = mean(val), n = n(), median = median(val), val = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
saveRDS(labs, 'data/summaries/meta_summaries/multiome_aucc_time.rds')

p3 = dat2 %>%
    ggplot(aes(x = reorder(method, val, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = val), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    # geom_jitter(aes(y = val), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -230, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('Runtime (min)', limits=c(-250, 620), breaks=seq(0,600,200)) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          aspect.ratio = 1.4)
p3

# check time
dat2$val = dat2$mem
labs = dat2 %>% 
    group_by(method, color) %>%
    summarise(mean = mean(val), n = n(), median = median(val), val = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
saveRDS(labs, 'data/summaries/meta_summaries/mutliome_aucc_mem.rds')

p4 = dat2 %>%
    ggplot(aes(x = reorder(method, val, stats::median), 
               color = method, fill = method)) +
    geom_boxplot(aes(y = val), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    # geom_jitter(aes(y = val), width = 0.2, height = 0, shape = 1, size = 0.5,
    # stroke = 0.3) +
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    geom_text(data = labs, aes(y = -7, label = label), size = 1.75, hjust = 0,
              color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous('Memory (GB)', limits=c(-8, 20), breaks=seq(0,20,5)) +
    coord_flip() +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.title.y = element_blank(),
          legend.position = 'none',
          aspect.ratio = 1.4)
p4

p5 = p1 | p2 | p3 | p4
ggsave("fig/Fig8/row-1.pdf", p5, width = 18, height = 9, units = "cm", useDingbats = FALSE)
 
alt_colors = colorRampPalette(c("#B8E600", "#218716"))(3) %>% setNames(c("t~test", "Wilcoxon~rank-sum~test", "Negative~binomial"))

# alternative implementations
dat3 = dat0 %>%
    filter(method %in% c("'ArchR::t test'", "t~test")) %>%
    mutate(
        method = fct_recode(method,
                            'Signac' = 't~test',
                            'ArchR' = "'ArchR::t test'"
        ) %>% as.character()
    )
dat3 %>% dplyr::select(method) %>% table()

t.test(
    dat3 %>% filter(method == 'ArchR') %>% arrange(cell_type) %>% pull(time),
    dat3 %>% filter(method == 'Signac') %>% arrange(cell_type) %>% pull(time),
    paired = T
)

t.test(
    dat3 %>% filter(method == 'ArchR') %>% arrange(cell_type) %>% pull(mem),
    dat3 %>% filter(method == 'Signac') %>% arrange(cell_type) %>% pull(mem),
    paired = T
)

labs = dat3 %>% 
    group_by(method) %>% 
    summarise(
        text_val = round(mean(time),3), time = -100) %>%
    ungroup()

pal2 = c('Signac' = unname(da_analysis_colors['t~test']),
         'ArchR' = unname(alt_colors['t~test']))

p6 = dat3 %>%
    mutate(
        method = factor(method, levels=c("Signac", "ArchR"))
    ) %>%
    ggplot(aes(x=method, y=time, color=method, fill=method)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
    geom_label(dat = labs,
               aes(label = text_val), color = 'black',
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0.5, vjust = 0,
               show.legend = FALSE) +
    scale_fill_manual(values = pal2) +
    scale_color_manual(values = pal2) +
    geom_signif(
        comparisons = list(c("Signac", "ArchR")),
        map_signif_level = TRUE,
        color = 'black',
        tip_length = 0,
        size = 0.3,
        textsize = 3,
        y_position = 500
    ) +
    ggtitle('t test') +
    boxed_theme() +
    scale_y_continuous('Runtime (min)', limits=c(-100, 700), breaks=seq(0,600,200)) +
    scale_x_discrete(labels = scales::parse_format()) +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 1.6,
        legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title = element_text(size=5)
    )
p6
ggsave('fig/Fig8/alt-t-test-time.pdf', p6, width=3, height=8, units='cm')

labs = dat3 %>% 
    group_by(method) %>% 
    summarise(
        text_val = round(mean(mem),3), mem = -2) %>%
    ungroup()

p7 = dat3 %>%
    mutate(
        method = factor(method, levels=c("Signac", "ArchR"))
    ) %>%
    ggplot(aes(x=method, y=mem, color=method, fill=method)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
    geom_label(dat = labs,
               aes(label = text_val), color = 'black',
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0.5, vjust = 0,
               show.legend = FALSE) +
    scale_fill_manual(values = pal2) +
    scale_color_manual(values = pal2) +
    geom_signif(
        comparisons = list(c("Signac", "ArchR")),
        map_signif_level = TRUE,
        color = 'black',
        tip_length = 0,
        size = 0.3,
        textsize = 3,
        y_position = 12
    ) +
    ggtitle('t test') +
    boxed_theme() +
    scale_y_continuous('Memory (GB)', limits=c(-2, 15), breaks=seq(0,10,5)) +
    scale_x_discrete(labels = scales::parse_format()) +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 1.6,
        legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title = element_text(size=5)
    )
p7
ggsave('fig/Fig8/alt-t-test-mem.pdf', p7, width=3, height=8, units='cm')

dat4 = dat0 %>%
    filter(method %in% c("'ArchR::Wilcoxon rank-sum test'", "Wilcoxon~rank-sum~test")) %>%
    mutate(
        method = fct_recode(method,
                            'Signac' = 'Wilcoxon~rank-sum~test',
                            'ArchR' = "'ArchR::Wilcoxon rank-sum test'"
        ) %>% as.character()
    )
dat4 %>% dplyr::select(method) %>% table()

t.test(
    dat4 %>% filter(method == 'ArchR') %>% arrange(cell_type) %>% pull(time),
    dat4 %>% filter(method == 'Signac') %>% arrange(cell_type) %>% pull(time),
    paired = T
)

t.test(
    dat4 %>% filter(method == 'ArchR') %>% arrange(cell_type) %>% pull(mem),
    dat4 %>% filter(method == 'Signac') %>% arrange(cell_type) %>% pull(mem),
    paired = T
)

labs = dat4 %>% 
    group_by(method) %>% 
    summarise(
        text_val = round(mean(time),3), time = -0.1) %>%
    ungroup()

pal2 = c('Signac' = unname(da_analysis_colors['Wilcoxon~rank-sum~test']),
         'ArchR' = unname(alt_colors['Wilcoxon~rank-sum~test']))

p14 = dat4 %>%
    mutate(
        method = factor(method, levels=c("Signac", "ArchR"))
    ) %>%
    ggplot(aes(x=method, y=time, color=method, fill=method)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
    geom_label(dat = labs,
               aes(label = text_val), color = 'black',
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0.5, vjust = 0,
               show.legend = FALSE) +
    scale_fill_manual(values = pal2) +
    scale_color_manual(values = pal2) +
    geom_signif(
        comparisons = list(c("Signac", "ArchR")),
        map_signif_level = TRUE,
        color = 'black',
        tip_length = 0,
        size = 0.3,
        textsize = 3,
        y_position = 0.55
    ) +
    ggtitle('Wilcoxon rank-sum test') + 
    boxed_theme() +
    scale_y_continuous('Runtime (min)', limits=c(-0.1, 0.7), breaks=seq(0,0.6,0.2)) +
    scale_x_discrete(labels = scales::parse_format()) +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 1.6,
        legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title = element_text(size=5)
    )
p14
ggsave('fig/Fig8/alt-wilcoxon-test-time.pdf', p14, width=3, height=8, units='cm')

labs = dat4 %>% 
    group_by(method) %>% 
    summarise(
        text_val = round(mean(mem),3), mem = -2) %>%
    ungroup()

p15 = dat4 %>%
    mutate(
        method = factor(method, levels=c("Signac", "ArchR"))
    ) %>%
    ggplot(aes(x=method, y=mem, color=method, fill=method)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
    geom_label(dat = labs,
               aes(label = text_val), color = 'black',
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0.5, vjust = 0,
               show.legend = FALSE) +
    scale_fill_manual(values = pal2) +
    scale_color_manual(values = pal2) +
    geom_signif(
        comparisons = list(c("Signac", "ArchR")),
        map_signif_level = TRUE,
        color = 'black',
        tip_length = 0,
        size = 0.3,
        textsize = 3
    ) +
    boxed_theme() +
    ggtitle('Wilcoxon rank-sum test') +
    scale_y_continuous('Memory (GB)', limits=c(-2, 17), breaks=seq(0,15,5)) +
    scale_x_discrete(labels = scales::parse_format()) +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 1.6,
        legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title = element_text(size=5)
    )
p15
ggsave('fig/Fig8/alt-wilcoxon-test-mem.pdf', p15, width=3, height=8, units='cm')

dat5 = dat0 %>%
    filter(
        method %in% c("Negative~binomial", "'glmGamPoi::Negative binomial'"),
        paper != 'Pliner_2018'
    ) %>%
    mutate(
        method = fct_recode(method,
                            'Signac' = 'Negative~binomial',
                            'glmGamPoi' = "'glmGamPoi::Negative binomial'"
        ) %>% as.character()
    )
dat5 %>% dplyr::select(method) %>% table()
t.test(
    dat5 %>% filter(method == 'glmGamPoi') %>% arrange(cell_type) %>% pull(time),
    dat5 %>% filter(method == 'Signac') %>% arrange(cell_type) %>% pull(time),
    paired = T
)

t.test(
    dat5 %>% filter(method == 'glmGamPoi') %>% arrange(cell_type) %>% pull(mem),
    dat5 %>% filter(method == 'Signac') %>% arrange(cell_type) %>% pull(mem),
    paired = T
)



labs = dat5 %>% 
    group_by(method) %>% 
    summarise(
        text_val = round(mean(time),3), time = -130) %>%
    ungroup()

pal2 = c('Signac' = unname(da_analysis_colors['Negative~binomial']),
         'glmGamPoi' = unname(alt_colors['Negative~binomial']))

p16 = dat5 %>%
    mutate(
        method = factor(method, levels=c("Signac", "glmGamPoi"))
    ) %>%
    ggplot(aes(x=method, y=time, color=method, fill=method)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
    geom_label(dat = labs,
               aes(label = text_val), color = 'black',
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0.5, vjust = 0,
               show.legend = FALSE) +
    scale_fill_manual(values = pal2) +
    scale_color_manual(values = pal2) +
    geom_signif(
        comparisons = list(c("Signac", "glmGamPoi")),
        map_signif_level = TRUE,
        color = 'black',
        tip_length = 0,
        size = 0.3,
        textsize = 3,
        y_position = 600
    ) +
    boxed_theme() +
    ggtitle('Negative binomial') +
    scale_y_continuous('Runtime (min)', limits=c(-130, 800), breaks=seq(0,600,200)) +
    scale_x_discrete(labels = scales::parse_format()) +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 1.6,
        legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title = element_text(size=5)
    )
p16
ggsave('fig/Fig8/alt-nb-test-time.pdf', p16, width=3, height=8, units='cm')

labs = dat5 %>% 
    group_by(method) %>% 
    summarise(
        text_val = round(mean(mem),3), mem = -5) %>%
    ungroup()

p17 = dat5 %>%
    mutate(
        method = factor(method, levels=c("Signac", "glmGamPoi"))
    ) %>%
    ggplot(aes(x=method, y=mem, color=method, fill=method)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
    geom_label(dat = labs,
               aes(label = text_val), color = 'black',
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0.5, vjust = 0,
               show.legend = FALSE) +
    scale_fill_manual(values = pal2) +
    scale_color_manual(values = pal2) +
    geom_signif(
        comparisons = list(c("Signac", "glmGamPoi")),
        map_signif_level = TRUE,
        color = 'black',
        tip_length = 0,
        size = 0.3,
        textsize = 3,
        y_position = 25
    ) +
    boxed_theme() +
    ggtitle('Negative binomial') +
    scale_y_continuous('Memory (GB)', limits=c(-5, 35), breaks=seq(0,30,10)) +
    scale_x_discrete(labels = scales::parse_format()) +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 1.6,
        legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title = element_text(size=5)
    )
p17
ggsave('fig/Fig8/alt-nb-test-mem.pdf', p17, width=3, height=8, units='cm')

bottom_row = p12 | p14 + theme(axis.title.y = element_blank()) | p16 + theme(axis.title.y = element_blank()) | p13 | p15 + theme(axis.title.y = element_blank()) | p17 + theme(axis.title.y = element_blank())
ggsave('fig/Fig8/bottom-row.pdf', bottom_row, width=15, height=9, units='cm')

