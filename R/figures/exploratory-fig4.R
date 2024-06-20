# Plot a summary figure with performance across all experiments.
setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(effsize)
library(upstartr)
library(ggpubr)
source("R/theme.R")

dat = readRDS('data/metadata/enhancer-promoter-peak_stats.rds') %>%
    type_convert() %>%
    mutate(
        peak_type = fct_recode(
            peak_type,
            'Enhancer'='enhancer',
            'Promoter'='promoter'
        ),
        pct_open = pct_open*100
    )

dat0 = dat %>%
    filter(data_type == 'sc')

color_pal = c(
    'Enhancer'="#FFDCA4",
    'Promoter'="#887AA1"
)

labs = dat0 %>%
    group_by(peak_type) %>%
    summarise(
        mean_lab = round(median(mean), 3),
        pct_open_lab = round(median(pct_open), 3),
        peak_size_lab = round(median(peak_size), 3)
    ) %>%
    ungroup()

p1 = dat0 %>%
    mutate(
        peak_type = factor(peak_type, levels=c("Enhancer", "Promoter"))
    ) %>%
    ggplot(aes(x=peak_type, y=mean, color=peak_type, fill=peak_type)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
    geom_label(dat = labs,
               aes(label = mean_lab, y= 0.18), color = 'black',
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0.5, vjust = 0,
               show.legend = FALSE) +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    # geom_signif(
    #     comparisons = list(c("Enhancer", "Promoter")),
    #     map_signif_level = TRUE,
    #     color = 'black',
    #     tip_length = 0,
    #     size = 0.3,
    #     textsize = 3
    # ) +
    boxed_theme() +
    scale_y_continuous('Read depth', limits=c(0, 0.2), breaks=seq(0, 0.18, 0.06)) +
    scale_x_discrete() +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 1.6,
        legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title = element_text(size=5)
    )
p1

p2 = dat0 %>%
    mutate(
        peak_type = factor(peak_type, levels=c("Enhancer", "Promoter"))
    ) %>%
    ggplot(aes(x=peak_type, y=pct_open, color=peak_type, fill=peak_type)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
    geom_label(dat = labs,
               aes(label = pct_open_lab, y= 13.2), color = 'black',
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0.5, vjust = 0,
               show.legend = FALSE) +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    # geom_signif(
    #     comparisons = list(c("Enhancer", "Promoter")),
    #     map_signif_level = TRUE,
    #     color = 'black',
    #     tip_length = 0,
    #     size = 0.3,
    #     textsize = 3,
    #     y_position = 0
    # ) +
    boxed_theme() +
    scale_y_continuous('% open', limits=c(0, 14.3), breaks=seq(0, 14, 4)) +
    scale_x_discrete() +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 1.6,
        legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title = element_text(size=5)
    )
p2

p3 = dat0 %>%
    mutate(
        peak_type = factor(peak_type, levels=c("Enhancer", "Promoter"))
    ) %>%
    ggplot(aes(x=peak_type, y=peak_size, color=peak_type, fill=peak_type)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
    geom_label(dat = labs,
               aes(label = peak_size_lab, y= 985), color = 'black',
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0.5, vjust = 0,
               show.legend = FALSE) +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    # geom_signif(
    #     comparisons = list(c("Enhancer", "Promoter")),
    #     map_signif_level = TRUE,
    #     color = 'black',
    #     tip_length = 0,
    #     size = 0.3,
    #     textsize = 3,
    #     y_position = 0
    # ) +
    boxed_theme() +
    scale_y_continuous('Peak width (bp)', limits=c(0, 1100)) +
    scale_x_discrete() +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = 1.6,
        legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title = element_text(size=5)
    )
p3

row = p1 | p2 | p3
ggsave('fig/Exp_Fig4/row.pdf', row, width=9, height=5, units='cm')

