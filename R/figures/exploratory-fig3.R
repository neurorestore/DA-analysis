# Plot a summary figure with performance across all experiments.
setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(effsize)
library(upstartr)
source("R/theme.R")

remove_quantile_per_method = function(df, val_col, group_cols, quantile_range=c(0.01, 0.99)) {
    df$idx = 1:nrow(df)
    group_grid = data.frame()
    for (group_col in group_cols) {
        if (nrow(group_grid) == 0) {
            group_grid = df[,group_col]
        } else {
            group_grid %<>% tidyr::crossing(., df[,group_col])    
        }
    }
    
    new_df = data.frame()
    for (i in 1:nrow(group_grid)) {
        tmp_row = group_grid[i,]
        tmp_df = tmp_row %>% left_join(df)
        quan_range = quantile(unlist(tmp_df[,val_col]), quantile_range)
        tmp_df %<>% mutate(
            filter_now = ifelse(get(val_col) < quan_range[2] & get(val_col) > quan_range[1], T, F)
        )
        # if remove all -- probably means all 0
        if (sum(tmp_df$filter_now) == 0) tmp_df$filter_now = T
        new_df %<>% rbind(tmp_df %>% dplyr::select(filter_now, idx))
    }
    new_df = df %>%
        left_join(new_df) %>%
        filter(filter_now) %>%
        dplyr::select(-filter_now, -idx)
    return(new_df)
}


###############################################################################-
## bulk logFC plot ####
###############################################################################-

# load data
dat = readRDS("data/summaries/false_discoveries/false_discoveries_logFC_absolute.rds") %>% 
    type_convert() %>%
    filter(is.na(latent_vars)) %>%
    filter(da_family != 'mixedmodel',
           pseudo_repl == FALSE)
# remove fixed columns
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

dat %<>% 
    mutate(expr = replace(percentile, is.na(percentile), '0.0')) %>% 
    mutate(expr_set = ifelse(grepl('pseudo', expr), 'Random', 'logFC'),
           expr = gsub("_pseudo", "", expr))

# re-code x-values and colors
dat %<>%
    mutate(method = paste0(da_method, '-', da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', da_family) %>% 
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
               as.character()
           ) %>% 
        # drop limma
        filter(!grepl('limma', method))

dat %<>% 
    mutate(logFC = gsub("^.*_", "", expr) %>% 
               replace(. == 'all', 0),
           logFC = case_when(
               logFC == 0 ~ 0,
               logFC == 0.5 ~ 1,
               logFC == 1 ~ 2,
               logFC > 1 ~ 3
           )) %>%
    mutate(logFC = as.factor(logFC))

dat$number_of_da_regions[is.na(dat$number_of_da_regions)] = 0

table(dat$logFC)
dplyr::count(dat, method)
labs = dat %>% 
    group_by(logFC, expr_set, method, color) %>%
    summarise(mean = mean(number_of_da_regions),
              n = n(), 
              median = median(number_of_da_regions), 
              number_of_da_regions = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
pal = c('logFC' = 'black',
        'Random' = 'grey75') 
p1 = dat %>%
    filter(logFC != 0 & logFC != 1) %>%
    # remove_quantile_per_method(., group_cols = c('method', 'logFC'), val_col='number_of_da_regions') %>%
    ggplot(aes(x = logFC, y = number_of_da_regions,
               color = expr_set, fill = expr_set)) +
    facet_wrap(~ method, scales = 'free', ncol = 7, labeller = label_parsed) +
    geom_boxplot(aes(y = number_of_da_regions), outlier.shape = NA, alpha = 0.4,
                 width = 0.7, size = 0.3) + 
    scale_color_manual('Filter by:', values = pal) +
    scale_fill_manual('Filter by:', values = pal) +
    scale_x_reordered(expression('Absolute fold-change filter')) +
    scale_y_continuous('# of DA peaks', breaks=pretty_breaks(n=3)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(legend.position = 'top',
          legend.key.width = unit(0.4, 'lines'),
          aspect.ratio = 0.7,
          # axis.text.x = element_text(angle=45, hjust=1)
          )
p1
ggsave("fig/Exp_Fig3/luecken-logFC-filter-overview.pdf", p1,
       width = 18, height = 15, units = "cm", useDingbats = FALSE)

enr = dat %>% 
    filter(logFC != 0 & logFC != 1) %>%
    group_by_at(vars(-expr_set, -percentile,
                     -number_of_regions, -number_of_da_regions,
                     -number_of_nominal_da_regions)) %>% 
    summarise(delta = ifelse(n() == 2, 
                             number_of_da_regions[expr_set == 'logFC'] - number_of_da_regions[expr_set == 'Random'],
                             number_of_da_regions),
              nn = n()) %>% 
    ungroup() 
pal = da_analysis_colors
p2 = enr %>%
    filter(logFC != 0) %>% 
    ggplot(aes(x = logFC, y = delta,
               color = method, fill = method)) +
    facet_wrap(~ method, scales = 'free', ncol = 7, labeller = label_parsed) +
    geom_boxplot(aes(y = delta), outlier.shape = NA, alpha = 0.4, width = 0.6,
                 size = 0.4) + 
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    # geom_text(data = labs, aes(y = -0.12, label = label), size = 1.75, hjust = 0,
    #           color = 'black') +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    scale_x_reordered(expression('Absolute fold-change filter')) +
    scale_y_continuous(expression(Delta~'# of DA peaks'), breaks=pretty_breaks(3)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(legend.position = 'none',
          aspect.ratio = 0.8,
          # axis.text.x = element_text(angle=45, hjust=1)
          )
p2
ggsave("fig/Exp_Fig3/luecken-logFC-filter-enrichment.pdf", p2,
       width = 18, height = 15, units = "cm", useDingbats = FALSE)

grid = distinct(enr, method, logFC)
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    df = filter(dat, method == current$method, logFC == current$logFC) %>% 
        mutate(expr_set = factor(expr_set, levels = c('logFC', 'Random')))
    eff = cohen.d(formula = number_of_da_regions ~ expr_set | Subject(label),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})
eff$d[is.na(eff$d)] = 0
labs = eff %>% 
    group_by(logFC) %>% 
    summarise(mean = mean(d), n = n(), median = median(d), d = median) %>%
    ungroup() %>%
    mutate(
        median = ifelse(is.na(median), 0, median),
        label = formatC(median, format = 'f', digits = 2)
        )

p3 = eff %>% 
    ggplot(aes(x = factor(logFC), y = d)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.6, size = 0.35,
                 fill = 'grey82', color = 'grey20') + 
    geom_label(data = labs, aes(y = Inf, label = label), size = 1.5, hjust = 0.5,
               vjust = 1, color = 'black', fill = NA, label.size = NA,
               label.padding = unit(0.6, 'lines')) +
    # scale_color_manual('', values = pal) +
    # scale_fill_manual('', values = pal) +
    scale_x_discrete(expression('Absolute fold-change filter')) +
    scale_y_continuous('Cohen\'s d\n(vs. random peak removal)') +
    ggtitle('Effect of log-fold change filtering on false discoveries\nin null comparisons of real scATAC-seq data') +
    coord_cartesian(ylim = c(0, 0.7)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(
        # axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1,
        plot.title = element_text(size=5)
    )
p3
ggsave("fig/Exp_Fig3/luecken-logFC-filter-d.pdf", p3,
       width = 4, height = 4, units = "cm", useDingbats = FALSE)

p4 = wrap_plots(p1,p2,ncol=1)
p4
ggsave("fig/Exp_Fig3/top-row.pdf", p4,
       width = 16, height = 14, units = "cm", useDingbats = FALSE)

