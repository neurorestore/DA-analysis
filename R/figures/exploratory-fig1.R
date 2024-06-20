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
dat = readRDS("data/summaries/concordance_res/matched_bulk/all_aucc_logFC_absolute_filter.rds") %>% 
    type_convert() %>% 
    # remove pseudoreplicates
    filter(!sc_pseudo_repl)
# remove fixed columns
keep = map_lgl(dat, ~ n_distinct(.x) > 1)
dat %<>% extract(, keep)

# restructure 'expr' column
dat %<>% mutate(expr_set = ifelse(grepl('pseudo', expr), 'Random', 'logFC') %>% 
                    fct_relevel('Random', 'logFC'),
                expr = gsub("_pseudo", "", expr))

# average over bulk methods
avg = dat %>% 
    group_by_at(vars(-bulk_da_method, -bulk_da_type,
                     -bulk_features, -overlap, -aucc)) %>% 
    summarise(aucc = mean(aucc),
              n = n()) %>% 
    ungroup()

    # re-code x-values and colors
avg %<>%
    mutate(method = paste0(sc_da_method, '-', sc_da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .),
           color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                          'singlecell', sc_da_family) %>% 
               fct_recode('single-cell' = 'singlecell',
                          'other' = 'non_libra') %>% 
               fct_relevel('single-cell', 'pseudobulk', 'other'))

# re-code filtering
avg %<>%
    mutate(method = fct_recode(method,
                               "'SnapATAC::findDAR'" = 'SnapATAC::findDAR',
                               'Permutation~test' = 'permutation',
                               't~test' = 't',
                               'LR[peaks]' = 'LR_peaks',
                               'LR[clusters]' = 'LR',
                               'Binomial' = 'binomial',
                               'Wilcoxon~rank-sum~test' = 'wilcox',
                               'Fisher~exact~test' = 'FET',
                               'Negative~binomial' = 'negbinom') %>% 
               as.character()) %>% 
    # drop limma
    filter(!grepl('limma', method))

avg1 = avg %>% 
    filter(!grepl('pseudo-replicates', method)) %>% 
    filter(!(paper == 'GonzalesBlas_2018' & !grepl('0 h', cell_type))) %>% 
    mutate(logFC = gsub("^.*_", "", expr) %>% 
               replace(. == 'all', 0))

dplyr::count(avg1, method)

labs = avg1 %>% 
    group_by(logFC, expr_set, method, color) %>%
    summarise(mean = mean(aucc), n = n(), median = median(aucc), aucc = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
pal = c('logFC' = 'black',
        'Random' = 'grey75') 

p1 = avg1 %>%
    filter(
        logFC != 1, logFC != 0
    ) %>%
    filter(
        !(method == 'Negative~binomial' & aucc > 0.1),
        !(method != 'Negative~binomial' & aucc > 0.01)
    ) %>%
    ggplot(aes(x = logFC, y = aucc,
           color = expr_set, fill = expr_set)) +
    facet_wrap(~ method, scales = 'free', ncol = 7, labeller = label_parsed) +
    geom_boxplot(aes(y = aucc), outlier.shape = NA, alpha = 0.4, width = 0.7,
                 size = 0.3) + 
    # geom_errorbar(data = labs, aes(ymin = mean, ymax = mean), width = 0.8) +
    # geom_text(data = labs, aes(y = -0.12, label = label), size = 1.75, hjust = 0,
    #           color = 'black') +
    scale_color_manual('Filter by:', values = pal) +
    scale_fill_manual('Filter by:', values = pal) +
    scale_x_reordered(expression('Absolute fold-change filter')) +
    scale_y_continuous('AUCC') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(legend.position = 'top',
          legend.key.width = unit(0.4, 'lines'),
          # aspect.ratio = 1.2,
          # axis.text.x = element_text(angle=45, hjust=1)
      )
p1
ggsave("fig/Exp_Fig1/bulk-concordance-absolute-logFC-filter-overview.pdf", p1,
       width = 16, height = 6, units = "cm", useDingbats = FALSE)

###############################################################################-
## AUCC: boxplot, enrichment over random ####
###############################################################################-

enr = avg1 %>%
    filter(
        logFC != 1, logFC != 0
    ) %>% 
    group_by_at(vars(-expr_set, -aucc)) %>% 
    summarise(delta = ifelse(n() == 2, 
                             aucc[expr_set == 'logFC'] - aucc[expr_set == 'Random'],
                             aucc),
              nn = n()) %>% 
    ungroup() 
pal = da_analysis_colors
p2 = enr %>%
    filter(logFC != '0.0') %>% 
    remove_quantile_per_method(., val_col='delta', group_cols = c('method', 'logFC')) %>%
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
    scale_y_continuous(expression(Delta~AUCC)) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(legend.position = 'none',
          # aspect.ratio = 0.8,
          # axis.text.x = element_text(angle=45, hjust=1),
          axis.title.y = element_text(size=6)
    )
p2
ggsave("fig/Exp_Fig1/bulk-concordance-absolute-logFC-filter-enrichment.pdf", p2,
       width = 16, height = 6, units = "cm", useDingbats = FALSE)

grid = distinct(enr, method, logFC) %>% filter(logFC != '0.0')
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    df = filter(avg1, method == current$method, logFC == current$logFC) %>% 
        mutate(expr_set = factor(expr_set, levels = c('logFC', 'Random')))
    eff = cohen.d(formula = aucc ~ expr_set | Subject(paste(paper, cell_type)),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})
labs = eff %>% 
    filter(!is.na(d)) %>%
    group_by(logFC) %>% 
    summarise(mean = mean(d), n = n(), median = median(d), d = median) %>%
    ungroup() %>%
    mutate(label = format(median, digits = 2, scientific=T))

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
    scale_y_continuous('Cohen\'s d\n(vs. random peak removal)', limits = c(-1e-04, 1e-04)) +
    ggtitle('Effect of fold-change filtering on concordance\nbetween scATAC-seq and matched bulk ATAC-seq') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(
        # axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1,
        plot.title = element_text(size=5, hjust=0.5)
        ) 
p3
ggsave("fig/Exp_Fig1/bulk-concordance-absolute-logFC-filter-d.pdf", p3,
       width = 4, height = 4, units = "cm", useDingbats = FALSE)

p4 = wrap_plots(p1,p2,ncol=1)
ggsave("fig/Exp_Fig1/top-row.pdf", p4,
       width = 16, height = 12, units = "cm", useDingbats = FALSE)

