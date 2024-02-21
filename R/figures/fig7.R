setwd("git/DA-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(effsize)
source("R/theme.R")

###############################################################################-
## Bulk downsampled ####
###############################################################################-

# load data
dat = readRDS("data/summaries/concordance_res/matched_bulk/all_aucc_downsampled.rds") %>% 
  type_convert() %>% 
  filter(!sc_pseudo_repl)
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
         color = sc_da_family,
         color = ifelse(method %in% c('binomial', 'FET', 'LRpeaks', 'permutation'),
                        'singlecell', color) %>% 
           fct_recode('single-cell' = 'singlecell',
                      'other' = 'non_libra') %>% 
           fct_relevel('single-cell', 'pseudobulk', 'other'))

# re-code filtering
avg %<>%
  mutate(filter = ifelse(sc_percent_filter, 
                         paste0(sc_min_cell_filter, '%'),
                         paste0(sc_min_cell_filter, ' cells')) %>% 
           fct_recode('1 cell' = '1 cells') %>% 
           fct_relevel('1 cell', '10 cells', '0.2%', '0.5%', '1%', '5%', '10%', '20%')) %>% 
  filter(!filter %in% c('10%', '20%')) %>% 
  mutate(method = fct_recode(method,
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
  filter(!grepl('limma', method))

avg0 = filter(avg, k == 5000, sc_min_cell_filter == 1, !sc_percent_filter)

grid = distinct(avg0, method, downsample_proportion) %>% 
  filter(downsample_proportion < 10000)
eff = pmap_dfr(grid, function(...) {
  current = tibble(...)
  df = filter(avg0, method == current$method, 
              downsample_proportion %in% c(10000, current$downsample_proportion))
  eff = cohen.d(formula = aucc ~ factor(downsample_proportion) | Subject(cell_type),
                data = df,
                paired = TRUE)$estimate
  mutate(current, d = eff)
})

labs = eff %>% 
  group_by(downsample_proportion) %>% 
  summarise(mean = mean(d), n = n(), median = median(d), d = median) %>%
  ungroup() %>%
  mutate(label = formatC(median, format = 'f', digits = 2))
p1 = eff %>% 
  ggplot(aes(x = factor(downsample_proportion), y = d)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.6, size = 0.35,
               fill = 'grey82', color = 'grey20') + 
  geom_text(data = labs, aes(y = 0.2, label = label), size = 1.75, hjust = 0.5,
            color = 'black') +
  # scale_color_manual('', values = pal) +
  # scale_fill_manual('', values = pal) +
  scale_x_discrete('# of fragments') +
  scale_y_continuous('Cohen\'s d (vs. 10,000 fragments)') +
  ggtitle('Effect of downsampling fragments on concordance\nbetween scATAC-seq and matched bulk ATAC-seq') +
  boxed_theme(size_sm = 5, size_lg = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1,
        plot.title = element_text(size=5, hjust=0.5)
        )
p1
ggsave("fig/Fig7/downsample-reads-d.pdf", p1,
       width = 6, height = 6, units = "cm", useDingbats = FALSE)

###############################################################################-
## Multiome downsampled ####
###############################################################################-

# load data
dat = readRDS("data/summaries/concordance_res/multiome/all_aucc_downsampled.rds") %>% 
    type_convert() 
# remove mixed models and pseudoreplicates
dat %<>% filter(!scatac_pseudo_repl, scatac_da_family != 'mixedmodel')
# average over bulk methods
avg = dat %>% 
    group_by_at(vars(-scrna_da_method, -bulk_da_type, -sc_features,
                     -bulk_features, -overlap, -aucc)) %>% 
    summarise(aucc = mean(aucc),
              n = n()) %>% 
    ungroup()

# re-code x-values and colors
avg %<>%
    mutate(method = paste0(scatac_da_method, '-', scatac_da_type) %>% 
               gsub("-singlecell|-binomial|-exact", "", .),
           method = ifelse(scatac_pseudo_repl,
                           paste(method, '(pseudo-replicates)'), method),
           method = gsub("snapatac", "SnapATAC::findDAR", method) %>% 
               gsub("fisher", "FET", .) %>% 
               gsub("-LR_.*|-perm.*$", "", .) %>% 
               gsub("LR_", "LR", .),
           color = ifelse(scatac_pseudo_repl, 'pseudo-replicates', scatac_da_family),
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
    filter(!grepl('limma', method))
# re-code filtering
avg %<>%
    mutate(filter = ifelse(scatac_percent_filter, 
                           paste0(scatac_min_cell_filter, '%'),
                           paste0(scatac_min_cell_filter, ' cells')) %>% 
               fct_recode('1 cell' = '1 cells') %>% 
               fct_relevel('1 cell', '10 cells', '0.2%', '0.5%', '1%', '5%', 
                           '10%', '20%')) %>% 
    filter(!filter %in% c('10%', '20%'))

avg0 = filter(avg, k == 500, scatac_min_cell_filter == 1, !scatac_percent_filter)

grid = distinct(avg0, method, downsample_proportion) %>% 
    filter(downsample_proportion < 1000)
eff = pmap_dfr(grid, function(...) {
    current = tibble(...)
    df = filter(avg0, method == current$method, 
                downsample_proportion %in% c(1000, current$downsample_proportion))
    eff = cohen.d(formula = aucc ~ factor(downsample_proportion) | Subject(cell_type),
                  data = df,
                  paired = TRUE)$estimate
    mutate(current, d = eff)
})
labs = eff %>% 
    group_by(downsample_proportion) %>% 
    summarise(mean = mean(d), n = n(), median = median(d), d = median) %>%
    ungroup() %>%
    mutate(label = formatC(median, format = 'f', digits = 2))
p2 = eff %>% 
    ggplot(aes(x = factor(downsample_proportion), y = d)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.6, size = 0.35,
                 fill = 'grey82', color = 'grey20') + 
    geom_text(data = labs, aes(y = 0.2, label = label), size = 1.75, hjust = 0.5,
              color = 'black') +
    ggtitle('Effect of downsampling cells on concordance\nbetween ATAC and RNA modalities of multiome data') +
    scale_x_discrete('# of cells') +
    scale_y_continuous('Cohen\'s d (vs. 1,000 cells)') +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = 1,
          plot.title = element_text(size=5)
          )
p2
ggsave("fig/Fig7/downsample-cells-d.pdf", p2,
       width = 6, height = 6, units = "cm", useDingbats = FALSE)

p3 = p1 + theme(aspect.ratio = 1.1, plot.margin = unit(c(0.01, 1, 1, 0.01), "cm")) | p2 + theme(aspect.ratio = 1.1, plot.margin = unit(c(0.01, 1, 1, 0.01), "cm"))
p3
ggsave("fig/Fig7/downsample-combined.pdf", p3,
       width = 12, height = 6, units = "cm", useDingbats = FALSE)
