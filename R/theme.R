library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(scales)
library(patchwork)
library(paletteer)
library(ggsci)
library(scico)
library(colorspace)
# library(drlib)
library(cowplot)
# library(nationalparkcolors)
# library(fishualize)
# library(BuenColors)
# library(ggstance)
# library(shadowtext)
library(ggrastr)

# Define theme
clean_theme = function(size_lg = 6, size_sm = 5) {
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = size_sm),
        axis.text.y = element_text(size = size_sm),
        axis.ticks.length.x = unit(0.15, 'lines'),
        axis.ticks.length.y = unit(0.15, 'lines'),
        axis.title.x = element_text(size = size_lg),
        axis.title.y = element_text(size = size_lg),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size = size_sm),
        # strip.background = element_rect(fill = "grey90", color = "grey90",
        #                                 size = 0),
        strip.background = element_blank(),
        axis.line.y = element_line(colour = "grey50"),
        axis.line.x = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "grey50"),
        legend.position = "top",
        legend.text = element_text(size = size_sm),
        legend.title = element_text(size = size_sm),
        legend.key.size = unit(0.6, "lines"),
        legend.margin = margin(rep(0, 4)),
        # legend.box.margin = ggplot2::margin(rep(0, 4), unit = 'lines'),
        # legend.box.spacing = ggplot2::margin(rep(0, 4)),
        legend.background = element_blank(),
        plot.title = element_text(size = size_lg, hjust = 0.5),
        axis.ticks.length = unit(2, 'pt'),)
}

grid_theme = function(size_lg = 6, size_sm = 5) {
  theme_bw() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = size_sm),
          axis.text.y = element_text(size = size_sm),
          axis.title.x = element_text(size = size_lg),
          axis.title.y = element_text(size = size_lg),
          strip.text = element_text(size = size_sm),
          strip.background = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks = element_line(colour = "grey50"),
          legend.position = "top",
          legend.text = element_text(size = size_sm),
          legend.title = element_text(size = size_sm),
          legend.key.size = unit(0.6, "lines"),
          legend.margin = margin(rep(0, 4)),
          # legend.box.margin = ggplot2::margin(rep(0, 4), unit = 'lines'),
          # legend.box.spacing = ggplot2::margin(rep(0, 4)),
          legend.background = element_blank(),
          plot.title = element_text(size = size_lg, hjust = 0.5),
          axis.ticks.length = unit(2, 'pt'),)
}

boxed_theme = function(size_lg = 6, size_sm = 5) {
  theme_bw() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = size_sm),
          axis.text.y = element_text(size = size_sm),
          axis.title.x = element_text(size = size_lg),
          axis.title.y = element_text(size = size_lg),
          panel.grid = element_blank(),
          strip.text = element_text(size = size_sm),
          strip.background = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks = element_line(colour = "grey50"),
          legend.position = "top",
          legend.text = element_text(size = size_sm),
          legend.title = element_text(size = size_sm),
          legend.key.size = unit(0.6, "lines"),
          legend.margin = margin(rep(0, 4)),
          # legend.box.margin = ggplot2::margin(rep(0, 4), unit = 'lines'),
          # legend.box.spacing = ggplot2::margin(rep(0, 4)),
          legend.background = element_blank(),
          plot.title = element_text(size = size_lg, hjust = 0.5),
          axis.ticks.length = unit(2, 'pt'),)
}

umap_theme = function(size_lg = 6, size_sm = 5) {
  boxed_theme(size_lg = size_lg, size_sm = size_sm) +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_text(hjust = 0),
          axis.title.x = element_text(hjust = 0))
}

darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}


lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  for (value in col){
    if (value > 255) {
      col[col == value] = 255
    }
  }
  col <- rgb(t(col), maxColorValue=255)
  col
}

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # remove zero
  l <- gsub("0e\\+00", "0", l)
  # remove one
  l <- gsub("^1e\\+00", "1", l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + from exponent
  l <- gsub("e\\+" ,"e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # remove 1 x 10^ (replace with 10^)
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

plot_pal = function(pal) {
  grid::grid.raster(pal, interpolate=F)
}

cubehelix = function(n_colors) {
  colours = c("#000000", "#1A1935", "#15474E", "#2B6F39", "#767B33", "#C17A6F",
              "#D490C6", "#C3C0F2")
  idxs = 0.3
  if (n_colors > 1)
    idxs = seq(0, 1, 1 / (n_colors - 1))
  colour_ramp(colours)(idxs)
}

kinney6 = c(
  "#c6c3bf",
  "#119e87",
  "#53bad3",
  "#559ed2",
  "#3b5687",
  "#e34d3b")

Gpal = c("#E30F17", "#0296E1", "#F49203", "#E5E5E5", "#E0A4D1") %>%
  # extended version
  ## https://www.sciencedirect.com/science/article/pii/S0896627316000106
  c('#2C8942', '#E8B820',
    # https://www.nature.com/articles/nature20118.pdf
    '#9e2062', '#a13027', '#57c3f0',
    # "Mechanisms Underlying the Neuromodulation...", Fig 1
    '#3b54a5', '#2b68b2', '#1c99b3', '#39ba92', '#6fc06b'
    )

colours.cafe447 = c('#ffb838', '#fee5a5', '#f7f6fee', '#486d87')
colours.cafe433 = c('#077893', '#e3deca', '#fcfaf1', '#ff9465')
colours.cafe425 = c('#2B5B6C', '#C7CFAC', '#FCFAF1', '#E34F33', '#FFC87E')
colours.cafe322 = c("#7bbaea", "#d46363", "#fbdaa7", "#fcfaf2", "#30598c")

winsorize = function(vec, limits) {
  lower_limit = limits[1]
  upper_limit = limits[2]
  if (!is.na(upper_limit))
    vec[vec > upper_limit] = upper_limit
  if (!is.na(lower_limit))
    vec[vec < lower_limit] = lower_limit
  return(vec)
}

nr_base_3 = c("#19D3C5", "#F0DF0D", "#FA525B")
nr_base_4 = c("#FA525B", "#372367", "#19D3C5", "#F0DF0D")
nr_base_5 = c("#FA525B", "#969899", "#372367", "#19D3C5", "#F0DF0D")
nr_base_6 = c("#088980", "#6C005F", "#FA525B", "#372367", "#19D3C5", "#F0DF0D")
nr_base_11 = rev(c("#FA525B", "#372367", "#19D3C5", "#F0DF0D", "#B8E600", "#FFB11E", "#218716", "#FE9E53","#A40A11", "#6C005F", "#088980"))
nr_pal_greys = c('#DBDBDD', '#BBBABA', '#969899', '#757679')
nr_heat_red = colorRampPalette(c("#000000", "#DBDBDD", "#FFFFFF", "#FA525B", "#A31B21"), interp = "spline")(100)
nr_heat_blue = colorRampPalette(c("#000000", "#DBDBDD", "#FFFFFF", "#19D3C5", "#088980"), interp = "spline")(100)
nr_heat_yellow = colorRampPalette(c("#000000", "#DBDBDD", "#FFFFFF", "#F0DF0D", "#F0DF0D"), interp = "linear")(100)
nr_base_11_light = rev(c('#F69DA4', '#887AA1', '#ADDEDA', '#F8F3AD', '#E0EBAE', '#FFDCA4', '#A7C2A2', '#FED5B1', '#D39986', '#B092AA', '#A5C4C1'))
nr_heat_red_no_white = colorRampPalette(c("#000000", "#DBDBDb", "#DBDBDb", "#DBDBDb", "#FA525B", "#A31B21"), interp = "linear")(100)
nr_heat_blue_no_white = colorRampPalette(c("#000000", "#DBDBDD", "#19D3C5", "#088980"), interp = "linear")(100)
nr_heat_blue_spatial = colorRampPalette(c("#FFFFFF", "#DBDBDD", "#19D3C5", "#088980"), interp = "linear")(100)
nr_heat_yellow_spatial = colorRampPalette(c("#DBDBDD","#FFFFFF", "#FFFFFF", "#F0DF0D", "#FFB11E"), interp = "linear")(100)
nr_heat_red_spatial = colorRampPalette(c("#DBDBDD","#DBDBDD", "#FFFFFF", "#FA525B", "#FD1E29"), interp = "spline")(100)
nr_tree_red = colorRampPalette(c("#757679", "#DBDBDD", "#FFFFFF", "#FA525B", "#A31B21"), interp = "linear")(100)
# semi-shuffle a NR palette
shuffle_nr = function(pal, nbin) {
  set.seed(2)
  tmp = data.frame(colors = pal,
                   bin = rep(seq(1, (length(pal)/nbin)),
                             each = (length(pal)/nbin))[1:length(pal)]
  ) %>%
    group_by(bin) %>%
    mutate(colors = sample(colors)) %>%
    pull(colors)
  return(tmp)
}
# create a bargraph function
nr_bargraph = function(data, metrics, cols, labels, scale_y = c(0, 10)){
  ggplot(data %>% filter(var == metrics), aes(x = Group, y = val, color = Group, fill = Group)) +
    stat_summary(geom = "bar", fun.y = "mean",
                 size = 0.2, width = 0.8, color = NA) +
    geom_jitter(width = .1, size = 1.5) +
    guides(fill = F, color = F) +
    scale_color_manual(values = darken(cols)) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(name = labels[2], expand = c(0,0), limits = scale_y) +
    scale_x_discrete(name = labels[1]) +
    clean_theme() +
    theme(axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5))
}

nr_base_3 = c("#19D3C5", "#F0DF0D", "#FA525B")
nr_base_4 = c("#FA525B", "#372367", "#19D3C5", "#F0DF0D")
nr_base_5 = c("#FA525B", "#969899", "#372367", "#19D3C5", "#F0DF0D")
nr_base_6 = c("#088980", "#6C005F", "#FA525B", "#372367", "#19D3C5", "#F0DF0D")
nr_base_11 = rev(c("#FA525B", "#372367", "#19D3C5", "#F0DF0D", "#B8E600", "#FFB11E", "#218716", "#FE9E53","#A40A11", "#6C005F", "#088980"))
nr_pal_greys = c('#DBDBDD', '#BBBABA', '#969899', '#757679')
nr_heat_red = colorRampPalette(c("#000000", "#DBDBDD", "#FFFFFF", "#FA525B", "#A31B21"), interp = "spline")(100)
nr_heat_blue = colorRampPalette(c("#000000", "#DBDBDD", "#FFFFFF", "#19D3C5", "#088980"), interp = "spline")(100)
nr_heat_yellow = colorRampPalette(c("#000000", "#DBDBDD", "#FFFFFF", "#F0DF0D", "#F0DF0D"), interp = "linear")(100)
nr_base_11_light = rev(c('#F69DA4', '#887AA1', '#ADDEDA', '#F8F3AD', '#E0EBAE', '#FFDCA4', '#A7C2A2', '#FED5B1', '#D39986', '#B092AA', '#A5C4C1'))
nr_heat_red_no_white = colorRampPalette(c("#000000", "#DBDBDb", "#DBDBDb", "#DBDBDb", "#FA525B", "#A31B21"), interp = "linear")(100)
nr_heat_blue_no_white = colorRampPalette(c("#000000", "#DBDBDD", "#19D3C5", "#088980"), interp = "linear")(100)
nr_heat_blue_spatial = colorRampPalette(c("#FFFFFF", "#DBDBDD", "#19D3C5", "#088980"), interp = "linear")(100)
nr_heat_yellow_spatial = colorRampPalette(c("#DBDBDD","#FFFFFF", "#FFFFFF", "#F0DF0D", "#FFB11E"), interp = "linear")(100)
nr_heat_red_spatial = colorRampPalette(c("#DBDBDD","#DBDBDD", "#FFFFFF", "#FA525B", "#FD1E29"), interp = "spline")(100)
nr_tree_red = colorRampPalette(c("#757679", "#DBDBDD", "#FFFFFF", "#FA525B", "#A31B21"), interp = "linear")(100)

x100 = function(x) x * 100

# da_analysis_colors = c(
#     # singlecell colors 
#     colorRampPalette(c(pals::stepped3()[1], pals::stepped3()[4], pals::stepped3()[13]))(8) %>% setNames(c("t~test", "Binomial", "Fisher~exact~test", "LR[peaks]", "LR[clusters]", "Permutation~test", "Negative~binomial", "Wilcoxon~rank-sum~test")),
#     # pseudobulk colors
#     colorRampPalette(pals::stepped3()[5:7])(4) %>% setNames(c("DESeq2-LRT", "DESeq2-Wald", "edgeR-LRT", "edgeR-QLF")),
#     # other colors
#     pals::stepped3()[11] %>% setNames("'SnapATAC::findDAR'")
#     )

da_analysis_colors = c(
    # singlecell colors 
    # colorRampPalette(c('#AAF2EC', '#02CBBC', '#088980', pals::stepped3()[1]))(8) %>% setNames(c("t~test", "Binomial", "Fisher~exact~test", "LR[peaks]", "LR[clusters]", "Permutation~test", "Negative~binomial", "Wilcoxon~rank-sum~test")),
    colorRampPalette(c(pals::stepped3()[1], pals::stepped3()[4], '#00BFAF'))(8) %>% setNames(c("t~test", "Binomial", "Fisher~exact~test", "LR[peaks]", "LR[clusters]", "Permutation~test", "Negative~binomial", "Wilcoxon~rank-sum~test")),
    # pseudobulk colors
    colorRampPalette(c('#FC999E', '#F04B54', '#A31B21'))(4) %>% setNames(c("DESeq2-LRT", "DESeq2-Wald", "edgeR-LRT", "edgeR-QLF")),
    # other colors
    '#218716' %>% setNames("'SnapATAC::findDAR'")
)
