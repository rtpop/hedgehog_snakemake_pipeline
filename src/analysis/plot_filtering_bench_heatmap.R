required_libraries <- c("data.table",
                        "ggplot2",
                        "ggrepel",
                        "ggpubr",
                        "dplyr",
                        "tidyr",
                        "purrr",
                        "RColorBrewer",
                        "stringr",
                        "magrittr",
                        "optparse",
                        "reshape2")
for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}

## options
options(stringsAsFactors = FALSE)

## --------------- ##
## Parse arguments ##
## --------------- ##
option_list <- list(
    optparse::make_option(c("-i", "--input"), type = "character", help = "Consolidated filtering benchmark dataframe file path"),
    optparse::make_option(c("-o", "--output"), type = "character", help = "Path to output directory"),
    optparse::make_option(c("-m", "--metric"), type = "character", default = "Modularity",
                          help = "Filtering method to plot (default: Modularity). Options: 'Modularity', 'density', 'edges'"),
    optparse::make_option(c("-f", "--filtering"), type = "character", default = "all",
                          help = "Filtering method to plot (default: all). Options: 'ELAND filtered PANDA', 'Prior', 'Unfiltered'")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

source('src/analysis/plot_filtering_bench_fn.R')

## --------- ##
## Variables ##
## --------- ##
FILE <- opt$input
OUT <- opt$output
METRIC <- opt$metric
FILTERING <- opt$filtering

p <- plot_heatmap(FILE, metric = METRIC, filtering = FILTERING)

ggsave(filename = OUT, plot = p, width = 5, height = 5)