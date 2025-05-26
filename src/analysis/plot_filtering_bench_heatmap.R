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
    optparse::make_option(c("-f", "--file"), type = "character", help = "Consolidated filtering benchmark dataframe file path"),
    optparse::make_option(c("-o", "--output"), type = "character", help = "Path to output directory")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

source('src/analysis/plot_filtering_bench_fn.R')

## --------- ##
## Variables ##
## --------- ##
FILE <- opt$file
OUT <- opt$output

p <- plot_heatmap(FILE)

ggsave(filename = OUT, plot = p, width = 10, height = 8)