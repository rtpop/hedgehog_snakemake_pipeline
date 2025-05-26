## loading libraries
required_libraries <- c("optparse",
                        "data.table",
                        "ggplot2")

for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}

## Options
options(stringsAsFactors = FALSE)

## --------------- ##
## Parse arguments ##
## --------------- ##
option_list <- list(
    optparse::make_option(c("-i", "--input"), type = "character", help = "Path to input file"),
    optparse::make_option(c("-o", "--output"), type = "character", help = "Path to output file")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


## --------- ##
## Variables ##
## --------- ##
INPUT <- opt$input
OUT_FILE <- opt$output

## source functions
source("src/analysis/community_plot_fn.R")

## Run
p <- plot_n_communities(INPUT, file_name = OUT_FILE)
