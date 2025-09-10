required_libraries <- c("optparse",
                        "data.table")

for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}

## Options
options(stringsAsFactors = FALSE)

## --------------- ##
## Parse arguments ##
## --------------- ##

option_list <- list(
    optparse::make_option(c("-i", "--input"), type = "character", help = "RData file containing GTEx data."),
    optparse::make_option(c("-e", "--edgelist"), type = "logical", default = TRUE, help = "If edges should be extracted for each tissue. If FALSE, expression will be extracted instead."),
    optparse::make_option(c("-d", "--out-dir"), type = "character", default = "data/", help = "Path to output directory.")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## --------- ##
## Variables ##
## --------- ##
DATA <- opt$`input`
EXTRACT_EDGES <- opt$`extract-edges`
OUT_DIR <- opt$`out-dir`

# source functions
source("src/process_data/process_gtex_fn.R")

get_gtex_data(data = DATA, extract_edges = EXTRACT_EDGES, out_dir = OUT_DIR)