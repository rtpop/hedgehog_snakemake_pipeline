## loading libraries
required_libraries <- c("optparse",
                        "data.table",
                        "pheatmap",
                        "stats",
                        "grid")

for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}

## Options
options(stringsAsFactors = FALSE)

## --------------- ##
## Parse arguments ##
## --------------- ##

option_list <- list(
    optparse::make_option(c("-f", "--file"), type = "character", help = "CSV file with community scores."),
    optparse::make_option(c("-o", "--output-dir"), type = "character", help = "Path to output directory.")
    optparse::make_option(c("-d", "--dist"), type = "character", help = "Path to a file with the distance matrix for clustering."),
    optparse::make_option(c("-b", "--binarise"), type = "logical", default = TRUE, help = "If TRUE, the scores will be binarised for plotting."),
    optparse::make_option(c("-l", "--log-transform"), type = "logical", default = TRUE, help = "If TRUE, the scores will be log transformed for plotting.")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## --------- ##
## Variables ##
## --------- ##
FILE <- opt$`file`
DIST <- opt$`dist`
BIN <- opt$`binarise`
LOG <- opt$`log-transform`
OUT_DIR <- opt$`output-dir`

# source functions
source("src/analysis/cluster_scores_fn.R")

cluster_scores(scores = FILE, distance_matrix = DIST, out_dir = OUT_DIR, binarise = BIN, log_transform = LOG)
