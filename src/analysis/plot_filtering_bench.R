# load libraries
required_libraries <- c("optparse",
                        "data.table",
                        "ggplot2",
                        "ggrepel",
                        "ggpubr",
                        "dplyr",
                        "tidyr",
                        "purrr",
                        "stringr",
                        "magrittr")
for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}

## options
options(stringsAsFactors = FALSE)

## --------------- ##
## Parse arguments ##
## --------------- ##

option_list <- list(
    optparse::make_option(c("-f", "--file"), type = "character", help = "Path to the file with the filtering benchmark results."),
    optparse::make_option(c("-o", "--output-dir"), type = "character", help = "Path to output directory"),
    optparse::make_option(c("-t", "--tissue-type"), type = "character", help = "Tissue type used in the filtering benchmark."),
    optparse::make_option(c("-d", "--data-frame"), type = "character", help = "File name for compiled benchmark data frame"),,
    optparse::make_option(c("-p", "--plot-type"), type = "character", default = "modularity", help = "Type of plot to generate. Options are 'modularity', 'density', or 'edges'."),
    optparse::make_option(c("-n", "--plot-title"), type = "character", default = "", help = "Title of the plot.")
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
# Check if the required arguments are provided
if (is.null(opt$file) || is.null(opt$output_dir) || is.null(opt$tissue_type)) {
    optparse::print_help(opt_parser)
    stop("Please provide the required arguments: -f, -o, and -t.")
}
# Check if the file exists
if (!file.exists(opt$file)) {
    stop(paste("The file", opt$file, "does not exist."))
}
# Check if the output directory exists, if not create it
if (!dir.exists(opt$output_dir)) {
    dir.create(opt$output_dir, recursive = TRUE)
}
# Check if the plot type is valid
valid_plot_types <- c("modularity", "density", "edges")
if (!opt$plot_type %in% valid_plot_types) {
    stop(paste("Invalid plot type. Choose from:", paste(valid_plot_types, collapse = ", ")))
}

## --------- ##
## Variables ##
## --------- ##
FILE_NAME <- opt$file
OUT_DIR <- opt$output_dir
TISSUE_TYPE <- opt$tissue_type
DATA_FRAME <- opt$data_frame
PLOT_TYPE <- opt$plot_type
PLOT_TITLE <- opt$plot_title

## source functions
source("src/analysis/plot_filtering_bench_fn.R")

## Run
if (!file.exists(DATA_FRAME)) {
    df <- prepare_filtering_bench(file_name = FILE_NAME, tissue_type = TISSUE_TYPE)
    data.table::fwrite(df, file = DATA_FRAME, sep = "\t", row.names = FALSE)
} else {
    df <- data.table::fread(DATA_FRAME, header = TRUE)

}

plot_filtering_bench(df = df, output_dir = OUT_DIR, plot_type = PLOT_TYPE, plot_title = PLOT_TITLE)

