# load libraries
required_libraries <- c("optparse",
                        "data.table",
                        "ggplot2",
                        "ggrepel",
                        "ggpubr",
                        "dplyr",
                        "tidyr",
                        "purrr",
                        "RColorBrewer",
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
    optparse::make_option(c("-f", "--files"), type = "character", help = "Comma-separated list of paths to the filtering benchmark result files."),
    optparse::make_option(c("-o", "--output-dir"), type = "character", help = "Path to output directory"),
    optparse::make_option(c("-t", "--tissue-type"), type = "character", help = "Tissue type used in the filtering benchmark."),
    optparse::make_option(c("-d", "--data-frame"), type = "character", help = "File name for compiled benchmark data frame"),
    optparse::make_option(c("-p", "--plot-type"), type = "character", default = "modularity", help = "Type of plot to generate. Options are 'modularity', 'density', or 'edges'."),
    optparse::make_option(c("-n", "--plot-title"), type = "character", default = "", help = "Title of the plot."),
    optparse::make_option(c("-l", "--plot-file"), type = "character", default = "", help = "Path to save the plot file."),
    optparse::make_option(c("-i", "--include-unfiltered"), type = "logical", default = TRUE, help = "Include the 'Unfiltered' network in the plot.")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Check if the required arguments are provided
if (is.null(opt$files) || is.null(opt$`output-dir`) || is.null(opt$`tissue-type`)) {
    optparse::print_help(opt_parser)
    stop("Please provide the required arguments: -f, -o, and -t.")
}

# Parse the input files
input_files <- strsplit(opt$files, ",")[[1]]

# Check if the output directory exists, if not create it
if (!dir.exists(opt$`output-dir`)) {
    dir.create(opt$`output-dir`, recursive = TRUE)
}

# Check if the plot type is valid
valid_plot_types <- c("modularity", "density", "edges", "all")
if (!opt$`plot-type` %in% valid_plot_types) {
    stop(paste("Invalid plot type. Choose from:", paste(valid_plot_types, collapse = ", ")))
}

## --------- ##
## Variables ##
## --------- ##
OUT_DIR <- opt$`output-dir`
TISSUE_TYPE <- opt$`tissue-type`
DATA_FRAME <- opt$`data-frame`
PLOT_TYPE <- opt$`plot-type`
PLOT_TITLE <- opt$`plot-title`
PLOT_FILE <- opt$`plot-file`
UNFILTERED <- opt$`include-unfiltered`

if (PLOT_TITLE == "") {
    PLOT_TITLE <- paste("Filtering Benchmark Results for", TISSUE_TYPE)
}

## source functions
source(file.path("src", "analysis", "plot_filtering_bench_fn.R"))

## Plot the consolidated data
plot_filtering_bench(df = consolidated_df, output_dir = OUT_DIR, plot_type = PLOT_TYPE, plot_title = PLOT_TITLE, plot_file = PLOT_FILE, include_unfiltered = UNFILTERED)

