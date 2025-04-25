required_libraries <- c("ggplot2", "data.table", "optparse")

# Load required libraries
for (lib in required_libraries) {
  suppressMessages(library(lib, character.only = TRUE))
}

## Options
options(stringsAsFactors = FALSE)

# --------------- #
# Parse arguments #
# --------------- #

# Define command line options
option_list <- list(
  optparse::make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to the input file with the results of the GO enrichment analysis"),
  optparse::make_option(c("-o", "--out-dir"), type = "character", default = NULL,
              help = "Path to output directory"),
  optparse::make_option(c("-f", "--file-name"), type = "character", default = NULL,
              help = "Name of the output file")
)
# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## --------- ##
## Variables ##
## --------- ##

INPUT_FILE <- opt$input
OUT_DIR <- opt$`out-dir`
FILE_NAME <- opt$`file-name`

# source functions
source("src/analysis/go_enrichment_plot_fn.r")

print(getwd())

plot_n_sig_terms(input_file = INPUT_FILE, out_dir = OUT_DIR, file_name = FILE_NAME)
