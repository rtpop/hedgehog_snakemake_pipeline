required_libraries <- c("data.table", "optparse")

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
    optparse::make_option(c("-o", "--output"), type = "character", help = "Path to output file")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# parse the input files
input_files <- strsplit(opt$files, ",")[[1]]

## --------- ##
## Variables ##
## --------- ##
OUT <- opt$output

## source functions
source(file.path("src", "utils", "consolidate_benchmark_fn.R"))

## Consolidate data
df_list <- lapply(input_files, function(file) {
    if (!file.exists(file)) {
        stop(paste("The file", file, "does not exist."))
    }
    prepare_filtering_bench(file_name = file, tissue_type = TISSUE_TYPE)
    })
    consolidated_df <- data.table::rbindlist(df_list)

# Save the consolidated data frame
data.table::fwrite(consolidated_df, file = DATA_FRAME, sep = "\t", row.names = FALSE)

consolidated_data <- consolidate_data(df_list, OUT)

# End of the script