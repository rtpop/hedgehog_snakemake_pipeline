suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
    optparse::make_option(c("-f", "--files"), type = "character", help = "Comma-separated list of input files"),
    optparse::make_option(c("-o", "--output"), type = "character", help = "Output file"),
    optparse::make_option(c("-t", "--tissue"), type = "character", help = "Tissue type")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

input_files <- strsplit(opt$files, ",")[[1]]
output_file <- opt$output
tissue_type <- opt$tissue

source(file.path("src", "utils", "consolidate_benchmark_fn.R"))

# Apply prepare_filtering_bench to each file and combine
df_list <- lapply(input_files, function(file) {
    prepare_filtering_bench(file_name = file, tissue_type = tissue_type)
})
consolidated_df <- data.table::rbindlist(df_list)

data.table::fwrite(consolidated_df, file = output_file, sep = "\t", row.names = FALSE)