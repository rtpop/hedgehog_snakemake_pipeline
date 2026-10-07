required_libraries <- c("data.table")

for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}

## options
options(stringsAsFactors = FALSE)

input_files <- as.character(snakemake@input$filtering_bench_dfs)
output_file <- snakemake@output$filtering_benchmark_consolidated
function_file <- snakemake@params$function_file

## source functions
source(function_file)

## Consolidate data
consolidate_data(input_files, output_file)

# End of the script