suppressPackageStartupMessages(library(data.table))

input_files <- as.character(snakemake@input$filtering_bench_list)
output_file <- snakemake@output$filtering_bench_df
tissue_type <- snakemake@params$tissue_type
function_file <- snakemake@params$function_file

source(function_file)

# Apply prepare_filtering_bench to each file and combine
df_list <- lapply(input_files, function(file) {
    prepare_filtering_bench(file_name = file, tissue_type = tissue_type)
})
consolidated_df <- data.table::rbindlist(df_list)

data.table::fwrite(consolidated_df, file = output_file, sep = "\t", row.names = FALSE)