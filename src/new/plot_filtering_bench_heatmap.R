required_libraries <- c("data.table", "ggplot2", "scales")
for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}

## options
options(stringsAsFactors = FALSE)

FILE <- snakemake@input$filtering_bench_df
OUT <- snakemake@output$benchmark_heatmap
METRIC <- snakemake@params$metric
FILTERING <- snakemake@params$filtering_method
SEPARATE_FILES <- snakemake@params$separate_files
FUNCTION_FILE <- snakemake@params$function_file

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(FUNCTION_FILE)

if (FILTERING == "all" && SEPARATE_FILES) {
  # Read the data to get available filtering methods
  df <- data.table::fread(FILE, header = TRUE)
  filtering_methods <- unique(df$Network)
  lapply(filtering_methods, function(filt) {
    p <- plot_heatmap(FILE, metric = METRIC, filtering = filt)
    out_file <- file.path(OUT, paste0("heatmap_", filt, ".pdf"))
    str(OUT)
    ggsave(filename = out_file, plot = p, width = 5, height = 5)
    
  })
} else {
  p <- plot_heatmap(FILE, metric = METRIC, filtering = FILTERING)
  out_file <- file.path(OUT, paste0("heatmap_", FILTERING, ".pdf"))
  ggsave(filename = out_file, plot = p, width = 5, height = 5)
}
