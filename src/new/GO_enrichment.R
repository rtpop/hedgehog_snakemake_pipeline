# Snakemake entry point: human biological-process enrichment of selected GMT sets.
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(data.table))
source(snakemake@params$function_file)

run_go_enrichment(
    gmt_file = snakemake@input$gmt,
    summary_file = snakemake@output$go_enrichment,
    details_dir = snakemake@output$details,
    auto_bg = snakemake@params$auto_bg,
    bg_file = if (length(snakemake@input$bg)) snakemake@input$bg[[1]] else NULL,
    save_all = snakemake@params$save_all,
    sig_thresh = snakemake@params$sig_thresh,
    statistic = snakemake@params$statistic,
    algorithm = snakemake@params$algorithm
)
