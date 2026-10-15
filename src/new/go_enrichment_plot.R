suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
source(snakemake@params$function_file)
plot_n_sig_terms(snakemake@input$go_enrichment, snakemake@output$go_enrichment_plot)
