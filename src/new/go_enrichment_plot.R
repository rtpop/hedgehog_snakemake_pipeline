suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))

# Consolidation counts each selected community once when n_sig_terms > 0.
# The tissue summary includes tissues with zero enriched communities.
results <- data.table::fread(snakemake@input$tissue_summary)
results$tissue <- factor(results$tissue, levels = results$tissue)
output_file <- snakemake@output$go_enrichment_plot

p <- ggplot(results, aes(x = tissue, y = n_enriched)) +
    geom_col(fill = "#287D8E") +
    geom_text(aes(label = n_enriched), vjust = -0.3) +
    theme_minimal(base_size = 11) +
    labs(
        x = "Tissue",
        y = "Number of enriched communities",
        title = "GO-enriched communities per tissue",
        subtitle = "Communities with at least one significant BP term"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
ggsave(output_file, p, width = 12, height = 6)
