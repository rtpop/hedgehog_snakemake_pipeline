plot_n_sig_terms <- function(input_file, output_file) {
    results <- data.table::fread(input_file)
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    if (!nrow(results)) {
        p <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
            label = "No selected communities") + ggplot2::theme_void()
    } else {
        p <- ggplot2::ggplot(results, ggplot2::aes(
            x = reorder(community, -n_sig_terms), y = n_sig_terms)) +
            ggplot2::geom_col(fill = "#287D8E") +
            ggplot2::theme_minimal(base_size = 11) +
            ggplot2::labs(x = "Selected community", y = "Significant BP terms",
                title = "GO enrichment of selected communities",
                subtitle = "BH adjustment over all tested terms within each community") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    }
    ggplot2::ggsave(output_file, p, width = 12, height = 6)
}
