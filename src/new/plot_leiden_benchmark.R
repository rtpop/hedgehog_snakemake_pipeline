suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

df <- fread(snakemake@input$summary)
required_columns <- c("tissue", "method", "resolution", "modularity")
if (!all(required_columns %in% names(df))) {
    stop("Leiden benchmark summary must contain: ", paste(required_columns, collapse = ", "))
}

# Plot the saved benchmark scores directly; each observation is one tissue.
df <- df[method %in% c("condor", "leiden")]
resolutions <- sort(unique(df$resolution))
omitted <- sum(!is.finite(df$modularity))
df <- df[is.finite(modularity)]
if (!nrow(df)) stop("No finite modularity scores available for the Leiden benchmark boxplot")
df[, resolution := factor(resolution, levels = resolutions)]
df[, method := factor(method, levels = c("condor", "leiden"), labels = c("CONDOR", "Leiden"))]

p <- ggplot(df, aes(x = resolution, y = modularity, fill = method)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.65,
                 outlier.shape = NA, alpha = 0.65) +
    geom_point(aes(group = method), shape = 21, size = 1.5, alpha = 0.65,
               position = position_jitterdodge(jitter.width = 0.12,
                                             dodge.width = 0.8, seed = 0)) +
    scale_fill_manual(values = c("CONDOR" = "steelblue", "Leiden" = "darkorange"), drop = FALSE) +
    labs(title = "CONDOR and Leiden benchmark",
         subtitle = "Modularity on the gene projection across tissues",
         x = "Resolution", y = "Modularity", fill = "Method",
         caption = paste0("Each point represents one tissue.",
                          if (omitted > 0) paste0(" Missing or nonfinite scores omitted: ", omitted, ".") else "")) +
    theme_minimal() +
    theme(legend.position = "top")

output <- snakemake@output$plot
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
ggsave(output, plot = p, width = 8, height = 5)
