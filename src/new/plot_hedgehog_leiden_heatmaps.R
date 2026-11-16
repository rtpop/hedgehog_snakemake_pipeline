## --------- ##
## Libraries ##
## --------- ##
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

options(stringsAsFactors = FALSE)

## ----------------------- ##
## Read comparison results ##
## ----------------------- ##
comparison_dirs <- as.character(snakemake@input$comparisons)
resource_logs <- as.character(unlist(snakemake@params$resource_logs))
tissues <- as.character(unlist(snakemake@params$tissues))
resolutions <- as.character(unlist(snakemake@params$resolutions))

results <- rbindlist(lapply(comparison_dirs, function(path) {
    summary <- fread(file.path(path, "summary.tsv"), na.strings = c("", "NA", "NaN"))
    if (!"whole_partition_modularity" %in% names(summary)) {
        stop("Comparison summary predates whole-partition scoring; regenerate comparison: ", path)
    }
    scores <- fread(file.path(path, "community_scores.tsv"), na.strings = c("", "NA", "NaN"))
    if (!all(c("method", "modularity_contribution") %in% names(scores))) {
        stop("Missing per-community modularity columns in ", path)
    }
    # Use the same arithmetic mean of individual community contributions for
    # both methods. Do not substitute Leiden's total partition modularity.
    # An empty method or any undefined contribution produces an NA cell.
    summary$mean_modularity_contribution <- vapply(summary$method, function(method_name) {
        values <- scores[method == method_name, modularity_contribution]
        if (!length(values) || any(!is.finite(values))) return(NA_real_)
        mean(values)
    }, numeric(1))
    summary
}))
required <- c("tissue", "method", "leiden_resolution", "mean_modularity_contribution",
              "whole_partition_modularity", "algorithm_seconds")
if (!all(required %in% names(results))) stop("Missing comparison summary columns")
if (length(resource_logs) != length(tissues)) stop("Expected one BiHiDeF log per tissue")

## ---------------------------- ##
## Read HEDGEHOG elapsed times ##
## ---------------------------- ##
# Logs measure the full BiHiDeF call, including its resolution search and output.
# They do not measure filtering or the subsequent selected-community GMT step.
# These are optional existing files, not workflow dependencies: plotting must
# never request a BiHiDeF run just to obtain a missing runtime measurement.
read_elapsed <- function(path) {
    if (!file.exists(path)) {
        message("No existing BiHiDeF runtime log; showing NA: ", path)
        return(NA_real_)
    }
    lines <- readLines(path, warn = FALSE)
    elapsed <- grep("^Elapsed time:", lines, value = TRUE)
    status <- grep("^Exit status:", lines, value = TRUE)
    if (length(status) != 1L || trimws(sub("^Exit status:", "", status)) != "0") {
        message("BiHiDeF log does not record a successful run; showing NA: ", path)
        return(NA_real_)
    }
    if (length(elapsed) != 1L) {
        message("Missing or duplicate BiHiDeF elapsed time; showing NA: ", path)
        return(NA_real_)
    }
    seconds <- suppressWarnings(as.numeric(sub(
        "^Elapsed time: ([0-9.eE+-]+) seconds$", "\\1", elapsed
    )))
    if (!is.finite(seconds) || seconds < 0) {
        message("Invalid BiHiDeF elapsed time; showing NA: ", path)
        return(NA_real_)
    }
    seconds
}
hedgehog_seconds <- setNames(vapply(resource_logs, read_elapsed, numeric(1)), tissues)
results[method == "hedgehog", algorithm_seconds := unname(hedgehog_seconds[tissue])]

# Each method must have exactly one cell per tissue and Leiden resolution.
if (anyDuplicated(results[, .(tissue, method, leiden_resolution)])) {
    stop("Duplicate tissue/method/resolution comparison rows")
}
results$tissue <- factor(results$tissue, levels = rev(tissues))
results$resolution <- factor(as.character(results$leiden_resolution), levels = resolutions)
results$method <- factor(results$method, levels = c("hedgehog", "leiden"),
                         labels = c("HEDGEHOG", "Leiden"))

## --------------------- ##
## Draw method heatmaps ##
## --------------------- ##
# Use the benchmark heatmap's white-to-steelblue tiles and method facets.
# Keep negative modularity values rather than clipping them to zero. Grey cells
# explicitly mark undefined/unavailable scores; they are never filled with zero.
plot_heatmap <- function(metric, title, legend, caption, output_file) {
    data <- copy(results)
    data$value <- as.numeric(data[[metric]])
    data$value[!is.finite(data$value)] <- NA_real_
    data$cell_label <- ifelse(is.na(data$value), "NA", "")
    plot <- ggplot(data, aes(x = resolution, y = tissue, fill = value)) +
        geom_tile() +
        geom_text(aes(label = cell_label), size = 2.5, colour = "grey35") +
        facet_wrap(~ method, nrow = 1, drop = FALSE) +
        scale_x_discrete(drop = FALSE) +
        scale_y_discrete(drop = FALSE) +
        scale_fill_gradient(low = "white", high = "steelblue", na.value = "grey85",
                            name = legend) +
        labs(title = title, x = "Leiden resolution", y = "Tissue", caption = caption) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_text(size = 12),
            legend.title = element_text(size = 12),
            strip.text = element_text(size = 12),
            panel.grid = element_blank(),
            plot.caption = element_text(hjust = 0)
        )
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    ggsave(output_file, plot, width = 10, height = max(5, length(tissues) * .22 + 2))
}

plot_heatmap(
    "mean_modularity_contribution", "HEDGEHOG and Leiden community modularity",
    "Mean modularity\ncontribution",
    paste(
        "Arithmetic mean of per-community weighted modularity contributions (gamma = 1).",
        "\nThis community score is not total partition modularity; community counts and sizes affect it.",
        "\nHEDGEHOG uses the same selected communities across columns. Grey NA: empty or undefined scores."
    ),
    snakemake@output$modularity
)

plot_heatmap(
    "whole_partition_modularity", "HEDGEHOG and Leiden whole-partition modularity",
    "Whole-partition\nmodularity",
    paste(
        "Fractional overlapping modularity (gamma = 1): equal membership across each gene's communities.",
        "\nUncovered genes receive distinct singletons; for Leiden this equals standard partition modularity.",
        "\nHEDGEHOG uses the same selected communities across columns. Grey NA: undefined scores."
    ),
    snakemake@output$partition_modularity
)

plot_heatmap(
    "algorithm_seconds", "HEDGEHOG and Leiden runtime", "Seconds",
    paste(
        "HEDGEHOG: full BiHiDeF call, repeated across columns; selection and filtering excluded.",
        "\nLeiden: partition detection at each resolution; projection and scoring excluded.",
        "\nThese timings measure different scopes and are not an equal-work speed comparison.",
        "\nGrey NA cells indicate unavailable or invalid recorded timings."
    ),
    snakemake@output$runtime
)
