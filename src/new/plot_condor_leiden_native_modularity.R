## --------- ##
## Libraries ##
## --------- ##
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
options(stringsAsFactors = FALSE)

## -------------------------- ##
## Read native modularities ##
## -------------------------- ##
tissues <- as.character(unlist(snakemake@params$tissues))
resolutions <- as.character(unlist(snakemake@params$resolutions))
benchmark <- fread(snakemake@input$filtering_benchmark)
required <- c("Network", "Modularity", "tissue_type", "resolution")
if (!all(required %in% names(benchmark))) stop("Missing filtering benchmark columns")

# Reuse the recorded CONDOR result; do not recompute or substitute selected GMTs.
condor <- benchmark[Network == "Prior filtered", .(
    tissue = tissue_type,
    resolution = as.numeric(resolution),
    method = "CONDOR: bipartite PANDA",
    modularity = as.numeric(Modularity)
)]
leiden <- rbindlist(lapply(as.character(snakemake@input$leiden_results), function(path) {
    data <- fread(file.path(path, "summary.tsv"), na.strings = c("", "NA", "NaN"))
    if (!all(c("tissue", "filtering", "method", "leiden_resolution", "modularity") %in% names(data))) {
        stop("Missing Leiden summary columns in ", path)
    }
    data <- data[method == "leiden"]
    if (!nrow(data) || any(data$filtering != "prior")) stop("Expected prior-filtered Leiden results in ", path)
    data[, .(
        tissue, resolution = as.numeric(leiden_resolution),
        method = "Leiden: gene projection", modularity = as.numeric(modularity)
    )]
}))
data <- rbindlist(list(condor, leiden))
data <- data[tissue %in% tissues & resolution %in% as.numeric(resolutions)]
if (anyDuplicated(data[, .(tissue, resolution, method)])) stop("Duplicate modularity cells")
if (nrow(data) != length(tissues) * length(resolutions) * 2L) stop("Missing native modularity cells")
data$modularity[!is.finite(data$modularity)] <- NA_real_
data$cell_label <- ifelse(is.na(data$modularity), "NA", "")
data$tissue <- factor(data$tissue, levels = rev(tissues))
data$resolution <- factor(as.character(data$resolution), levels = resolutions)
data$method <- factor(data$method, levels = c("CONDOR: bipartite PANDA", "Leiden: gene projection"))

## --------------------- ##
## Draw and save tiles ##
## --------------------- ##
plot <- ggplot(data, aes(x = resolution, y = tissue, fill = modularity)) +
    geom_tile() +
    geom_text(aes(label = cell_label), size = 2.5, colour = "grey35") +
    facet_wrap(~ method, nrow = 1) +
    scale_fill_gradient(low = "white", high = "steelblue", na.value = "grey85", name = "Modularity") +
    labs(
        title = "Prior-filtered CONDOR and Leiden native modularity",
        x = "Algorithm resolution", y = "Tissue",
        caption = paste(
            "CONDOR: stored bipartite modularity from the Prior filtered benchmark at each resolution.",
            "\nLeiden: CPM detection at each resolution, scored with weighted projection modularity at gamma = 1.",
            "\nDifferent graphs and null models: scores are not a common-objective ranking. Grey NA: undefined scores."
        )
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12), strip.text = element_text(size = 12),
        panel.grid = element_blank(), plot.caption = element_text(hjust = 0)
    )
output_file <- snakemake@output$heatmap
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
ggsave(output_file, plot, width = 10, height = max(5, length(tissues) * .22 + 2))
