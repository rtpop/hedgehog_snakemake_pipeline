## --------- ##
## Libraries ##
## --------- ##
suppressPackageStartupMessages(library(ggplot2))

options(stringsAsFactors = FALSE)

## -------------------- ##
## Workflow parameters ##
## -------------------- ##
community_files <- as.character(snakemake@input$selected_communities)
tissues <- as.character(unlist(snakemake@params$tissues))
min_jaccard <- as.numeric(snakemake@params$min_jaccard)
output_file <- snakemake@output$jaccard_boxplot

if (length(tissues) < 2L || anyDuplicated(tissues)) {
    stop("Configure at least two distinct tissues for the Jaccard boxplot")
}
if (length(community_files) != length(tissues)) {
    stop("Expected one selected-community GMT file per tissue")
}
if (length(min_jaccard) != 1L || !is.finite(min_jaccard) ||
    min_jaccard < 0 || min_jaccard > 1) {
    stop("jaccard_boxplot_min_jaccard must be between 0 and 1")
}

## ------------------------ ##
## Read selected gene sets ##
## ------------------------ ##
# GMT columns are community ID, description, then gene symbols. The description
# must not be included when calculating shared genes or the union of gene sets.
read_communities <- function(path) {
    lines <- readLines(path, warn = FALSE)
    sets <- list()
    for (line in lines[nzchar(trimws(lines))]) {
        fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
        if (length(fields) < 3L || !nzchar(fields[1])) {
            stop("Invalid selected-community GMT record in ", path)
        }
        if (fields[1] %in% names(sets)) {
            stop("Duplicate community ID in ", path, ": ", fields[1])
        }
        genes <- unique(trimws(fields[-c(1, 2)]))
        genes <- genes[nzchar(genes)]
        if (!length(genes)) stop("Empty gene set in ", path, ": ", fields[1])
        sets[[fields[1]]] <- genes
    }
    sets
}

communities <- setNames(lapply(community_files, read_communities), tissues)

## ----------------------------- ##
## Compare cross-tissue pairs   ##
## ----------------------------- ##
# Compare every community in one tissue with every community in the other.
# Keep the original script's inclusive Jaccard cutoff before computing boxes.
# These are final BiHiDeF sets, not resolution-specific CONDOR partitions.
tissue_pairs <- utils::combn(sort(tissues), 2, simplify = FALSE)
pair_labels <- vapply(tissue_pairs, paste, character(1), collapse = " vs ")
comparisons <- lapply(seq_along(tissue_pairs), function(index) {
    pair <- tissue_pairs[[index]]
    sets_a <- communities[[pair[1]]]
    sets_b <- communities[[pair[2]]]
    values <- unlist(lapply(sets_a, function(genes_a) {
        vapply(sets_b, function(genes_b) {
            length(intersect(genes_a, genes_b)) / length(union(genes_a, genes_b))
        }, numeric(1))
    }), use.names = FALSE)
    retained <- values[values >= min_jaccard]
    data.frame(
        pair = rep(pair_labels[index], length(retained)),
        jaccard_index = as.numeric(retained)
    )
})
df_pairs_filt <- do.call(rbind, comparisons)
df_pairs_filt$pair <- factor(df_pairs_filt$pair, levels = pair_labels)

## -------------------- ##
## Draw and save boxes ##
## -------------------- ##
# Preserve the original boxplot styling and dimensions. Retain empty tissue-pair
# categories so a pair with no overlaps above the cutoff is still represented.
p_box <- ggplot(df_pairs_filt, aes(x = pair, y = jaccard_index)) +
    geom_boxplot(fill = "grey80", color = "black") +
    theme_minimal() +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
        title = sprintf("Jaccard boxplots (J \u2265 %s)", min_jaccard),
        x = "Tissue pair",
        y = "Jaccard index"
    )

if (!nrow(df_pairs_filt)) {
    p_box <- p_box + annotate(
        "text", x = 1, y = 0.5,
        label = "No community pairs meet the Jaccard cutoff"
    )
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
ggsave(output_file, p_box, width = 7, height = 4, dpi = 300)
