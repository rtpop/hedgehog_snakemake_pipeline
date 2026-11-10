required_libraries <- c("dplyr", "ggplot2")

for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}

options(stringsAsFactors = FALSE)

community_files <- as.character(snakemake@input$selected_communities)
output_dir <- snakemake@output$jaccard_heatmaps
tissues <- as.character(unlist(snakemake@params$tissues))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_communities <- function(file, tissue) {
    if (!file.exists(file)) {
        stop("Missing BiHiDef community file: ", file)
    }

    lines <- readLines(file, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]

    if (length(lines) == 0) {
        return(data.frame())
    }

    communities <- lapply(lines, function(line) {
        fields <- strsplit(trimws(line), "\\s+")[[1]]
        if (length(fields) < 2) {
            return(NULL)
        }

        gene_start <- 2
        if (
            length(fields) > 2 &&
            (grepl("^[0-9]+$", fields[[2]]) || tolower(fields[[2]]) == "genes")
        ) {
            gene_start <- 3
        }

        genes <- unique(fields[gene_start:length(fields)])
        genes <- genes[nzchar(genes)]
        if (length(genes) == 0) {
            return(NULL)
        }

        data.frame(
            tissue = tissue,
            community = fields[[1]],
            genes = I(list(genes)),
            stringsAsFactors = FALSE
        )
    })

    dplyr::bind_rows(communities)
}

community_data <- dplyr::bind_rows(
    lapply(seq_along(community_files), function(index) {
        read_communities(community_files[[index]], tissues[[index]])
    })
)

if (nrow(community_data) == 0) {
    stop("No BiHiDef communities were found.")
}

calculate_jaccard <- function(genes_a, genes_b) {
    union_size <- length(union(genes_a, genes_b))
    if (union_size == 0) {
        return(0)
    }
    length(intersect(genes_a, genes_b)) / union_size
}

for (tissue_pair in utils::combn(tissues, 2, simplify = FALSE)) {
    tissue_a <- tissue_pair[[1]]
    tissue_b <- tissue_pair[[2]]
    communities_a <- community_data[community_data$tissue == tissue_a, ]
    communities_b <- community_data[community_data$tissue == tissue_b, ]

    if (nrow(communities_a) == 0 || nrow(communities_b) == 0) {
        message("No communities found for ", tissue_a, " and ", tissue_b)
        next
    }

    comparisons <- do.call(rbind, lapply(seq_len(nrow(communities_a)), function(index_a) {
        do.call(rbind, lapply(seq_len(nrow(communities_b)), function(index_b) {
            jaccard_index <- calculate_jaccard(
                communities_a$genes[[index_a]],
                communities_b$genes[[index_b]]
            )
            data.frame(
                community1 = communities_a$community[[index_a]],
                community2 = communities_b$community[[index_b]],
                jaccard_index = jaccard_index,
                stringsAsFactors = FALSE
            )
        }))
    }))

    comparisons$community1 <- factor(
        comparisons$community1, levels = rev(unique(communities_a$community))
    )
    comparisons$community2 <- factor(
        comparisons$community2, levels = unique(communities_b$community)
    )

    output_file <- file.path(
        output_dir,
        paste0(tissue_a, "_", tissue_b, "_jaccard_heatmap.pdf")
    )
    # Split large matrices across pages so every community ID stays readable.
    row_page <- (nlevels(comparisons$community1) - as.integer(comparisons$community1)) %/% 30L
    col_page <- (as.integer(comparisons$community2) - 1L) %/% 30L
    pages <- split(comparisons, list(col_page, row_page), drop = TRUE)

    grDevices::pdf(output_file, width = 12, height = 11)
    tryCatch({
        for (page_number in seq_along(pages)) {
            plot <- ggplot2::ggplot(
                pages[[page_number]],
                ggplot2::aes(x = community2, y = community1, fill = jaccard_index)
            ) +
                ggplot2::geom_tile(color = "grey85", linewidth = 0.2) +
                ggplot2::scale_fill_gradient(
                    low = "white",
                    high = "darkred",
                    limits = c(0, 1),
                    breaks = seq(0, 1, by = 0.25)
                ) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                    panel.grid = ggplot2::element_blank()
                ) +
                ggplot2::labs(
                    title = paste("BiHiDeF community overlap |", tissue_a, "vs", tissue_b),
                    subtitle = sprintf(
                        "%d of %d pairs have Jaccard > 0.5 across all pages | Page %d of %d",
                        sum(comparisons$jaccard_index > 0.5),
                        nrow(comparisons), page_number, length(pages)
                    ),
                    x = paste(tissue_b, "community ID"),
                    y = paste(tissue_a, "community ID"),
                    fill = "Jaccard index",
                    caption = paste(
                        "Shared genes / genes in either community.",
                        "0 = no shared genes. 1 = identical gene sets."
                    )
                )
            print(plot)
        }
    }, finally = grDevices::dev.off())
}
