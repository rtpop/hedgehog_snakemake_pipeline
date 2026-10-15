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

    comparisons$pair_id <- paste(
        tissue_a,
        comparisons$community1,
        tissue_b,
        comparisons$community2,
        sep = ":"
    )
    comparisons$pair_id <- factor(comparisons$pair_id, levels = comparisons$pair_id)

    plot <- ggplot2::ggplot(
        comparisons,
        ggplot2::aes(x = 1, y = pair_id, fill = jaccard_index)
    ) +
        ggplot2::geom_tile(color = NA) +
        ggplot2::scale_fill_gradient(
            low = "white",
            high = "darkred",
            limits = c(0, 1)
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank()
        ) +
        ggplot2::labs(
            title = paste("Jaccard per BiHiDeF community pairing for", tissue_a, "and", tissue_b),
            subtitle = paste(
                "Pairs with Jaccard > 0.5:",
                sum(comparisons$jaccard_index > 0.5, na.rm = TRUE)
            ),
            fill = "Jaccard"
        )

    output_file <- file.path(
        output_dir,
        paste0(tissue_a, "_", tissue_b, "_jaccard_heatmap.pdf")
    )
    ggplot2::ggsave(output_file, plot, width = 10, height = 10, dpi = 300)
}
