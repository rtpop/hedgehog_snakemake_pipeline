rank_enriched_tf_qscores <- function(go, matches, regulators, top_n) {
    if (length(top_n) != 1L || !is.numeric(top_n) || !is.finite(top_n) ||
        top_n < 1 || top_n != floor(top_n)) stop("top_n must be a positive integer")
    if (!all(c("community", "n_sig_terms") %in% names(go))) stop("Missing GO summary columns")
    if (!all(c("selected_community", "condor_community", "resolution") %in% names(matches))) {
        stop("Missing matched-target columns")
    }
    if (!all(c("reg", "community", "resolution", "qscore") %in% names(regulators))) {
        stop("Missing regulator qscore columns: reg, community, resolution, qscore required")
    }
    for (field in c("resolution", "qscore")) {
        text <- trimws(as.character(regulators[[field]]))
        value <- suppressWarnings(as.numeric(text))
        missing <- is.na(text) | tolower(text) %in% c("", "na", "nan")
        if (any(!missing & !is.finite(value))) stop("Invalid regulator ", field)
        value[missing] <- NA_real_
        regulators[[field]] <- value
    }
    matches$resolution <- as.numeric(matches$resolution)
    if (any(!is.finite(regulators$resolution) | regulators$resolution <= 0) ||
        any(!is.finite(matches$resolution) | matches$resolution <= 0)) stop("Invalid resolution")
    if (anyNA(regulators$reg) || anyNA(regulators$community) ||
        any(!nzchar(regulators$reg)) || any(!nzchar(regulators$community))) stop("Missing TF or community ID")
    if (anyDuplicated(regulators[c("resolution", "reg")])) stop("Duplicate TF at a resolution")
    n_terms <- suppressWarnings(as.numeric(go$n_sig_terms))
    enriched <- unique(as.character(go$community[is.finite(n_terms) & n_terms > 0]))
    mapped <- unique(matches[matches$selected_community %in% enriched,
                             c("selected_community", "condor_community", "resolution"), drop = FALSE])
    joined <- merge(mapped, regulators, by.x = c("condor_community", "resolution"),
                    by.y = c("community", "resolution"), sort = FALSE)
    top <- data.frame(selected_community = character(), condor_community = character(),
                      resolution = numeric(), reg = character(), qscore = numeric(), rank = integer())
    coverage <- data.frame(selected_community = character(), condor_community = character(),
                           resolution = numeric(), n_tfs = integer(), n_undefined = integer(),
                           n_plotted = integer(), status = character())
    for (community in enriched) {
        mappings <- mapped[mapped$selected_community == community, , drop = FALSE]
        if (!nrow(mappings)) {
            coverage <- rbind(coverage, data.frame(selected_community = community,
                condor_community = NA_character_, resolution = NA_real_, n_tfs = 0L,
                n_undefined = 0L, n_plotted = 0L, status = "no_exact_target_match"))
        }
        for (i in seq_len(nrow(mappings))) {
            mapping <- mappings[i, ]
            rows <- joined[joined$selected_community == community &
                           joined$condor_community == mapping$condor_community &
                           joined$resolution == mapping$resolution, , drop = FALSE]
            finite <- rows[is.finite(rows$qscore), , drop = FALSE]
            finite <- finite[order(-finite$qscore, finite$reg), , drop = FALSE]
            finite <- head(finite, top_n)
            if (nrow(finite)) {
                finite$rank <- seq_len(nrow(finite))
                top <- rbind(top, finite[, names(top), drop = FALSE])
            }
            coverage <- rbind(coverage, data.frame(selected_community = community,
                condor_community = mapping$condor_community, resolution = mapping$resolution,
                n_tfs = nrow(rows), n_undefined = sum(!is.finite(rows$qscore)), n_plotted = nrow(finite),
                status = if (!nrow(rows)) "no_regulator_scores" else if (!nrow(finite)) "all_scores_undefined" else "plotted"))
        }
    }
    list(top = top, coverage = coverage, n_enriched = length(enriched))
}

read_enriched_go_terms <- function(go, details_dir, sig_thresh) {
    terms <- data.frame(selected_community = character(), GO.ID = character(), Term = character(),
                         p_value = numeric(), p_adj = numeric())
    communities <- as.character(go$community[!is.na(go$n_sig_terms) & as.numeric(go$n_sig_terms) > 0])
    for (community in communities) {
        path <- file.path(details_dir, paste0(community, ".txt"))
        if (!file.exists(path)) stop("Missing GO term details for ", community,
                                     "; generate GO enrichment with save_all: true")
        rows <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
        if (!all(c("GO.ID", "Term", "p_value", "p_adj") %in% names(rows))) {
            stop("Missing GO term columns in ", path)
        }
        rows <- rows[is.finite(rows$p_adj) & rows$p_adj < sig_thresh,
                     c("GO.ID", "Term", "p_value", "p_adj"), drop = FALSE]
        rows <- rows[order(rows$p_adj, rows$p_value, rows$GO.ID), , drop = FALSE]
        if (nrow(rows)) {
            rows$selected_community <- community
            terms <- rbind(terms, rows[, names(terms), drop = FALSE])
        }
    }
    terms
}

plot_enriched_tf_qscores <- function(result, output, tissue, top_n, go_terms = NULL) {
    plotted <- which(result$coverage$status == "plotted")
    n_panels <- length(plotted)
    n_columns <- max(1L, ceiling(sqrt(n_panels)))
    n_rows <- max(1L, ceiling(n_panels / n_columns))
    # Grow the single page with the panel count to keep TF labels readable.
    term_labels <- lapply(plotted, function(i) {
        community <- result$coverage$selected_community[i]
        terms <- go_terms[go_terms$selected_community == community, , drop = FALSE]
        if (is.null(terms) || !nrow(terms)) return("No significant GO term details available")
        unlist(lapply(seq_len(nrow(terms)), function(j) {
            strwrap(paste0(terms$GO.ID[j], "  ", terms$Term[j],
                          " (BH-adjusted p = ", format(terms$p_adj[j], digits = 3), ")"), width = 55)
        }), use.names = FALSE)
    })
    panel_height <- max(4, 2 + top_n * .25, 1.5 + max(c(0, lengths(term_labels))) * .17)
    grDevices::pdf(output, width = 11 * n_columns, height = panel_height * n_rows + .6,
                   onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
    panels <- matrix(0L, nrow = n_rows, ncol = 2 * n_columns)
    if (n_panels) {
        for (k in seq_len(n_panels)) {
            row <- (k - 1L) %/% n_columns + 1L
            col <- 2 * ((k - 1L) %% n_columns) + 1L
            panels[row, col:(col + 1L)] <- c(2L * k - 1L, 2L * k)
        }
    } else panels[1, 1:2] <- 1L
    graphics::layout(panels, widths = rep(c(6, 5), n_columns))
    graphics::par(oma = c(0, 0, 3, 0), cex = 1)
    if (!n_panels) {
        graphics::plot.new()
        graphics::text(.5, .5, if (!result$n_enriched) "No GO-enriched communities" else
                       "No matched finite TF qscores; see the coverage TSV")
    }
    for (k in seq_along(plotted)) {
        i <- plotted[k]
        mapping <- result$coverage[i, ]
        rows <- result$top[result$top$selected_community == mapping$selected_community &
                           result$top$condor_community == mapping$condor_community &
                           result$top$resolution == mapping$resolution, , drop = FALSE]
        rows <- rows[order(rows$rank, decreasing = TRUE), ]
        graphics::par(mar = c(5, max(8, min(20, max(nchar(rows$reg)) * .6)), 5, 2))
        graphics::barplot(rows$qscore, names.arg = rows$reg, horiz = TRUE, las = 1,
                          col = "steelblue", border = NA, xlab = "TF qscore",
                          cex.names = min(1, 30 / nrow(rows)),
                          main = paste("HEDGEHOG community", mapping$selected_community))
        graphics::abline(v = 0, col = "grey50")
        graphics::par(mar = c(1, 0, 3, 1))
        graphics::plot.new()
        graphics::title(main = "Significant GO terms", cex.main = 1)
        labels <- term_labels[[k]]
        graphics::text(rep(0, length(labels)), seq(1, 0, length.out = length(labels) + 1L)[seq_along(labels)],
                       labels = labels, adj = c(0, 1), cex = .8)
    }
    graphics::mtext(paste("Top TF qscores in GO-enriched communities -", tissue),
                    outer = TRUE, side = 3, line = 1, font = 2)
}
