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

plot_enriched_tf_qscores <- function(result, output, tissue, top_n) {
    grDevices::pdf(output, width = 10, height = max(6, 2 + min(top_n, 40) * .25), onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::plot.new()
    graphics::title(main = paste("Top TF qscores in GO-enriched communities -", tissue))
    graphics::text(.05, .85, adj = c(0, 1), labels = paste0(
        "GO-enriched selected communities: ", result$n_enriched, "\n",
        "Requested TFs per community / resolution: ", top_n, "\n\n",
        "Enriched means at least one significant GO term in the existing GO summary.\n",
        "TFs are ranked by descending original qscore, separately at each exact target-set match.\n",
        "Ties are resolved by TF identifier; undefined qscores are excluded.\n",
        "Unmatched communities and missing scores are listed in the coverage TSV."))
    for (i in which(result$coverage$status == "plotted")) {
        mapping <- result$coverage[i, ]
        rows <- result$top[result$top$selected_community == mapping$selected_community &
                           result$top$condor_community == mapping$condor_community &
                           result$top$resolution == mapping$resolution, , drop = FALSE]
        rows <- rows[order(rows$rank, decreasing = TRUE), ]
        graphics::par(mar = c(5, max(8, min(20, max(nchar(rows$reg)) * .6)), 5, 2))
        graphics::barplot(rows$qscore, names.arg = rows$reg, horiz = TRUE, las = 1,
                          col = "steelblue", border = NA, xlab = "TF qscore",
                          cex.names = min(1, 30 / nrow(rows)),
                          main = paste(tissue, mapping$selected_community, sep = " | "))
        graphics::abline(v = 0, col = "grey50")
        graphics::mtext(paste("CONDOR community", mapping$condor_community,
                             "| resolution", mapping$resolution), side = 3, line = .5)
    }
}
