read_selected_gmt <- function(path) {
    lines <- readLines(path, warn = FALSE)
    sets <- list()
    for (line in lines[nzchar(trimws(lines))]) {
        fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
        if (length(fields) < 3 || !nzchar(fields[1])) stop("Invalid GMT record")
        if (fields[1] %in% names(sets)) stop("Duplicate selected community: ", fields[1])
        genes <- unique(trimws(fields[-c(1, 2)]))
        genes <- genes[nzchar(genes)]
        if (!length(genes)) stop("Empty selected community: ", fields[1])
        sets[[fields[1]]] <- genes
    }
    sets
}

gene_signature <- function(genes) {
    genes <- sort(unique(genes))
    paste0(nchar(genes), ":", genes, collapse = "|")
}

match_selected_qscores <- function(sets, scores) {
    required <- c("tar", "community", "qscore", "resolution")
    if (!all(required %in% names(scores))) stop("Missing qscore CSV columns")
    scores$tar <- as.character(scores$tar)
    scores$community <- as.character(scores$community)
    if (anyNA(scores$tar) || anyNA(scores$community) ||
        any(!nzchar(scores$tar)) || any(!nzchar(scores$community))) stop("Missing gene or community ID")
    resolution <- suppressWarnings(as.numeric(as.character(scores$resolution)))
    if (any(!is.finite(resolution) | resolution <= 0)) stop("Invalid resolution")
    text <- trimws(as.character(scores$qscore))
    missing <- is.na(text) | tolower(text) %in% c("", "na", "nan")
    values <- suppressWarnings(as.numeric(text))
    if (any(!missing & !is.finite(values))) stop("Invalid or infinite qscore")
    values[missing] <- NA_real_
    scores$resolution <- resolution
    scores$qscore <- values
    if (anyDuplicated(scores[c("resolution", "tar")])) stop("Duplicate gene score at a resolution")
    signatures <- vapply(sets, gene_signature, character(1))
    # IDs alone never determine a match: compare every gene in the full group.
    groups <- split(seq_len(nrow(scores)), interaction(scores$resolution, scores$community,
                                                       drop = TRUE, lex.order = TRUE))
    mapped <- list()
    for (rows in groups) {
        signature <- gene_signature(scores$tar[rows])
        for (selected in names(signatures)[signatures == signature]) {
            mapped[[length(mapped) + 1L]] <- data.frame(
                selected_community = selected, condor_community = scores$community[rows],
                resolution = scores$resolution[rows], tar = scores$tar[rows],
                qscore = scores$qscore[rows], stringsAsFactors = FALSE)
        }
    }
    empty <- data.frame(selected_community = character(), condor_community = character(),
                        resolution = numeric(), tar = character(), qscore = numeric())
    mapped <- if (length(mapped)) do.call(rbind, mapped) else empty
    coverage <- data.frame(selected_community = as.character(names(sets)), n_genes = lengths(sets),
                            n_exact_matches = integer(length(sets)), n_score_rows = integer(length(sets)),
                            n_undefined = integer(length(sets)), status = rep("no_exact_match", length(sets)))
    for (i in seq_along(sets)) {
        rows <- mapped[mapped$selected_community == names(sets)[i], , drop = FALSE]
        coverage$n_exact_matches[i] <- nrow(unique(rows[c("resolution", "condor_community")]))
        coverage$n_score_rows[i] <- nrow(rows)
        coverage$n_undefined[i] <- sum(!is.finite(rows$qscore))
        if (nrow(rows)) coverage$status[i] <- "matched"
    }
    list(scores = mapped, coverage = coverage)
}

plot_selected_qscores <- function(result, output, tissue) {
    grDevices::pdf(output, width = 12, height = 6, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::plot.new()
    graphics::title(main = paste("Selected community qscores -", tissue))
    graphics::text(.05, .8, adj = c(0, 1), labels = paste0(
        "Selected communities: ", nrow(result$coverage), "\n",
        "Exact target-membership matches: ", sum(result$coverage$status == "matched"), "\n",
        "Unmatched communities: ", sum(result$coverage$status != "matched"), "\n\n",
        "Matches require identical complete gene sets, not identical numeric IDs.\n",
        "Scores retain original CONDOR regulator assignments and resolution.\n",
        "Histograms pool gene-resolution observations; genes can appear repeatedly.\n",
        "Undefined scores are counted and excluded from numeric plots.\n",
        "The current bihidef export omits the initial minimum resolution.\n",
        "The coverage CSV lists every selected community, including unmatched sets."))
    for (selected in result$coverage$selected_community[result$coverage$status == "matched"]) {
        rows <- result$scores[result$scores$selected_community == selected, , drop = FALSE]
        finite <- rows[is.finite(rows$qscore), , drop = FALSE]
        graphics::par(mfrow = c(1, 2), mar = c(5, 5, 4, 2), oma = c(0, 0, 3, 0))
        if (nrow(finite)) {
            graphics::hist(finite$qscore, breaks = min(40, max(1, ceiling(sqrt(nrow(finite))))),
                col = "#287D8E", border = "white", main = "Matched score distribution",
                xlab = "Qscore", ylab = "Gene-resolution observations")
            graphics::mtext(paste(sum(!is.finite(rows$qscore)), "undefined scores excluded"), side = 3, cex = .8)
            groups <- split(finite$qscore, finite$resolution)
            x <- as.numeric(names(groups))
            quantiles <- vapply(groups, stats::quantile, numeric(3), probs = c(.25, .5, .75), names = FALSE)
            graphics::plot(x, quantiles[2, ], log = "x", pch = 19, col = "#287D8E",
                ylim = range(quantiles), xlab = "CONDOR resolution", ylab = "Qscore",
                main = "Median and interquartile range")
            graphics::segments(x, quantiles[1, ], x, quantiles[3, ], col = "#287D8E", lwd = 2)
        } else {
            for (i in 1:2) {
                graphics::plot.new()
                graphics::text(.5, .5, "All matched scores are undefined")
            }
        }
        graphics::mtext(paste(tissue, selected, sep = " | "), outer = TRUE, font = 2)
    }
}

run_selected_qscores <- function(gmt, csv, mapped_file, coverage_file, plot_file, tissue) {
    sets <- read_selected_gmt(gmt)
    # Read IDs as character to preserve leading zeroes. Numeric columns are validated above.
    scores <- utils::read.csv(csv, colClasses = "character", check.names = FALSE,
                              na.strings = character(), fileEncoding = "UTF-8-BOM")
    result <- match_selected_qscores(sets, scores)
    for (path in c(mapped_file, coverage_file, plot_file)) {
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    }
    utils::write.csv(result$scores, mapped_file, row.names = FALSE, na = "NaN")
    utils::write.csv(result$coverage, coverage_file, row.names = FALSE, na = "NaN")
    plot_selected_qscores(result, plot_file, tissue)
    invisible(result)
}
