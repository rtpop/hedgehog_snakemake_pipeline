# Human gene symbols, BP ontology, and default classic/Fisher testing are
# preserved from the old workflow. Background files use column 2 as before.
read_gmt <- function(path) {
    lines <- readLines(path, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    sets <- list()
    for (line in lines) {
        fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
        if (length(fields) < 3L || !nzchar(fields[1])) {
            stop("Invalid GMT record: expected name, description, and genes")
        }
        if (fields[1] %in% names(sets)) stop("Duplicate GMT community: ", fields[1])
        genes <- unique(trimws(fields[-c(1, 2)]))
        sets[[fields[1]]] <- genes[nzchar(genes)]
    }
    sets
}

empty_go_summary <- function() {
    data.frame(community = character(), n_sig_terms = integer(),
               min_adj_p = numeric(), top_go_id = character(),
               top_go_term = character(), status = character())
}

enrich_community <- function(gene_set, background, mapping, comm_name,
                             sig_thresh, statistic, algorithm) {
    summary <- data.frame(community = comm_name, n_sig_terms = 0L,
                          min_adj_p = NA_real_, top_go_id = NA_character_,
                          top_go_term = NA_character_, status = "no_annotated_genes")
    terms <- data.frame(GO.ID = character(), Term = character(),
                        p_value = numeric(), p_adj = numeric())
    if (!any(intersect(gene_set, background) %in% names(mapping))) {
        return(list(summary = summary, terms = terms))
    }
    selection <- factor(as.integer(background %in% gene_set), levels = c(0, 1))
    names(selection) <- background
    go_object <- methods::new("topGOdata", ontology = "BP", allGenes = selection,
                              annot = topGO::annFUN.gene2GO, gene2GO = mapping)
    if (!length(topGO::usedGO(go_object))) {
        summary$status <- "no_testable_terms"
        return(list(summary = summary, terms = terms))
    }
    result <- topGO::runTest(go_object, algorithm = algorithm, statistic = statistic)
    # Adjust raw numeric scores over ALL tested terms, not GenTable's formatted
    # strings or its default top-ten subset. Adjustment is per community.
    raw_p <- topGO::score(result)
    adjusted <- stats::p.adjust(raw_p, method = "BH")
    names(adjusted) <- names(raw_p)
    terms <- topGO::GenTable(go_object, test = result, orderBy = "test",
                             ranksOf = "test", topNodes = length(raw_p))
    terms$p_value <- unname(raw_p[terms$GO.ID])
    terms$p_adj <- unname(adjusted[terms$GO.ID])
    terms <- terms[order(terms$p_adj, terms$p_value, terms$GO.ID), , drop = FALSE]
    summary$n_sig_terms <- sum(terms$p_adj < sig_thresh, na.rm = TRUE)
    if (any(is.finite(terms$p_adj))) {
        summary$min_adj_p <- terms$p_adj[1]
        summary$top_go_id <- terms$GO.ID[1]
        summary$top_go_term <- terms$Term[1]
        summary$status <- "tested"
    } else {
        summary$status <- "no_finite_p_values"
    }
    list(summary = summary, terms = terms)
}

run_go_enrichment <- function(gmt_file, summary_file, details_dir,
                              auto_bg = TRUE, bg_file = NULL, save_all = TRUE,
                              sig_thresh = 0.05, statistic = "fisher", algorithm = "classic") {
    if (!is.finite(sig_thresh) || sig_thresh <= 0 || sig_thresh > 1) {
        stop("sig_thresh must be in (0, 1]")
    }
    gene_sets <- read_gmt(gmt_file)
    dir.create(dirname(summary_file), recursive = TRUE, showWarnings = FALSE)
    dir.create(details_dir, recursive = TRUE, showWarnings = FALSE)
    if (!length(gene_sets)) {
        data.table::fwrite(empty_go_summary(), summary_file, sep = "\t", na = "NA")
        return(invisible(empty_go_summary()))
    }
    if (isTRUE(auto_bg)) {
        background <- unique(unlist(gene_sets, use.names = FALSE))
    } else {
        if (is.null(bg_file)) stop("bg_file is required when auto_bg is FALSE")
        annotation <- data.table::fread(bg_file, header = FALSE, colClasses = "character")
        if (ncol(annotation) < 2L) stop("Background file must have gene symbols in column 2")
        background <- unique(trimws(annotation[[2]]))
        background <- background[!is.na(background) & nzchar(background)]
    }
    if (!length(background)) stop("The gene background is empty")
    go_mapping <- topGO::annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
    gene_mapping <- topGO::inverseList(go_mapping)
    if (!any(background %in% names(gene_mapping))) {
        stop("No background genes map to human BP annotations; expected human gene symbols")
    }
    summaries <- lapply(names(gene_sets), function(name) {
        result <- enrich_community(gene_sets[[name]], background, gene_mapping,
                                    name, sig_thresh, statistic, algorithm)
        if (isTRUE(save_all)) {
            # Keep normal HiDeF IDs unchanged, and prevent IDs becoming paths.
            if (grepl("[/\\\\]", name) || name %in% c(".", "..")) stop("Invalid community ID: ", name)
            data.table::fwrite(result$terms, file.path(details_dir, paste0(name, ".txt")),
                               sep = "\t", na = "NA")
        }
        result$summary
    })
    summary <- do.call(rbind, summaries)
    data.table::fwrite(summary, summary_file, sep = "\t", na = "NA")
    invisible(summary)
}
