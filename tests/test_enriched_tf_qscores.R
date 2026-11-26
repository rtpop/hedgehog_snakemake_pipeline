# Run from repository root: Rscript tests/test_enriched_tf_qscores.R
source("src/new/enriched_tf_qscores_fn.R")
go <- data.frame(community = c("A", "B", "C"), n_sig_terms = c(2, 0, 1))
matches <- data.frame(selected_community = c("A", "A", "A", "B"),
                      condor_community = c("007", "007", "9", "007"), resolution = c(1, 1, 2, 1))
regulators <- data.frame(reg = c("TF1", "TF2", "TF3", "TF4", "TF1", "TF2"),
                         community = c("007", "007", "007", "007", "9", "9"),
                         resolution = c(1, 1, 1, 1, 2, 2), qscore = c(.5, .8, .8, NA, -.1, -.2))
result <- rank_enriched_tf_qscores(go, matches, regulators, 2)
stopifnot(nrow(result$top) == 4, identical(result$top$reg, c("TF2", "TF3", "TF1", "TF2")),
          all(result$top$selected_community == "A"), result$coverage$n_undefined[1] == 1,
          result$coverage$status[3] == "no_exact_target_match")
empty <- rank_enriched_tf_qscores(go[FALSE, ], matches, regulators, 10)
stopifnot(nrow(empty$top) == 0, nrow(empty$coverage) == 0)
undefined <- regulators
undefined$qscore <- NA_real_
missing <- rank_enriched_tf_qscores(go, matches, undefined, 10)
stopifnot(nrow(missing$top) == 0, missing$coverage$status[1] == "all_scores_undefined")
stopifnot(inherits(try(rank_enriched_tf_qscores(go, matches, regulators, 0), silent = TRUE), "try-error"))
pdf <- tempfile(fileext = ".pdf")
for (example in list(result, empty, missing)) {
    plot_enriched_tf_qscores(example, pdf, "test tissue", 2)
    stopifnot(file.info(pdf)$size > 0)
    pdf_lines <- readLines(pdf, warn = FALSE, skipNul = TRUE)
    stopifnot(sum(grepl("/Type /Page\\b", pdf_lines, perl = TRUE)) == 1L)
}
unlink(pdf)
cat("Passed enrichment filtering, exact mapping, per-resolution ranking, ties, missing scores and PDF checks.\n")
