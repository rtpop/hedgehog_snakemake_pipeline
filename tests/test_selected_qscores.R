# Run from the repository root: Rscript tests/test_selected_qscores.R
source("src/new/selected_qscores_fn.R")
sets <- list(final_A = c("g1", "g2"), final_B = "g4")
scores <- data.frame(tar = c("g1", "g2", "g1", "g2", "g1", "g2", "g3", "g1"),
    community = c("007", "007", "9", "9", "007", "007", "007", "8"),
    qscore = c(.3, .2, .4, NA, .1, .1, .3, .5), resolution = c(1, 1, 2, 2, 3, 3, 3, 4))
result <- match_selected_qscores(sets, scores)
stopifnot(nrow(result$scores) == 4,
          identical(unique(result$scores$selected_community), "final_A"),
          setequal(result$scores$condor_community, c("007", "9")),
          result$coverage$n_exact_matches[1] == 2,
          result$coverage$n_undefined[1] == 1,
          result$coverage$status[2] == "no_exact_match")
duplicate <- try(match_selected_qscores(sets, rbind(scores, scores[1, ])), silent = TRUE)
stopifnot(inherits(duplicate, "try-error"))
empty <- match_selected_qscores(list(), scores)
stopifnot(nrow(empty$scores) == 0, nrow(empty$coverage) == 0,
          "selected_community" %in% names(empty$coverage))
aliases <- match_selected_qscores(list(A = c("g1", "g2"), B = c("g2", "g1")), scores)
stopifnot(nrow(aliases$scores) == 8)
undefined <- match_selected_qscores(list(C = "g1"),
    data.frame(tar = "g1", community = "0", resolution = 1, qscore = "NaN"))
stopifnot(undefined$coverage$n_undefined == 1)
tmp <- tempfile(fileext = ".pdf")
for (example in list(result, empty, undefined)) {
    plot_selected_qscores(example, tmp, "test tissue")
    stopifnot(file.info(tmp)$size > 0)
}
unlink(tmp)
cat("Passed exact matching, ID preservation, NaN, empty selection, duplicate detection, and PDF rendering checks.\n")
