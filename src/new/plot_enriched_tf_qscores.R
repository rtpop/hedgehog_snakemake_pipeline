source(snakemake@params$function_file)
read_csv <- function(path) utils::read.csv(path, colClasses = "character", check.names = FALSE,
                                         na.strings = character(), fileEncoding = "UTF-8-BOM")
go <- utils::read.delim(snakemake@input$go_summary, colClasses = "character", check.names = FALSE)
result <- rank_enriched_tf_qscores(go, read_csv(snakemake@input$matched_targets),
                                  read_csv(snakemake@input$regulator_scores), snakemake@params$top_n)
terms <- read_enriched_go_terms(go, snakemake@input$go_details, snakemake@params$sig_thresh)
for (path in c(snakemake@output$plot, snakemake@output$top_tfs, snakemake@output$coverage,
              snakemake@output$go_terms)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
utils::write.table(result$top, snakemake@output$top_tfs, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
utils::write.table(result$coverage, snakemake@output$coverage, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
utils::write.table(terms, snakemake@output$go_terms, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
plot_enriched_tf_qscores(result, snakemake@output$plot, snakemake@wildcards$tissue_type,
                         snakemake@params$top_n, terms)
