# Community analyses

The active Snakefile uses `script:` entry points in `src/new`. The original
scripts under `src/old` remain available to the archived workflow.

## GO enrichment

`go_enrichment` reads the selected GMT and writes the configured `go_results`
summary under `<output_dir>/<tissue>/go_enrichment/`. Per-community term tables
are tracked in its `community_terms` directory when `save_all` is true.
`plot_n_enriched_communities` writes `go_enrichment_counts.pdf`.

The defaults restored from `config_old.yaml` use human gene symbols,
biological-process terms, and topGO classic/Fisher testing. With `auto_bg: false`,
the `derive_go_background` rule reads ALL unique target symbols (column 2) from
the tissue's original, unfiltered PANDA edge list. It writes
`<output_dir>/<tissue>/go_enrichment/network_background.tsv`, with a row number
and gene symbol, without a header. This includes targets outside selected
communities and targets removed by edge filtering. No external annotation file
is required. Missing target identifiers cause an explicit error. topGO performs
the human BP annotation mapping using this supplied universe.
With `auto_bg: true`, the union of the selected GMT genes supplies
the background and no background file is required. These backgrounds answer
different enrichment questions; the configuration choice is deliberate.

BH correction now uses the numeric p-values for **all** tested GO terms within
each community, rather than the formatted default top-ten `GenTable` output.
Consequently, adjusted p-values and significant-term counts can differ from
the old script. The configured significance threshold uses a strict `<` test.
An empty GMT produces a header-only summary and an explanatory plot. Unannotated
communities are retained with zero significant terms, missing p-values, and a
status column. The detailed result tables retain raw and adjusted p-values.

## Qscore distributions

`run_bihidef` now declares `qscore_tar.csv` and `qscore_reg.csv` as outputs.
The `plot_selected_community_qscores` rule runs entirely in R and consumes the target CSV and the
selected GMT. Under the corresponding `bihidef_run` directory it writes:

- `selected_community_qscores.pdf`: an overview and one distribution page per
  matched selected community; each page includes the median and interquartile
  range by resolution.
- `selected_community_qscores.csv`: matched gene-level scores with both
  `selected_community` and `condor_community` IDs, plus the original resolution.
- `selected_community_qscore_coverage.csv`: every selected community, its match
  count, score-row count, undefined-score count, and match status.

Only identical complete target gene sets match. A matching numeric ID, subset,
or overlap is insufficient. These are existing CONDOR scores evaluated with the
original regulator memberships and resolution, not newly calculated scores for
the final consensus hierarchy. Final communities may have no exact match.
The CSV exporter currently omits the initial minimum resolution, so coverage
means coverage of the exported resolutions. No approximate matches are made.
NaN scores are retained in the mapped CSV and counted in the coverage report;
they are excluded from numeric plots. Histograms pool gene-resolution records,
so a gene can contribute once per matching resolution. Regulator scores are not
mixed into the selected gene-community plots.

The R container already lists topGO, org.Hs.eg.db, data.table and ggplot2, and
the qscore matching and plotting scripts use base R. No new runtime packages or container
rebuild are required for these workflow scripts if those images match their
Dockerfiles. Adding the CSV outputs may make Snakemake rerun existing bihidef jobs
if either CSV is missing; do not create placeholder score files.

Validate on the execution host with `snakemake --cores 1 --dry-run`, then run
the desired targets using the usual container flags. Matching regression tests:
`Rscript tests/test_selected_qscores.R`.

## Community consolidation and counts

`rule all` now includes `consolidate_community_stats` and `plot_n_communities`.
The consolidation rule consumes the selected GMT, selection statistics, and GO
summary for every configured tissue. Their existing producer rules supply any
missing inputs automatically. Tissue labels come directly from `tissue` in
`config.yaml`, so nested or custom output directories are supported.

The following config values are paths relative to `output_dir`:

- `community_stats_consolidated`: `community_stats_consolidated.tsv`, one row per
  selected community with its unique gene count, significant GO term count and
  tissue-level selection statistics.
- `community_tissue_summary`: `community_tissue_summary.tsv`, one row per tissue
  including tissues with zero selected communities. `n_enriched` counts selected
  communities with at least one significant GO term, using the existing GO output.
- `community_plot`: `community_counts.pdf`, total, selected and GO-enriched counts
  across tissues, generated from the tissue summary.

Run the normal workflow or target `plot_n_communities` / `consolidate_community_stats`
directly with the usual Snakemake container flags. The Python script uses only the
standard library, and the plot uses base R. The rules use the existing
`python_container` and `r_container`; no additional packages or images are needed
for these steps. Archived `Snakefile_old` and `config_old.yaml` are unchanged.

Consolidation fails on missing or duplicate GO community IDs, incomplete input
lists, or inconsistent selected counts rather than silently dropping records.
Regression checks: `python -m unittest discover -s tests -p test_community_consolidation.py -v`.
