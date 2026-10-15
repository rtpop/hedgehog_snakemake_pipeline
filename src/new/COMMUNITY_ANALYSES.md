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
`<data_dir>/<bg_file>` must contain background symbols in its second column
(no header). With `auto_bg: true`, the union of the selected GMT genes supplies
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
The new `plot_selected_community_qscores` rule consumes the target CSV and the
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
the Python container lists matplotlib/numpy. No new runtime packages or container
rebuild are required for these workflow scripts if those images match their
Dockerfiles. Adding the CSV outputs may make Snakemake rerun existing bihidef jobs
if either CSV is missing; do not create placeholder score files.

Validate on the execution host with `snakemake --cores 1 --dry-run`, then run
the desired targets using the usual container flags. Matching regression tests:
`python -m unittest discover -s tests -p 'test_selected_qscores.py'`.
