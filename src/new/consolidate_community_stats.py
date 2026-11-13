"""Consolidate selected community and GO statistics with explicit tissue labels."""
import csv
from pathlib import Path

STAT_KEYS = {
    "Maximum number of genes per community": "max_genes",
    "Minimum number of genes per community": "min_genes",
    "Selected communities": "n_selected",
    "Total communities": "n_total",
}
SUMMARY_COLUMNS = ["tissue", "max_genes", "min_genes", "n_selected", "n_enriched", "n_total"]
COMMUNITY_COLUMNS = ["tissue", "community", "size", "n_sig_terms"] + SUMMARY_COLUMNS[1:]


def parse_stats_file(path):
    stats = {}
    for line in Path(path).read_text().splitlines():
        key, separator, value = line.partition(":")
        if separator and key.strip() in STAT_KEYS:
            number = int(value.strip())
            if number < 0:
                raise ValueError(f"Negative statistic in {path}: {key}")
            stats[STAT_KEYS[key.strip()]] = number
    if set(stats) != set(STAT_KEYS.values()):
        raise ValueError(f"Missing community statistics in {path}")
    return stats


def parse_gmt_file(path):
    communities = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 3 or not fields[0] or fields[0] in communities:
            raise ValueError(f"Invalid or duplicate GMT community in {path}")
        genes = {gene.strip() for gene in fields[2:] if gene.strip()}
        if not genes:
            raise ValueError(f"Empty GMT community in {path}")
        communities[fields[0]] = len(genes)
    return communities


def read_go_summary(path):
    result = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"community", "n_sig_terms"}.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing GO summary columns in {path}")
        for row in reader:
            name = row["community"]
            if not name or name in result:
                raise ValueError(f"Missing or duplicate GO community in {path}")
            count = int(row["n_sig_terms"])
            if count < 0:
                raise ValueError(f"Negative significant term count in {path}")
            result[name] = count
    return result


def write_table(path, rows, columns):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def consolidate(tissues, stats_files, gmt_files, go_files, output, summary_output):
    if len(set(tissues)) != len(tissues):
        raise ValueError("Tissue labels must be unique")
    if not (len(tissues) == len(stats_files) == len(gmt_files) == len(go_files)):
        raise ValueError("Expected one statistics, GMT and GO file per tissue")
    records, summaries = [], []
    for tissue, stats_file, gmt_file, go_file in zip(tissues, stats_files, gmt_files, go_files):
        stats = parse_stats_file(stats_file)
        communities = parse_gmt_file(gmt_file)
        go = read_go_summary(go_file)
        if stats["n_selected"] != len(communities):
            raise ValueError(f"Selected count differs from GMT for {tissue}")
        if stats["n_total"] < stats["n_selected"]:
            raise ValueError(f"Total community count is smaller than selected count for {tissue}")
        if set(go) != set(communities):
            raise ValueError(f"GO and selected GMT community IDs differ for {tissue}")
        summary = dict(tissue=tissue, **stats, n_enriched=sum(n > 0 for n in go.values()))
        summaries.append(summary)
        for name, size in communities.items():
            records.append(dict(summary, community=name, size=size, n_sig_terms=go[name]))
    write_table(output, records, COMMUNITY_COLUMNS)
    write_table(summary_output, summaries, SUMMARY_COLUMNS)


if __name__ == "__main__":
    consolidate(list(snakemake.params.tissues), list(snakemake.input.community_stats),
                list(snakemake.input.selected_communities), list(snakemake.input.go_summary),
                snakemake.output.consolidated_community_stats, snakemake.output.tissue_summary)
