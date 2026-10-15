"""Match selected HiDeF gene sets to complete CONDOR memberships, then plot.

Scores retain their original resolution and regulator context. A match means
identical TARGET membership, not a newly calculated score for the hierarchy.
CSV processing uses the standard library so matching can be tested separately
from matplotlib. Two passes avoid holding the full qscore table in memory.
"""

import csv
import math
from collections import defaultdict
from pathlib import Path


def read_gmt(path):
    communities = {}
    with open(path, encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3 or not fields[0]:
                raise ValueError(f"Invalid GMT record at line {line_number}")
            if fields[0] in communities:
                raise ValueError(f"Duplicate selected community: {fields[0]}")
            genes = frozenset(gene.strip() for gene in fields[2:] if gene.strip())
            if not genes:
                raise ValueError(f"Empty selected community: {fields[0]}")
            communities[fields[0]] = genes
    return communities


def score_rows(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"tar", "community", "qscore", "resolution"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Qscore CSV must contain {sorted(required)}")
        for row in reader:
            gene = row["tar"]
            community = row["community"]
            resolution = float(row["resolution"])
            if not gene or not community or not math.isfinite(resolution) or resolution <= 0:
                raise ValueError("Qscores require a gene, a community, and a positive finite resolution")
            text = row["qscore"].strip()
            value = math.nan if text.lower() in {"", "na", "nan"} else float(text)
            if math.isinf(value):
                raise ValueError(f"Infinite qscore for {gene} at resolution {resolution}")
            yield gene, community, resolution, value


def exact_matches(communities, qscore_file):
    by_gene = defaultdict(set)
    for selected, genes in communities.items():
        for gene in genes:
            by_gene[gene].add(selected)
    # A candidate is discarded as soon as a gene rules out every selected set.
    candidates, memberships, rejected = {}, {}, set()
    for gene, community, resolution, _ in score_rows(qscore_file):
        key = (resolution, community)
        if key in rejected:
            continue
        possible = by_gene.get(gene, set())
        possible = candidates[key] & possible if key in candidates else set(possible)
        if not possible:
            candidates.pop(key, None)
            memberships.pop(key, None)
            rejected.add(key)
            continue
        candidates[key] = possible
        members = memberships.setdefault(key, set())
        if gene in members:
            raise ValueError(f"Duplicate qscore for {gene} in community {community} at {resolution}")
        members.add(gene)
    return {
        key: sorted(selected for selected in options if memberships[key] == communities[selected])
        for key, options in candidates.items()
        if any(memberships[key] == communities[selected] for selected in options)
    }


def write_mapped_scores(communities, qscore_file, matches, output_file, coverage_file):
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    Path(coverage_file).parent.mkdir(parents=True, exist_ok=True)
    distributions = defaultdict(lambda: defaultdict(list))
    with open(output_file, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["selected_community", "condor_community", "resolution", "tar", "qscore"])
        for gene, community, resolution, value in score_rows(qscore_file):
            for selected in matches.get((resolution, community), []):
                writer.writerow([selected, community, resolution, gene,
                                 value if math.isfinite(value) else "NaN"])
                distributions[selected][(resolution, community)].append(value)
    with open(coverage_file, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["selected_community", "n_genes", "n_exact_matches",
                         "n_score_rows", "n_undefined", "status"])
        for selected, genes in communities.items():
            groups = distributions.get(selected, {})
            values = [value for group in groups.values() for value in group]
            writer.writerow([selected, len(genes), len(groups), len(values),
                             sum(not math.isfinite(value) for value in values),
                             "matched" if groups else "no_exact_match"])
    return distributions


def plot_distributions(communities, distributions, output_file, tissue):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(output_file) as pdf:
        fig, ax = plt.subplots(figsize=(11, 7))
        ax.axis("off")
        n_matched = sum(bool(distributions.get(name)) for name in communities)
        ax.text(.05, .9, f"Selected community qscores — {tissue}", fontsize=18, weight="bold")
        ax.text(.05, .72,
                f"Selected communities: {len(communities)}\n"
                f"With exact target-membership matches: {n_matched}\n"
                f"Without matches in the exported resolutions: {len(communities) - n_matched}",
                fontsize=13, linespacing=1.8)
        ax.text(.05, .38,
                "Only identical gene sets are matched; numeric community IDs are not compared.\n"
                "Scores retain their original CONDOR regulator assignments and resolution.\n"
                "Histograms pool gene–resolution observations, so genes can appear repeatedly.\n"
                "NaN scores are counted but excluded from numeric distributions.\n"
                "The current bihidef export omits the initial minimum resolution.\n"
                "The coverage CSV includes every selected community, including unmatched sets.",
                fontsize=11, linespacing=1.6)
        pdf.savefig(fig)
        plt.close(fig)
        for selected, genes in communities.items():
            groups = distributions.get(selected, {})
            if not groups:
                continue  # Explicitly recorded in coverage CSV and PDF overview.
            all_values = [v for values in groups.values() for v in values]
            finite = [v for v in all_values if math.isfinite(v)]
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            fig.suptitle(f"{tissue} | {selected} | {len(genes)} genes", fontsize=14)
            if finite:
                axes[0].hist(finite, bins=min(40, max(1, math.ceil(math.sqrt(len(finite))))),
                             color="#287D8E", edgecolor="white")
            else:
                axes[0].text(.5, .5, "All matched scores are undefined", ha="center",
                             transform=axes[0].transAxes)
            axes[0].set(xlabel="Qscore", ylabel="Gene–resolution observations",
                        title=f"{len(groups)} exact matches; {len(all_values) - len(finite)} NaN scores")
            points = []
            for (resolution, _), values in sorted(groups.items()):
                values = [v for v in values if math.isfinite(v)]
                if values:
                    low, median, high = np.quantile(values, [.25, .5, .75])
                    points.append((resolution, low, median, high))
            if points:
                x, low, median, high = np.array(points).T
                axes[1].errorbar(x, median, yerr=[median - low, high - median],
                                 fmt="o", color="#287D8E", capsize=3)
                axes[1].set_xscale("log")
            else:
                axes[1].text(.5, .5, "No finite scores", ha="center", transform=axes[1].transAxes)
            axes[1].set(xlabel="CONDOR resolution", ylabel="Qscore",
                        title="Median and interquartile range at each match")
            for ax in axes:
                ax.spines[["top", "right"]].set_visible(False)
            fig.tight_layout(rect=(0, 0, 1, .94))
            pdf.savefig(fig)
            plt.close(fig)


def main(snakemake):
    communities = read_gmt(snakemake.input.communities)
    matches = exact_matches(communities, snakemake.input.qscores)
    distributions = write_mapped_scores(
        communities, snakemake.input.qscores, matches,
        snakemake.output.matched_scores, snakemake.output.coverage,
    )
    plot_distributions(communities, distributions, snakemake.output.plot,
                       snakemake.wildcards.tissue_type)


if __name__ == "__main__":
    main(snakemake)
