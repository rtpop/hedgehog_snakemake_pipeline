"""Score saved Leiden memberships and selected HEDGEHOG sets without detection."""
from pathlib import Path
from time import perf_counter
import importlib.util

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix


def projection_helpers():
    """Compatibility helpers for the existing optional prior-only benchmark.

    That benchmark imports this source module directly, so __file__ points to
    src/new rather than Snakemake's temporary copy of this scoring entry point.
    """
    path = Path(__file__).with_name("leiden_gene_projection.py")
    spec = importlib.util.spec_from_file_location("projection", path)
    projection = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(projection)
    return projection


def load_network(path, delimiter):
    return projection_helpers().load_network(path, delimiter)


def gene_projection(frame):
    return projection_helpers().gene_projection(frame)


def load_communities(path, genes):
    """Read final selected GMT sets without flattening overlapping communities."""
    communities = {}
    universe = set(genes)
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 3 or not fields[0] or any(not g for g in fields[2:]):
            raise ValueError(f"Invalid GMT record at line {number}")
        name, members = fields[0], set(fields[2:])
        if name in communities:
            raise ValueError(f"Duplicate GMT community: {name}")
        unknown = members - universe
        if unknown:
            raise ValueError(f"Community {name} has genes outside the network: {sorted(unknown)[:10]}")
        communities[name] = members
    return communities


def community_scores(graph, communities, method, resolution):
    # Each set's contribution to standard weighted modularity (gamma=1).
    # Overlapping contributions must not be summed into partition modularity.
    total = sum(graph.es["weight"])
    # Standard modularity is undefined for signed weights. CPM supports them.
    modularity_defined = total > 0 and all(w >= 0 for w in graph.es["weight"])
    degree = dict(zip(graph.vs["name"], graph.strength(weights="weight")))
    adjacency = graph.get_adjacency_sparse(attribute="weight")
    indices = {gene: i for i, gene in enumerate(graph.vs["name"])}
    rows = []
    for name, members in communities.items():
        selected = [indices[g] for g in members]
        internal = float(adjacency[selected][:, selected].sum()) / 2
        contribution = (internal / total - (sum(degree[g] for g in members) / (2 * total)) ** 2
                        if modularity_defined else np.nan)
        rows.append(dict(method=method, community=name, n_genes=len(members),
                         modularity_contribution=contribution,
                         cpm_contribution=2 * internal - resolution * len(members) * (len(members) - 1)))
    return rows


def whole_partition_modularity(graph, communities):
    """Score a fractional cover at gamma=1, including uncovered singletons.

    A gene in k communities receives membership 1/k in each. Uncovered genes
    each receive their own singleton. For disjoint exhaustive communities this
    is exactly standard weighted partition modularity.
    """
    weights = graph.es["weight"]
    total = sum(weights)
    if total <= 0 or any(weight < 0 for weight in weights):
        return np.nan

    genes = graph.vs["name"]
    indices = {gene: index for index, gene in enumerate(genes)}
    memberships = [[] for _ in genes]
    for column, members in enumerate(communities.values()):
        for gene in members:
            memberships[indices[gene]].append(column)

    # Every row sums to one. Add one distinct column for each uncovered gene;
    # these singletons are used only for whole-partition scoring, not GMT output.
    row_indices, column_indices, values = [], [], []
    next_column = len(communities)
    for row, columns in enumerate(memberships):
        if not columns:
            columns = [next_column]
            next_column += 1
        row_indices.extend([row] * len(columns))
        column_indices.extend(columns)
        values.extend([1 / len(columns)] * len(columns))
    membership = coo_matrix(
        (values, (row_indices, column_indices)),
        shape=(len(genes), next_column),
    ).tocsr()

    # Q = sum_c [alpha_c^T A alpha_c / (2W)
    #            - (sum_i strength_i * alpha_ic / (2W))^2].
    # Sparse operations avoid constructing a dense gene-by-gene modularity matrix.
    adjacency = graph.get_adjacency_sparse(attribute="weight")
    observed = float(membership.multiply(adjacency @ membership).sum()) / (2 * total)
    strengths = np.asarray(graph.strength(weights="weight"), dtype=float)
    community_strengths = np.asarray(membership.T @ strengths).ravel()
    expected = np.square(community_strengths / (2 * total)).sum()
    return float(observed - expected)


def best_matches(source, target, method):
    rows = []
    for name, members in source.items():
        best_name, best_score, intersection = None, -1, 0
        for other, candidates in target.items():
            shared = len(members & candidates)
            score = shared / len(members | candidates)
            if score > best_score:
                best_name, best_score, intersection = other, score, shared
        rows.append(dict(method=method, community=name, matched_community=best_name,
                         jaccard=max(0, best_score), intersection_genes=intersection))
    return rows


def load_leiden_membership(path, genes):
    """Require exactly one saved community label per projected gene."""
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if not {"gene", "community"}.issubset(frame.columns):
        raise ValueError("Leiden membership must contain gene and community columns")
    if frame.gene.duplicated().any() or (frame[["gene", "community"]] == "").any().any():
        raise ValueError("Leiden membership contains duplicate genes or empty identifiers")
    if set(frame.gene) != set(genes):
        raise ValueError("Saved Leiden memberships do not match the filtered network genes")
    members = dict(zip(frame.gene, frame.community))
    labels = {label: index for index, label in enumerate(dict.fromkeys(frame.community))}
    partition = [labels[members[gene]] for gene in genes]
    communities = {}
    for gene in genes:
        communities.setdefault(members[gene], set()).add(gene)
    return communities, partition


def compare(network, communities, membership_file, detection_file, projection_script,
            output, delimiter, resolution, tissue, filtering):
    output = Path(output)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing scores: {output}")
    spec = importlib.util.spec_from_file_location("projection", projection_script)
    projection = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(projection)
    frame = projection.load_network(network, delimiter)
    start = perf_counter()
    graph = projection.gene_projection(frame)
    scoring_projection_seconds = perf_counter() - start
    selected = load_communities(communities, graph.vs["name"])
    leiden, partition = load_leiden_membership(membership_file, graph.vs["name"])

    # Both the original combined summary and the new detection-only summary
    # contain a Leiden row. Read its timing without rerunning find_partition.
    detection = pd.read_csv(detection_file, sep="\t")
    detection = detection.loc[detection.method == "leiden"]
    if len(detection) != 1:
        raise ValueError("Expected one Leiden detection summary row")
    metadata = detection.iloc[0]
    if (str(metadata.tissue) != tissue or str(metadata.filtering) != filtering
            or not np.isclose(float(metadata.leiden_resolution), resolution)):
        raise ValueError("Saved Leiden detection metadata does not match the scoring job")
    leiden_seconds = metadata.get("algorithm_seconds", np.nan)
    projection_seconds = metadata.get("projection_seconds", np.nan)
    scores = (community_scores(graph, selected, "hedgehog", resolution)
              + community_scores(graph, leiden, "leiden", resolution))
    modularity_defined = (sum(graph.es["weight"]) > 0
                          and all(w >= 0 for w in graph.es["weight"]))
    rows = []
    for method, sets in (("hedgehog", selected), ("leiden", leiden)):
        covered = set().union(*sets.values()) if sets else set()
        counts = [sum(g in members for members in sets.values()) for g in covered]
        rows.append(dict(tissue=tissue, filtering=filtering, method=method,
                         leiden_resolution=resolution,
                         modularity=(graph.modularity(partition, weights="weight")
                                     if method == "leiden" and modularity_defined else np.nan),
                         modularity_defined=modularity_defined,
                         whole_partition_modularity=whole_partition_modularity(graph, sets),
                         whole_partition_scoring="fractional_membership_uncovered_singletons_gamma1",
                         cpm_quality=(sum(row["cpm_contribution"] for row in scores if row["method"] == "leiden")
                                      if method == "leiden" else np.nan),
                         algorithm_seconds=leiden_seconds if method == "leiden" else np.nan,
                         projection_seconds=projection_seconds,
                         scoring_projection_seconds=scoring_projection_seconds, n_genes=graph.vcount(),
                         n_projection_edges=graph.ecount(), n_communities=len(sets),
                         n_covered_genes=len(covered), n_uncovered_genes=graph.vcount()-len(covered),
                         n_overlapping_genes=sum(c > 1 for c in counts)))
    output.mkdir(parents=True, exist_ok=False)
    for method, sets in (("hedgehog", selected), ("leiden", leiden)):
        pd.DataFrame([(g, name) for name, members in sets.items() for g in sorted(members)],
                     columns=["gene", "community"]).to_csv(
            output / f"{method}_gene_membership.tsv", sep="\t", index=False, mode="x")
    pd.DataFrame(scores, columns=["method", "community", "n_genes", "modularity_contribution", "cpm_contribution"]).to_csv(
        output / "community_scores.tsv", sep="\t", index=False, mode="x")
    pd.DataFrame(best_matches(selected, leiden, "hedgehog") + best_matches(leiden, selected, "leiden"),
                 columns=["method", "community", "matched_community", "jaccard", "intersection_genes"]).to_csv(
        output / "community_matches.tsv", sep="\t", index=False, mode="x")
    pd.DataFrame(rows).to_csv(output / "summary.tsv", sep="\t", index=False, mode="x")


if __name__ == "__main__":
    compare(snakemake.input.network, snakemake.input.communities,
            snakemake.input.membership, snakemake.input.detection_summary,
            snakemake.input.projection_script, snakemake.output.results,
            snakemake.params.delimiter, float(snakemake.params.resolution),
            snakemake.params.tissue, snakemake.params.filtering)
