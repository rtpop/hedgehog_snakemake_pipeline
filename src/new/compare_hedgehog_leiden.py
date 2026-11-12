"""Compare selected HEDGEHOG gene sets with Leiden on a gene projection."""
from pathlib import Path
from time import perf_counter

import igraph as ig
import leidenalg as la
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix


def load_network(path, delimiter):
    # filter_panda writes headerless edge lists. Preserve identifiers literally.
    frame = pd.read_csv(path, sep=delimiter, header=None,
                        names=["source", "target", "weight"],
                        dtype={"source": str, "target": str}, keep_default_na=False)
    frame["weight"] = pd.to_numeric(frame["weight"], errors="raise")
    if frame.empty or (frame[["source", "target"]] == "").any().any():
        raise ValueError("Network must contain nonempty TF and gene identifiers")
    if not np.isfinite(frame.weight).all():
        raise ValueError("Network weights must be finite")
    return frame


def gene_projection(frame):
    genes = sorted(frame.target.unique())
    tfs = sorted(frame.source.unique())
    matrix = coo_matrix((frame.weight.to_numpy(), (
        frame.source.map({name: i for i, name in enumerate(tfs)}),
        frame.target.map({name: i for i, name in enumerate(genes)}))),
        shape=(len(tfs), len(genes))).tocsr()
    projection = (matrix.T @ matrix).tocsr()
    projection.setdiag(0)
    projection.eliminate_zeros()
    edges = projection.tocoo()
    upper = edges.row < edges.col
    graph = ig.Graph(n=len(genes), edges=list(zip(
        edges.row[upper].tolist(), edges.col[upper].tolist())), directed=False)
    graph.vs["name"] = genes
    graph.es["weight"] = edges.data[upper].tolist()
    if not graph.ecount():
        raise ValueError("Gene projection has no edges; modularity is undefined")
    return graph


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


def compare(network, communities, output, delimiter, resolution, tissue, filtering):
    output = Path(output)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing comparison: {output}")
    frame = load_network(network, delimiter)
    start = perf_counter()
    graph = gene_projection(frame)
    projection_seconds = perf_counter() - start
    selected = load_communities(communities, graph.vs["name"])
    start = perf_counter()
    result = la.find_partition(
        graph, la.CPMVertexPartition, resolution_parameter=resolution,
        weights="weight", seed=0)
    partition = result.membership
    leiden_seconds = perf_counter() - start
    leiden = {}
    for gene, label in zip(graph.vs["name"], partition):
        leiden.setdefault(str(label), set()).add(gene)
    scores = (community_scores(graph, selected, "hedgehog", resolution)
              + community_scores(graph, leiden, "leiden", resolution))
    modularity_defined = all(w >= 0 for w in graph.es["weight"])
    rows = []
    for method, sets in (("hedgehog", selected), ("leiden", leiden)):
        covered = set().union(*sets.values()) if sets else set()
        counts = [sum(g in members for members in sets.values()) for g in covered]
        rows.append(dict(tissue=tissue, filtering=filtering, method=method,
                         leiden_resolution=resolution,
                         modularity=(graph.modularity(partition, weights="weight")
                                     if method == "leiden" and modularity_defined else np.nan),
                         modularity_defined=modularity_defined,
                         cpm_quality=result.quality() if method == "leiden" else np.nan,
                         algorithm_seconds=leiden_seconds if method == "leiden" else np.nan,
                         projection_seconds=projection_seconds, n_genes=graph.vcount(),
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
    compare(snakemake.input.network, snakemake.input.communities, snakemake.output.results,
            snakemake.params.delimiter, float(snakemake.params.resolution),
            snakemake.params.tissue,
            snakemake.params.filtering)
