"""Compare HEDGEHOG's BRIM partition with Leiden CPM on a gene projection."""
from pathlib import Path
from time import perf_counter

import igraph as ig
import leidenalg as la
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from hedgehog.filter_panda import condor


def load_network(path, delimiter):
    # filter_panda writes headerless edge lists. Preserve identifiers literally.
    frame = pd.read_csv(path, sep=delimiter, header=None,
                        names=["source", "target", "weight"],
                        dtype={"source": str, "target": str}, keep_default_na=False)
    frame["weight"] = pd.to_numeric(frame["weight"], errors="raise")
    if frame.empty or (frame[["source", "target"]] == "").any().any():
        raise ValueError("Network must contain nonempty TF and gene identifiers")
    if not np.isfinite(frame.weight).all() or (frame.weight < 0).any():
        raise ValueError("This weighted benchmark requires finite, nonnegative weights")
    if frame.weight.sum() <= 0:
        raise ValueError("Network must have positive total weight")
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


def hedgehog_membership(frame, resolution, comm_mult):
    # Same calls as hedgehog.filter_panda.calculate_modularity, retaining tar_memb.
    model = condor.condor_object(dataframe=frame.copy(), silent=True)
    model.initial_community(resolution=resolution)
    maximum = int(model.tar_memb.community.nunique() * comm_mult)
    model.brim(resolution=resolution, c=maximum)
    # CONDOR adds exactly one tar_ prefix to target identifiers.
    members = {str(name).removeprefix("tar_"): community
               for name, community in zip(model.tar_memb.tar, model.tar_memb.community)}
    if set(members) != set(frame.target):
        raise ValueError("HEDGEHOG target memberships do not match network genes")
    return members


def compare(network, output, delimiter, resolution, comm_mult, tissue, filtering):
    output = Path(output)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing comparison: {output}")
    frame = load_network(network, delimiter)
    start = perf_counter()
    graph = gene_projection(frame)
    projection_seconds = perf_counter() - start
    start = perf_counter()
    members = hedgehog_membership(frame, resolution, comm_mult)
    hedgehog_seconds = perf_counter() - start
    # Compact arbitrary community labels for igraph.
    labels = {label: i for i, label in enumerate(dict.fromkeys(members.values()))}
    hedgehog_partition = [labels[members[gene]] for gene in graph.vs["name"]]
    start = perf_counter()
    leiden_partition = la.find_partition(
        graph, la.CPMVertexPartition, resolution_parameter=resolution,
        weights="weight", seed=0).membership
    leiden_seconds = perf_counter() - start
    rows = []
    for method, partition, seconds in (
        ("hedgehog", hedgehog_partition, hedgehog_seconds),
        ("leiden", leiden_partition, leiden_seconds),
    ):
        rows.append(dict(tissue=tissue, filtering=filtering, method=method,
                         resolution=resolution,
                         modularity=graph.modularity(partition, weights="weight"),
                         algorithm_seconds=seconds, projection_seconds=projection_seconds,
                         n_genes=graph.vcount(), n_projection_edges=graph.ecount(),
                         n_communities=len(set(partition))))
    output.mkdir(parents=True, exist_ok=False)
    for method, partition in (("hedgehog", hedgehog_partition), ("leiden", leiden_partition)):
        pd.DataFrame({"gene": graph.vs["name"], "community": partition}).to_csv(
            output / f"{method}_gene_membership.tsv", sep="\t", index=False, mode="x")
    pd.DataFrame(rows).to_csv(output / "summary.tsv", sep="\t", index=False, mode="x")


if __name__ == "__main__":
    compare(snakemake.input.network, snakemake.output.results,
            snakemake.params.delimiter, float(snakemake.params.resolution),
            int(snakemake.params.comm_mult), snakemake.params.tissue,
            snakemake.params.filtering)
