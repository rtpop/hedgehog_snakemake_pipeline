"""Shared network reader and weighted gene projection for the Leiden benchmark."""
import igraph as ig
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


