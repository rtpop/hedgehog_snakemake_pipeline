"""Leiden on the prior-filtered projection when the main comparison uses another filter."""
import importlib.util
from pathlib import Path

import leidenalg as la
import numpy as np
import pandas as pd

# Reuse exactly the network reader and weighted projection of the main benchmark.
spec = importlib.util.spec_from_file_location("comparison", snakemake.input.comparison_script)
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)

output = Path(snakemake.output.results)
if output.exists():
    raise FileExistsError(f"Refusing to overwrite existing comparison: {output}")
frame = comparison.load_network(snakemake.input.network, snakemake.params.delimiter)
graph = comparison.gene_projection(frame)
resolution = float(snakemake.params.resolution)
partition = la.find_partition(
    graph, la.CPMVertexPartition, resolution_parameter=resolution,
    weights="weight", seed=0,
)
defined = sum(graph.es["weight"]) > 0 and all(weight >= 0 for weight in graph.es["weight"])
row = dict(
    tissue=snakemake.params.tissue, filtering="prior", method="leiden",
    leiden_resolution=resolution,
    modularity=graph.modularity(partition.membership, weights="weight") if defined else np.nan,
)
output.mkdir(parents=True, exist_ok=False)
pd.DataFrame([row]).to_csv(output / "summary.tsv", sep="\t", index=False, mode="x")
