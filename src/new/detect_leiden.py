"""Run Leiden detection once and save memberships and timing for later scoring."""
import importlib.util
from pathlib import Path
from time import perf_counter

import leidenalg as la
import pandas as pd

spec = importlib.util.spec_from_file_location("projection", snakemake.params.projection_script)
projection = importlib.util.module_from_spec(spec)
spec.loader.exec_module(projection)

frame = projection.load_network(snakemake.input.network, snakemake.params.delimiter)
start = perf_counter()
graph = projection.gene_projection(frame)
projection_seconds = perf_counter() - start
resolution = float(snakemake.params.resolution)
start = perf_counter()
result = la.find_partition(
    graph, la.CPMVertexPartition, resolution_parameter=resolution,
    weights="weight", seed=0,
)
algorithm_seconds = perf_counter() - start

membership_file = Path(snakemake.output.membership)
summary_file = Path(snakemake.output.detection_summary)
membership_file.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame({"gene": graph.vs["name"], "community": result.membership}).to_csv(
    membership_file, sep="\t", index=False,
)
pd.DataFrame([dict(
    tissue=snakemake.params.tissue, filtering=snakemake.params.filtering,
    method="leiden", leiden_resolution=resolution, algorithm_seconds=algorithm_seconds,
    projection_seconds=projection_seconds, cpm_quality=result.quality(),
    n_genes=graph.vcount(), n_projection_edges=graph.ecount(),
)]).to_csv(summary_file, sep="\t", index=False)
