"""Concatenate tissue summaries without changing their scores or row order."""
import csv
from pathlib import Path


columns = ["tissue", "method", "resolution", "modularity"]
output = Path(snakemake.output.summary)
output.parent.mkdir(parents=True, exist_ok=True)

with output.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    for path in snakemake.input.summaries:
        with open(path, newline="") as source:
            reader = csv.DictReader(source, delimiter="\t")
            if reader.fieldnames != columns:
                raise ValueError(f"Unexpected Leiden summary columns in {path}")
            writer.writerows(reader)
