"""Derive the GO universe from all targets in a headerless PANDA edge list."""
import csv
from pathlib import Path


def target_genes(handle, delimiter=","):
    genes = set()
    for line_number, row in enumerate(csv.reader(handle, delimiter=delimiter), 1):
        if not row or all(not field.strip() for field in row):
            continue
        if len(row) < 2:
            raise ValueError(f"Network row {line_number} has no target column")
        gene = row[1].strip()
        if not gene or gene.upper() in {"NA", "NAN", "NULL"}:
            raise ValueError(f"Network row {line_number} has a missing target gene")
        genes.add(gene)
    if not genes:
        raise ValueError("Cannot derive a GO background from an empty network")
    return sorted(genes)


def write_background(genes, handle):
    # The enrichment reader expects symbols in column 2, with no header.
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerows(enumerate(genes, start=1))


def main(snakemake):
    with open(snakemake.input.network, encoding="utf-8-sig", newline="") as handle:
        genes = target_genes(handle, snakemake.params.delimiter)
    output = Path(snakemake.output.background)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        write_background(genes, handle)
    print(f"GO background: {len(genes)} unique targets from {snakemake.input.network}")


if __name__ == "__main__":
    main(snakemake)
