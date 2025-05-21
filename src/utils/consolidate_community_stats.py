import sys
import os
import pandas as pd

def parse_stats_file(stats_file):
    """Parse a stats file and return a dictionary of stats."""
    stats = {}
    with open(stats_file) as f:
        lines = f.readlines()
        
    # Skip the first line (date/time)
    for line in lines[1:]:
        if ':' in line:
            key, value = line.strip().split(':', 1)
            stats[key.strip()] = value.strip()      
    return stats

def parse_gmt_file(gmt_file):
    """Parse a GMT file and return a list of (community_name, size) tuples."""
    comms = []
    with open(gmt_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                comm_name = parts[0]
                size = len(parts) - 2
                comms.append((comm_name, size))
    return comms

def infer_tissue_from_path(path):
    """Infer tissue name from file path: o"""
    parts = path.split(os.sep)
    try:
        idx = parts.index('output')
        tissue = parts[idx + 1]
    except Exception:
        tissue = "unknown"
    return tissue

def main():
    stats_files = sys.argv[1].split(",")
    gmt_files = sys.argv[2].split(",")
    out_file = sys.argv[3]
    
    all_records = []
    for stats_file, gmt_file in zip(stats_files, gmt_files):
        tissue = infer_tissue_from_path(stats_file)
        stats = parse_stats_file(stats_file)
        comms = parse_gmt_file(gmt_file)
        for comm_name, size in comms:
            record = {
                "tissue": tissue,
                "community": comm_name,
                "size": size
            }
            record.update(stats)
            all_records.append(record)
    
    df = pd.DataFrame(all_records)
    
    # Rename columns
    df = df.rename(columns={
        "Maximum number of genes per community": "max_genes",
        "Minimum number of genes per community": "min_genes",
        "Number of selected communities": "n_selected",
        "Total number of communities": "n_total"
    })
    # Reorder columns
    df = df[["tissue", "community", "siize", "max_genes", "min_genes", "n_selected", "n_total"]]
    df.to_csv(out_file, sep="\t", index=False)

if __name__ == "__main__":
    main()