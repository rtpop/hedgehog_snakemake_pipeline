from hedgehog import process_bihidef

communities = snakemake.input.communities
selected_communities = snakemake.output.selected_communities
max_genes = snakemake.params.max_genes
min_genes = snakemake.params.min_genes
stats = snakemake.output.stats

def main():
    
    # Select communities
    comms = process_bihidef.select_communities(communities, min_genes, max_genes, stats)
    
    # Save selected communities to a gmt file    
    process_bihidef.gmt_from_bihidef(comms, selected_communities)
    
if __name__ == "__main__":
    main()
