# import libraries
import pandas as pd
import networkx as nx
import hedgehog

prior_file = snakemake.input.prior
panda_file = snakemake.input.panda
panda_filtered = snakemake.input.panda_filtered
filtering_bench = snakemake.output.filtering_bench
delimiter = snakemake.params.delimiter
resolution = float(snakemake.params.resolution)
max_communities = int(snakemake.params.max_communities)


def main():
    
    # Load the unfiltered PANDA network
    print("Loading PANDA network")
    panda = pd.read_csv(panda_file, delimiter=delimiter)
    
    # Load filtered networks
    print("Loading filtered networks")
    prior_fil = pd.read_csv(panda_filtered[0], delimiter=delimiter)
    hedgehog_fil = pd.read_csv(panda_filtered[1], delimiter=delimiter)

    # Calculate modularity for all networks
    modularity_hedgehog = hedgehog.filter_panda.calculate_modularity(hedgehog_fil, resolution=resolution, comm_mult=max_communities)
    modularity_prior = hedgehog.filter_panda.calculate_modularity(prior_fil, resolution=resolution, comm_mult=max_communities)
    modularity_unfiltered = hedgehog.filter_panda.calculate_modularity(panda, resolution=resolution, comm_mult=max_communities)

    # Convert dataframes to networkx graphs
    hedgehog_graph = nx.from_pandas_edgelist(hedgehog_fil, source=hedgehog_fil.columns[0], target=hedgehog_fil.columns[1], edge_attr=hedgehog_fil.columns[2])
    prior_graph = nx.from_pandas_edgelist(prior_fil, source=prior_fil.columns[0], target=prior_fil.columns[1], edge_attr=prior_fil.columns[2])
    unfiltered_graph = nx.from_pandas_edgelist(panda, source=panda.columns[0], target=panda.columns[1], edge_attr=panda.columns[2])
    
    # Calculate density and number of edges for all networks
    density_hedgehog = nx.density(hedgehog_graph)
    density_prior = nx.density(prior_graph)
    density_unfiltered = nx.density(unfiltered_graph)

    num_edges_hedgehog = hedgehog_graph.number_of_edges()
    num_edges_prior = prior_graph.number_of_edges()
    num_edges_unfiltered = unfiltered_graph.number_of_edges()
    
    # Calculate the number of unique TFs and genes for all networks
    unique_tfs_hedgehog = hedgehog_fil.iloc[:, 0].nunique()
    unique_genes_hedgehog = hedgehog_fil.iloc[:, 1].nunique()

    unique_tfs_prior = prior_fil.iloc[:, 0].nunique()
    unique_genes_prior = prior_fil.iloc[:, 1].nunique()
    
    unique_tfs_unfiltered = panda.iloc[:, 0].nunique()
    unique_genes_unfiltered = panda.iloc[:, 1].nunique()
    
    # Create a DataFrame to store the results
    results = pd.DataFrame({
        'Network': ['HEDGEHOG filtered PANDA', 'Prior filtered', 'Unfiltered'],
        'Modularity': [modularity_hedgehog, modularity_prior, modularity_unfiltered],
        'Density': [density_hedgehog, density_prior, density_unfiltered],
        'Number of Edges': [num_edges_hedgehog, num_edges_prior, num_edges_unfiltered],
        'Unique TFs': [unique_tfs_hedgehog, unique_tfs_prior, unique_tfs_unfiltered],
        'Unique Genes': [unique_genes_hedgehog, unique_genes_prior, unique_genes_unfiltered]
    })
    
    # Save the results to a CSV file
    results.to_csv(filtering_bench, index=False)

if __name__ == '__main__':
    main()
