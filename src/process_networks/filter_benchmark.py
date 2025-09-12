# import libraries
import argparse
import pandas as pd
import networkx as nx
import hedgehog

def parse_arguments():
    parser = argparse.ArgumentParser(description="benchmark filtering methods")

    # Parse args for the main function
    parser.add_argument('--prior_file', type=str, help='Path to edge list file.')
    parser.add_argument('--panda', type=str, help='Path to edge list file.')
    parser.add_argument('--filtered_net', type=str, nargs="+", help='Path to filtered network file.')
    parser.add_argument('--output_file', type=str, help='Path to output file.')
    parser.add_argument('--delimiter', type=str, help='Delimiter used in the edge list files')
    parser.add_argument('--resolution', type=float, help='Resolution for modularity calculation')
    parser.add_argument('--max_communities', type=float, help='Maximum number of communities')

    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # Load the unfiltered PANDA network
    print("Loading PANDA network")
    panda = pd.read_csv(args.panda, delimiter=args.delimiter)
    
    # Load filtered networks
    print("Loading filtered networks")
    prior_fil = pd.read_csv(args.filtered_net[0], delimiter=args.delimiter)
    hedgehog_fil = pd.read_csv(args.filtered_net[1], delimiter=args.delimiter)

    # Calculate modularity for all networks
    modularity_hedgehog = hedgehog.filter_panda.calculate_modularity(hedgehog_fil, resolution=args.resolution, comm_mult=args.max_communities)
    modularity_prior = hedgehog.filter_panda.calculate_modularity(prior_fil, resolution=args.resolution, comm_mult=args.max_communities)
    modularity_unfiltered = hedgehog.filter_panda.calculate_modularity(panda, resolution=args.resolution, comm_mult=args.max_communities)

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
        'Network': ['ELAND filtered PANDA', 'Prior filtered', 'Unfiltered'],
        'Modularity': [modularity_hedgehog, modularity_prior, modularity_unfiltered],
        'Density': [density_hedgehog, density_prior, density_unfiltered],
        'Number of Edges': [num_edges_hedgehog, num_edges_prior, num_edges_unfiltered],
        'Unique TFs': [unique_tfs_hedgehog, unique_tfs_prior, unique_tfs_unfiltered],
        'Unique Genes': [unique_genes_hedgehog, unique_genes_prior, unique_genes_unfiltered]
    })
    
    # Save the results to a CSV file
    results.to_csv(args.output_file, index=False)

if __name__ == '__main__':
    main()
