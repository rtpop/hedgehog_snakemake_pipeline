# import libraries
import argparse
import pandas as pd
import eland

def parse_arguments():
    parser = argparse.ArgumentParser(description="benchmark filtering methods")

    # Parse args for the main function
    parser.add_argument('prior_file', type=str, help='Path to edge list file.')
    parser.add_argument('input_file', type=str, help='Path to edge list file.')
    parser.add_argument('filtered-net', type=str, help='Path to filtered network file.')
    parser.add_argument('output_file', type=str, help='Path to output file.')
    parser.add_argument('--delimiter', type=str, help='Delimiter used in the edge list files') 
    
    return parser.parse_args()

def main():
    args = parse_arguments()

    # load filtered network
    eland_fil = pd.read_csv(args.filtered-net, delimiter = args.delimiter)
    
    # load panda network as df
    panda = pd.read_csv(args.input_file, delimiter = args.delimiter)
    
    # filter to prior edges
    prior_fil = eland.filter_panda.filter_panda(args.prior_file, args.input_file, args.output_file, delimiter = args.delimiter, prior_only = True)
    
    # filter based on top n edges
    # n is the number of edges of eland_fil
    top_fil = panda.nlargest(len(eland_fil), 2)
    
    # calculate modularity for the two networks
    modularity_eland = eland.filter_panda.calculate_modularity(eland_fil)
    modularity_prior = eland.filter_panda.calculate_modularity(prior_fil)
    modularity_top = eland.filter_panda.calculate_modularity(top_fil)
    
    # save modularity values
    with open(args.output_file, 'a') as f:
        f.write(f"Modularity of the ELAND filtered PANDA network: {modularity_eland}\n")
        f.write(f"Modularity of the prior filtered network: {modularity_prior}\n")
        f.write(f"Modularity of the top filtered network: {modularity_top}\n")