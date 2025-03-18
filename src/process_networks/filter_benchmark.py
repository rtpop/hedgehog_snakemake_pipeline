# import libraries
import argparse
import pandas as pd
import eland
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description="benchmark filtering methods")

    # Parse args for the main function
    parser.add_argument('--prior_file', type=str, help='Path to edge list file.')
    parser.add_argument('--panda_edgelist', type=str, help='Path to edge list file.')
    parser.add_argument('--filtered_net', type=str, help='Path to filtered network file.')
    parser.add_argument('--output_file', type=str, help='Path to output file.')
    parser.add_argument('--delimiter', type=str, help='Delimiter used in the edge list files') 
    
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # load filtered network
    print("Loading filtered network")
    eland_fil = pd.read_csv(args.filtered_net, delimiter=args.delimiter)
    
    # load panda network as df
    print("Loading PANDA network")
    panda = pd.read_csv(args.panda_edgelist, delimiter=args.delimiter)
    
    # filter to prior edges
    prior_fil = eland.filter_panda.filter_panda(args.prior_file, args.panda_edgelist, delimiter=args.delimiter, prior_only=True)
    
    # filter based on top n edges
    # n is the number of edges of eland_fil
    top_fil = panda.iloc[:, 2].nlargest(len(eland_fil)).index
    top_fil = panda.loc[top_fil]
        
    # calculate modularity for the two networks
    modularity_eland = eland.filter_panda.calculate_modularity(eland_fil)
    modularity_prior = eland.filter_panda.calculate_modularity(prior_fil)
    modularity_top = eland.filter_panda.calculate_modularity(top_fil)
    
    # save modularity values
    with open(args.output_file, 'a') as f:
        f.write(f"Modularity of the ELAND filtered PANDA network: {modularity_eland}\n")
        f.write(f"Modularity of the prior filtered network: {modularity_prior}\n")
        f.write(f"Modularity of the top filtered network: {modularity_top}\n")

if __name__ == '__main__':
    main()
