# Import libraries
import argparse
import pandas as pd
import eland

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process and filter panda network.")

    # Parse args for the main function
    parser.add_argument('input_file', type=str, help='Path to edge list file.')
    parser.add_argument('prior_file', type=str, help='Path to prior file.')
    parser.add_argument('output_file', type=str, help='Path to output file.')
    parser.add_argument('delimiter', type=str, help='Delimiter used in the edge list files')
    
    return parser.parse_args()

def main():
    args = parse_arguments()
    
        # Define new file names with "_updated_sep" added
    updated_prior_file = args.prior_file.replace(".txt", "_updated_sep.txt")
    updated_panda_file = args.input_file.replace(".txt", "_updated_sep.txt")
    
    # Process prior and panda files so they have the same separator
    prior = pd.read_csv(args.prior_file, sep="\t", header=None)
    prior.to_csv(updated_prior_file, sep=args.delimiter, index=False, header=False)
    
    panda = pd.read_csv(args.input_file, sep=" ", header=1)
    panda.to_csv(updated_panda_file, sep=args.delimiter, index=False, header=False)

    # process panda result into edgelist
    eland.process_edge_list(args.input_file, args.output_file, sep = args.delimiter)

if __name__ == "__main__":
    main()
